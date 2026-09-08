#!/usr/bin/env python3
"""Compare short serial and two-rank MPI solver runs and restarts."""

import argparse
import math
import os
import select
import signal
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = REPO_ROOT / "runSimulation.sh"
SMOKE_PARAMS = REPO_ROOT / "testCases" / "smoke.params"
BUILD_FILES = ("bretherton", "bretherton.c", "build.meta")
INVOCATION_TAG = "BRETHERTON_PARALLEL_SMOKE_INVOCATION"


def decode_output(output: str | bytes | None) -> str:
    """Normalize partial subprocess output for an on-disk receipt."""
    if output is None:
        return ""
    if isinstance(output, bytes):
        return output.decode("utf-8", errors="replace")
    return output


def invocation_processes(process_group: int, invocation_tag: str) -> set[int]:
    """Find non-zombie processes by process group or inherited unique tag."""
    members = set()
    expected_tag = f"{INVOCATION_TAG}={invocation_tag}".encode()
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        process_id = int(entry.name)
        try:
            status_lines = (entry / "status").read_text(encoding="utf-8").splitlines()
            process_state = next(
                line.split()[1] for line in status_lines if line.startswith("State:")
            )
            if process_state == "Z":
                continue
            in_process_group = os.getpgid(process_id) == process_group
            try:
                process_environment = (entry / "environ").read_bytes().split(b"\0")
                has_invocation_tag = expected_tag in process_environment
            except PermissionError:
                has_invocation_tag = False
            if in_process_group or has_invocation_tag:
                members.add(process_id)
        except (FileNotFoundError, ProcessLookupError, PermissionError, StopIteration):
            continue
    return members


def path_is_within(candidate: Path, root: Path) -> bool:
    """Return whether candidate is root or one of its descendants."""
    try:
        candidate.relative_to(root)
    except ValueError:
        return False
    return True


def open_verified_group(
    process_group: int, work_root: Path, invocation_tag: str
) -> tuple[list[tuple[int, int]], list[str]]:
    """Open pidfds and prove ownership of every process-group member."""
    initial_members = invocation_processes(process_group, invocation_tag)
    evidence = [f"process_group={process_group}",
                f"initial_members={sorted(initial_members)}"]
    if not initial_members:
        raise RuntimeError("no live process-group members were available for ownership proof")

    expected_uid = os.getuid()
    expected_tag = f"{INVOCATION_TAG}={invocation_tag}".encode()
    handles: list[tuple[int, int]] = []
    try:
        for process_id in sorted(initial_members):
            process_fd = os.pidfd_open(process_id)
            handles.append((process_id, process_fd))

            process_root = Path(f"/proc/{process_id}")
            process_uid = process_root.stat().st_uid
            process_cwd = Path(os.readlink(process_root / "cwd")).resolve()
            process_environment = (process_root / "environ").read_bytes().split(b"\0")
            if process_uid != expected_uid:
                raise RuntimeError(
                    f"pid {process_id} has uid {process_uid}, expected {expected_uid}"
                )
            if not path_is_within(process_cwd, work_root):
                raise RuntimeError(
                    f"pid {process_id} cwd {process_cwd} is outside {work_root}"
                )
            if expected_tag not in process_environment:
                raise RuntimeError(f"pid {process_id} lacks the invocation tag")
            evidence.append(
                f"pid={process_id} uid={process_uid} cwd={process_cwd} tag=verified"
            )

        final_members = invocation_processes(process_group, invocation_tag)
        evidence.append(f"final_members={sorted(final_members)}")
        if final_members != initial_members:
            raise RuntimeError(
                "process-group membership changed during ownership proof: "
                f"{sorted(initial_members)} -> {sorted(final_members)}"
            )
        return handles, evidence
    except Exception:
        for _, process_fd in handles:
            os.close(process_fd)
        raise


def signal_verified_group(
    handles: list[tuple[int, int]], process_group: int, invocation_tag: str
) -> list[str]:
    """Terminate only the processes held by previously verified pidfds."""
    evidence = []
    try:
        for process_id, process_fd in handles:
            try:
                signal.pidfd_send_signal(process_fd, signal.SIGTERM)
                evidence.append(f"pid={process_id} signal=TERM")
            except ProcessLookupError:
                evidence.append(f"pid={process_id} exited_before_TERM")

        pending_fds = {process_fd for _, process_fd in handles}
        term_deadline = time.monotonic() + 5.0
        while pending_fds and time.monotonic() < term_deadline:
            readable, _, _ = select.select(
                list(pending_fds), [], [], max(0.0, term_deadline - time.monotonic())
            )
            pending_fds.difference_update(readable)

        for process_id, process_fd in handles:
            if process_fd not in pending_fds:
                continue
            try:
                signal.pidfd_send_signal(process_fd, signal.SIGKILL)
                evidence.append(f"pid={process_id} signal=KILL")
            except ProcessLookupError:
                evidence.append(f"pid={process_id} exited_before_KILL")

        deadline = time.monotonic() + 5.0
        while (invocation_processes(process_group, invocation_tag) and
               time.monotonic() < deadline):
            time.sleep(0.05)
        survivors = invocation_processes(process_group, invocation_tag)
        evidence.append(f"survivors={sorted(survivors)}")
        if survivors:
            raise RuntimeError(f"verified process-group members survived: {sorted(survivors)}")
        return evidence
    finally:
        for _, process_fd in handles:
            os.close(process_fd)


def write_params(destination: Path, case_no: int, **updates: str) -> None:
    """Write a smoke parameter file with selected key replacements."""
    values = {"CaseNo": str(case_no), **{key: str(value) for key, value in updates.items()}}
    lines = []
    replaced = set()
    for line in SMOKE_PARAMS.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#") and "=" in stripped:
            key = stripped.split("=", 1)[0].strip()
            if key in values:
                lines.append(f"{key}={values[key]}")
                replaced.add(key)
                continue
        lines.append(line)
    lines.extend(f"{key}={value}" for key, value in values.items() if key not in replaced)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def invoke(params: Path, output_root: Path, extra: list[str], timeout: int,
           *, expected_status: int = 0) -> str:
    """Run the case runner and return its combined diagnostic stream."""
    command = ["bash", str(RUNNER), str(params), *extra]
    work_root = output_root.parent.resolve(strict=True)
    invocation_tag = uuid.uuid4().hex
    environment = os.environ.copy()
    environment["OUTPUT_ROOT"] = str(output_root)
    environment[INVOCATION_TAG] = invocation_tag
    stage = "build" if "--build-only" in extra else "run"
    receipt = output_root.parent / f"{output_root.name}-{params.stem}-{stage}.log"
    process = subprocess.Popen(
        command,
        cwd=work_root,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )
    try:
        output, _ = process.communicate(timeout=timeout + 30)
    except subprocess.TimeoutExpired as timeout_error:
        partial_output = decode_output(timeout_error.output)
        timeout_evidence = [
            "",
            "# Harness timeout ownership evidence",
            f"invocation_tag={invocation_tag}",
            f"launcher_pid={process.pid}",
            f"work_root={work_root}",
        ]
        try:
            handles, ownership_evidence = open_verified_group(
                process.pid, work_root, invocation_tag
            )
            timeout_evidence.extend(ownership_evidence)
        except Exception as ownership_error:
            timeout_evidence.append(f"ownership_proof=FAILED: {ownership_error}")
            receipt.write_text(
                partial_output + "\n" + "\n".join(timeout_evidence) + "\n",
                encoding="utf-8",
            )
            raise RuntimeError(
                "runner timed out; ownership proof failed, so no process was signalled; "
                f"evidence retained in {receipt}"
            ) from timeout_error

        try:
            timeout_evidence.extend(
                signal_verified_group(handles, process.pid, invocation_tag)
            )
            output, _ = process.communicate(timeout=5)
        except Exception as cleanup_error:
            timeout_evidence.append(f"verified_cleanup=FAILED: {cleanup_error}")
            receipt.write_text(
                partial_output + "\n" + "\n".join(timeout_evidence) + "\n",
                encoding="utf-8",
            )
            raise RuntimeError(
                "runner timed out; verified cleanup did not complete; "
                f"evidence retained in {receipt}"
            ) from timeout_error

        receipt.write_text(
            decode_output(output) + "\n" + "\n".join(timeout_evidence) + "\n",
            encoding="utf-8",
        )
        raise RuntimeError(
            f"runner exceeded {timeout + 30} seconds; its verified processes were "
            f"terminated and evidence was retained in {receipt}"
        ) from timeout_error

    receipt.write_text(output, encoding="utf-8")
    if process.returncode != expected_status:
        raise RuntimeError(
            f"command returned {process.returncode}, expected {expected_status}: "
            f"{' '.join(command)}\n"
            f"{output}"
        )
    if "CFL must be <=" in output:
        raise AssertionError(f"VOF CFL warning in {receipt}")
    return output


def copy_build(source_case: Path, destination_case: Path) -> None:
    """Copy a verified case-local build without copying runtime output."""
    destination_case.mkdir(parents=True, exist_ok=True)
    for name in BUILD_FILES:
        shutil.copy2(source_case / name, destination_case / name)


def last_metrics(case_dir: Path) -> list[float]:
    """Return the final numerical log row and require a single header."""
    log_files = list(case_dir.glob("c*-log"))
    if len(log_files) != 1:
        raise AssertionError(f"expected one case log in {case_dir}, found {log_files}")
    lines = log_files[0].read_text(encoding="utf-8").splitlines()
    if sum(line.startswith("# CaseNo") for line in lines) != 1:
        raise AssertionError(f"expected one log header in {log_files[0]}")
    rows = [line for line in lines if line and not line.startswith("#")]
    if not rows:
        raise AssertionError(f"no numerical rows in {log_files[0]}")
    return [float(value) for value in rows[-1].split()[:8]]


def assert_metrics_close(serial: list[float], mpi: list[float], label: str) -> None:
    """Compare state columns while allowing reduction-order roundoff."""
    for index, name in ((2, "t"), (3, "ke"), (4, "volume"),
                        (5, "front"), (6, "rear"), (7, "film")):
        if not math.isclose(serial[index], mpi[index], rel_tol=5e-4, abs_tol=1e-7):
            raise AssertionError(
                f"{label} {name} mismatch: serial={serial[index]:.9g}, "
                f"mpi={mpi[index]:.9g}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rankfile", required=True, type=Path,
                        help="Open MPI rankfile containing at least ranks 0 and 1")
    parser.add_argument("--work-root", type=Path,
                        help="retain test cases beneath this directory")
    parser.add_argument("--timeout", type=int, default=600,
                        help="timeout in seconds for each runner invocation")
    args = parser.parse_args()

    if (sys.platform != "linux" or not hasattr(os, "pidfd_open") or
            not hasattr(signal, "pidfd_send_signal")):
        parser.error("guarded timeout handling requires Linux pidfd support")

    rankfile = args.rankfile.expanduser().resolve(strict=True)
    if args.timeout < 1:
        parser.error("--timeout must be positive")

    temporary_root = None
    if args.work_root:
        work_root = args.work_root.expanduser().resolve()
        if work_root.exists() and any(work_root.iterdir()):
            parser.error("--work-root must be empty; retain old test evidence separately")
        work_root.mkdir(parents=True, exist_ok=True)
    else:
        temporary_root = Path(tempfile.mkdtemp(prefix="bretherton-parallel-smoke-"))
        work_root = temporary_root

    serial_output = work_root / "serial-cases"
    mpi_output = work_root / "mpi-cases"
    mpi_checkpoint_output = work_root / "mpi-checkpoint-cases"
    stop_output_root = work_root / "mpi-stop-cases"
    params_root = work_root / "params"
    serial_output.mkdir(parents=True, exist_ok=True)
    mpi_output.mkdir(parents=True, exist_ok=True)
    mpi_checkpoint_output.mkdir(parents=True, exist_ok=True)
    stop_output_root.mkdir(parents=True, exist_ok=True)
    params_root.mkdir(parents=True, exist_ok=True)

    serial_params = params_root / "serial.params"
    mpi_params = params_root / "mpi.params"
    write_params(serial_params, 9901)
    write_params(mpi_params, 9901)

    invoke(serial_params, serial_output, ["--build-only"], args.timeout)
    invoke(serial_params, serial_output, ["--no-build"], args.timeout)
    invoke(mpi_params, mpi_output,
           ["--ranks", "2", "--rankfile", str(rankfile), "--build-only"],
           args.timeout)
    invoke(mpi_params, mpi_output,
           ["--ranks", "2", "--rankfile", str(rankfile),
            "--mpi-timeout", str(args.timeout), "--no-build"],
           args.timeout)

    serial_fresh_metrics = last_metrics(serial_output / "9901")
    mpi_fresh_metrics = last_metrics(mpi_output / "9901")
    assert_metrics_close(serial_fresh_metrics, mpi_fresh_metrics, "fresh run")
    for case_dir in (serial_output / "9901", mpi_output / "9901"):
        if (not (case_dir / "restart").is_file() or
                not list((case_dir / "intermediate").glob("snapshot-*"))):
            raise AssertionError(f"missing restart or snapshot output in {case_dir}")

    mpi_checkpoint_params = params_root / "mpi-checkpoint-restart.params"
    write_params(mpi_checkpoint_params, 9901, tmax="0.06")
    copy_build(mpi_output / "9901", mpi_checkpoint_output / "9901")
    shutil.copy2(mpi_output / "9901" / "restart",
                 mpi_checkpoint_output / "9901" / "restart")
    shutil.copy2(mpi_output / "9901" / "c9901-log",
                 mpi_checkpoint_output / "9901" / "c9901-log")
    mpi_checkpoint_restart_output = invoke(
        mpi_checkpoint_params,
        mpi_checkpoint_output,
        ["--ranks", "2", "--rankfile", str(rankfile),
         "--mpi-timeout", str(args.timeout), "--no-build"],
        args.timeout,
    )
    if (
        "Restart file found - simulation will resume from checkpoint."
        not in mpi_checkpoint_restart_output
    ):
        raise AssertionError("MPI checkpoint restart was not detected by the runner")
    mpi_checkpoint_restart_metrics = last_metrics(
        mpi_checkpoint_output / "9901"
    )
    if mpi_checkpoint_restart_metrics[2] <= mpi_fresh_metrics[2]:
        raise AssertionError(
            "MPI checkpoint restart did not advance beyond its checkpoint time"
        )

    serial_restart_params = params_root / "serial-restart.params"
    mpi_restart_params = params_root / "mpi-restart.params"
    write_params(serial_restart_params, 9901, tmax="0.06")
    write_params(mpi_restart_params, 9901, tmax="0.06")
    shutil.copy2(serial_output / "9901" / "restart",
                 mpi_output / "9901" / "restart")
    invoke(serial_restart_params, serial_output, ["--no-build"], args.timeout)
    invoke(mpi_restart_params, mpi_output,
           ["--ranks", "2", "--rankfile", str(rankfile),
            "--mpi-timeout", str(args.timeout), "--no-build"],
           args.timeout)
    assert_metrics_close(last_metrics(serial_output / "9901"),
                         last_metrics(mpi_output / "9901"), "restart")
    assert_metrics_close(last_metrics(mpi_output / "9901"),
                         mpi_checkpoint_restart_metrics,
                         "serial-owned/MPI-owned checkpoint restart")

    stop_params = params_root / "mpi-collective-stop.params"
    write_params(stop_params, 9902, Ldomain="7")
    copy_build(mpi_output / "9901", stop_output_root / "9902")
    stop_output = invoke(
        stop_params,
        stop_output_root,
        ["--ranks", "2", "--rankfile", str(rankfile),
         "--mpi-timeout", str(args.timeout), "--no-build"],
        args.timeout,
        expected_status=1,
    )
    if "Front tip reached the outlet buffer" not in stop_output:
        raise AssertionError("MPI early-stop case did not exercise the outlet guard")
    last_metrics(stop_output_root / "9902")

    print("PASS: serial/MPI fresh run, cross-mode and MPI checkpoint restarts, "
          "outputs, and collective early stop")
    print(f"work root: {work_root}")
    if temporary_root is not None:
        shutil.rmtree(temporary_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
