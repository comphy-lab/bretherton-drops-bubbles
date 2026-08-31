/**
# getFacets.c

Extract `f` interface facets from a saved snapshot.

## Purpose

Restore a snapshot written by `simulationCases/bretherton.c` and emit
line segments from `output_facets(f, ...)` on standard output, for
interface-shape and film-thickness analysis.

## Build Example

```bash
qcc -O2 -Wall -disable-dimensions postProcess/getFacets.c -o getFacets -lm
```
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

scalar f[];
char filename[256];

/**
## main()

Usage: `./getFacets snapshot`

#### Arguments

- `snapshot`: Basilisk dump/snapshot file to restore.

#### Returns

`0` after writing facet segments to standard output.
*/
int main (int a, char const *arguments[])
{
  if (a < 2) {
    fprintf (ferr, "Usage: %s snapshot\n", arguments[0]);
    return 1;
  }

  snprintf (filename, sizeof(filename), "%s", arguments[1]);
  restore (file = filename);

  output_facets (f, stdout);
  fflush (stdout);
  return 0;
}
