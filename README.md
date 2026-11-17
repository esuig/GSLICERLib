# GSlicerLib

A single-header C++ library that slices a triangle mesh with a set of parallel
planes and reconstructs, for each plane, the closed profiles (contours) of the
cross-section. 

> **Project status:** this library is at version 0.0.1, an early release. It
> is still being tested and will be further developed, so its API may change
> in future versions. Feedback and bug reports are welcome.
>
> This library has been tested only on Windows and compiled with MSVC 2022 so far
>
> A test suite for the library will be introduced soon

## Features

- **Slicing along any direction**: the slicing planes are defined by a start
  point and a normal, so the mesh can be sliced along any axis, not only Z.
- **Closed, oriented profiles**: the intersection segments are joined into
  closed profiles; outer contours come out counterclockwise and holes
  clockwise (seen from the tip of the plane normal), so the two can be told
  apart without further processing.
- **Robust handling of degenerate cases**: vertices and edges lying exactly on
  a slicing plane, triangles lying in a plane, and non-manifold crossings
  (several contours touching at one point) are handled explicitly.
- **Efficient**: each triangle is tested only against the planes its bounding
  box can reach, and segments are joined through shared mesh edges instead of
  by geometric search.
- **Plain data in and out**: the mesh is passed as flat arrays of coordinates
  and indices, and the profiles are returned as flat coordinate arrays, so the
  library works with any mesh or geometry types.

## Getting started

Copy `GSlicerLib.h` into your project and include it:

```cpp
#include "GSlicerLib.h"
```

Requires C++17. Everything is in the `gslicerlib` namespace.

## Input

- **Vertices**: a flat vector of coordinates, `x, y, z` for each vertex.
- **Triangles**: a flat vector of vertex indices, three per triangle.
- **Orientation**: triangles must be consistently oriented, counterclockwise
  when seen from outside the mesh (the usual convention, e.g. in STL files).
  The orientation of the profiles is derived from it.
- **Watertight mesh**: for closed profiles, the mesh should be closed. Holes or
  gaps in the mesh produce profiles that cannot be closed, and slicing fails
  with error code 2.
- **Slicing planes**: the first plane passes through the given start point and
  is perpendicular to the given normal (which does not need to be unit
  length). The following planes are placed every `layer height` along the
  normal. Only the part of the mesh in front of the first plane (on the side
  the normal points to) is sliced.

## Output

The profiles are returned in a three-level container:

| Level         | Content                                                                                 |
| ------------- | --------------------------------------------------------------------------------------- |
| Outer vector  | One entry per slicing plane that produced profiles, in slicing order.                   |
| Middle vector | The profiles of that plane.                                                             |
| Inner vector  | One profile: flat coordinates `x, y, z, ...`, closed (first point repeated at the end). |

Planes that produce no profiles (e.g. planes before the mesh starts) are
omitted, so the position of a layer in the outer vector is not the plane
index: read the height of a layer from the coordinates of its profiles.

TODO: add an empty vector for planes not containing anything

Outer contours are oriented counterclockwise and holes clockwise, when seen
from the tip of the plane normal looking back along it.

## Return values

| Value | Meaning                                                                                   |
| ----- | ----------------------------------------------------------------------------------------- |
| `0`   | Slicing completed (the output is empty if the mesh lies entirely behind the first plane). |
| `1`   | Invalid input: fewer than 3 coordinates or indices, or layer height not positive.         |
| `2`   | A profile could not be closed, typically because the mesh is not watertight.              |

## License

Released under the [Boost Software License 1.0](LICENSE).

SPDX: `BSL-1.0`

## References

- P. J. Schneider, D. H. Eberly, *Geometric Tools for Computer Graphics*,
  Morgan Kaufmann, 2003: segment-plane intersection (p. 484) and
  triangle-plane intersection configurations (section 11.5).
