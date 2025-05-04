The `traceStreamline` function computes streamlines of vector fields and n-vector fields on meshes.

`#include "geometrycentral/surface/streamlines.h"`

??? func "`#!cpp std::vector<SurfacePoint> traceStreamline(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, SurfacePoint pStart, const FaceData<Vector2>& field, size_t nSym = 1, TraceStreamlineOptions opt = TraceStreamlineOptions());`"

    Trace a streamline (path following a vector field) on a surface mesh.

    - `mesh`: The input mesh
    - `geom`: the input geometry (as always, a `VertexPositionGeometry` is valid input)
    - `pStart`: Starting point on the surface
    - `field`: Vector field defined per face. (When working with n-vector fields, each face should contain one representative vector)
    - `nSym`: Symmetry order of the field (1 for ordinary vector fields, >1 for n-vector fields)
    - `opt`: Options controlling the tracing behavior

    Returns a list of points along the streamline path.

!!! warning
    When `nSym > 1`, the field should contain a representative vector (one of the nSym vectors in the face), not the "power vector" used by the [direction field methods](/surface/algorithms/direction_fields). A power vector can be converted into a representative vector by calling `powerVector.pow(1. / nSym)`.

??? func "`#!cpp std::vector<std::vector<SurfacePoint>> traceManyStreamlines(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, const FaceData<Vector2>& field, size_t nSym = 1, TraceStreamlineOptions opt = TraceStreamlineOptions(), PoissonDiskOptions pdOpt = PoissonDiskOptions());`"

    Trace multiple streamlines starting from points sampled via Poisson disk sampling.

    - `mesh`, `geom`, `field`, `nSym`, `opt`: Same as `traceStreamline`
    - `pdOpt`: Options for Poisson disk sampling to determine starting points

    Returns a collection of streamline paths.

This convenience function automatically samples starting points across the mesh and traces streamlines from each. Note that the input must now be a `VertexPositionGeometry`, as the vertex positions are used when performing Poisson disk sampling.

## Helper Types

### TraceStreamlineOptions

Options are passed in to the streamline tracing routines via a `TraceStreamlineOptions` object:

| Field | Default | Description |
|-------|---------|-------------|
| `#!cpp size_t maxSegments` | `10` | Maximum number of segments in the streamline |
| `#!cpp double maxLen` | `∞` | Maximum length of the streamline along the mesh |
| `#!cpp FaceData<int>* nVisits` | `nullptr` | If set, counts face visits during tracing |
| `#!cpp const FaceData<int>* maxVisits` | `nullptr` | If set along with `nVisits`, limits visits per face |
| `#!cpp int verbosity` | `0` | Debug output level (>0 for debug info) |

## Example

```cpp
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/streamlines.h"
#include "geometrycentral/surface/meshio.h"

// Load mesh
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;
std::unique_ptr<CornerData<Vector2>> parameterization;
std::tie(mesh, geom, param) = readParameterizedManifoldSurfaceMesh(filename);

// Create a vector field
int nSym = 1;
FaceData<Vector2> vectorField = computeSmoothestFaceDirectionField(*geom, nSym);

// Set some options
TraceStreamlineOptions options;
options.maxSegments = 100;
options.maxLen = 1.0;

// Trace streamlines
auto streamlines = traceManyStreamlines(*mesh, *geom, vectorField, 1, options);

// Save streamlines to an svg texture
SvgCurveOptions svgOptions;
svgOptions.curveColor = "#0000ff";
draw_mesh_curves_to_svg(*mesh, *geom, *parameterization, streamlines,
                       "streamlines.svg", svgOptions);
```

## Visualization

### `draw_mesh_curves_to_svg`

??? func "`#!cpp void draw_mesh_curves_to_svg(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const CornerData<Vector2>& param, const std::vector<std::vector<SurfacePoint>>& curves, std::string filepath, SvgCurveOptions opt = SvgCurveOptions());`"

    Export streamlines to an SVG file using a UV parameterization.

    - `mesh`, `geom`: The surface mesh and its geometry
    - `param`: UV coordinates for each corner
    - `curves`: Collection of curves to draw (e.g., output from `traceManyStreamlines`)
    - `filepath`: Output SVG file path
    - `opt`: Rendering options

### SvgCurveOptions

Controls SVG rendering appearance:

| Field | Default | Description |
|-------|---------|-------------|
| `#!cpp double worldspaceStrokeWidth` | `0.05` | Stroke width relative to mesh scale |
| `#!cpp std::string backgroundColor` | `"#fff"` | SVG background color |
| `#!cpp std::string curveColor` | `"#000000"` | Default curve color |
| `#!cpp std::function<std::string(size_t, size_t)> curveColorFunction` | `nullptr` | Custom color function (overrides `curveColor`) |
| `#!cpp double imageSize` | `500` | SVG viewBox size |
