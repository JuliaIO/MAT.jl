# Types and conversions

MAT.jl uses the following type conversions from MATLAB types to Julia types when reading:

| MATLAB    | Julia |
| -------- | ------- |
| numeric array (`double`, `single`, `int8` to `uint64`) | `Array{T}` (`Float64`, `Float32`, `Int8` to `UInt64`) |
| complex numeric array | `Array{Complex{T}}` |
| `logical` array | `Array{Bool}` |
| sparse matrix (numeric, logical or complex) | `SparseArrays.SparseMatrixCSC` |
| char array with one row | `String` |
| char array with several rows | `Vector{String}`, one string per row |
| cell array | `Array{Any}` |
| `struct`  | `Dict{String,Any}`    |
| `struct` array | `MAT.MatlabStructArray`     |
| function handle | `MAT.FunctionHandle`, which holds MATLAB's description of the handle; it cannot be called |
| old class object    | `MAT.MatlabClassObject`    |
| new (opaque) class    | `MAT.MatlabOpaque`    |
| array of new class objects | `Array{MAT.MatlabOpaque}` |

Note that single element arrays are typically converted to scalars in Julia, because MATLAB cannot distinguish between scalars and `1x1` sized arrays. v4 files are the exception: there a scalar reads as a `1x1` matrix.

A few more details:

- Empty arrays keep their element type and size: `zeros(0, 3)` reads as a 0×3 `Matrix{Float64}`. An empty char array `''` reads as `String[]`.
- A single character reads as a `Char` from a v7.3 file and as a one-character `String` from a v4, v6 or v7 file. This applies inside cell arrays and structs too.
- The rows of a char array keep the trailing spaces MATLAB pads them with.
- A handle object that a file refers to more than once reads as one `MatlabOpaque`, so the references stay identical (`===`).

A few of the `MatlabOpaque` classes are automatically converted upon reading:

| MATLAB    | Julia |
| -------- | ------- |
| `string`    | `String`, or an `Array{String}` for a string array    |
| `datetime`    | `Dates.DateTime` (`TimeZones.ZonedDateTime` if zoned and TimeZones.jl is loaded)    |
| `duration`    | `Dates.Millisecond`    |
| `categorical`    | `PooledArrays.PooledArray`    |
| `table`    | `MAT.MatlabTable` (or any other table) |
| `timetable`    | `MAT.MatlabTable` (or any other table), with the row times as its first column |

- A `datetime` keeps its instant; how a time zone is read is described under [Time zones](@ref) below. The display format is not kept.
- The categories of a `categorical` become the pool of the `PooledArray`, in MATLAB's order. Whether it is ordinal is not kept.
- A `table` keeps its column names and data. Its row names and properties, such as units and descriptions, are not kept. Pass `table = DataFrame` (or another Tables.jl type) to `matread` or `matopen` to read tables, and timetables, as that type.
- Pass `convert_opaque = false` to `matread` or `matopen` to read these classes as `MatlabOpaque` instead.

A timetable's row times come back as `Dates.DateTime` or `Dates.Millisecond`, in a
column named after the row dimension (`Time` unless renamed). A regular timetable, which
MATLAB stores as a start time and a sample rate or time step, has its row times
generated; one stepped in calendar units (months, for example) has no fixed rate and is
left as a `MatlabOpaque`, with a warning. Row times are rounded to whole milliseconds,
the resolution of `Dates.DateTime`.

Other MATLAB classes, such as `calendarDuration`, `containers.Map` and `dictionary`, and your own classes, read as a `MatlabOpaque` holding the properties MATLAB saved. Enumeration values are not supported: reading one returns the saved metadata, with a warning.

## Time zones

MATLAB stores a `datetime` that has a time zone as its UTC instant together with the
zone's name. With [TimeZones.jl](https://github.com/JuliaTime/TimeZones.jl) loaded (Julia
1.9 and later, through a package extension), it is read as a `ZonedDateTime` in that zone:
an IANA zone such as `Europe/London`, `UTC`, or a fixed offset such as `+05:30`.

```julia
using MAT, TimeZones
matread("file.mat")["t"]   # 2022-07-20T12:00:00+01:00, in Europe/London
```

Without it, a zoned `datetime` is read as a `DateTime` holding the UTC instant, with a
warning. A `datetime` without a time zone is always a `DateTime` of the stored wall-clock
time. Zoned row times of a timetable follow the same rules. The `UTCLeapSeconds` zone counts
leap seconds, which neither type can represent, so such a `datetime` is left as a
`MatlabOpaque`, with a warning.
