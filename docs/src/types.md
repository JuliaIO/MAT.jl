# Types and conversions

MAT.jl uses the following type conversions from MATLAB types to Julia types:

| MATLAB    | Julia |
| -------- | ------- |
| numerical array | `Array{T}` |
| cell array | `Array{Any}` |
| char array | `String` |
| `struct`  | `Dict{String,Any}`    |
| `struct` array | `MAT.MatlabStructArray`     |
| old class object    | `MAT.MatlabClassObject`    |
| new (opaque) class    | `MAT.MatlabOpaque`    |

A few of the `MatlabOpaque` classes are automatically converted upon reading:

| MATLAB    | Julia |
| -------- | ------- |
| `string`    | `String`    |
| `datetime`    | `Dates.DateTime`    |
| `duration`    | `Dates.Millisecond`    |
| `category`    | `PooledArrays.PooledArray`    |
| `table`    | `MAT.MatlabTable` (or any other table) |
| `timetable`    | `MAT.MatlabTable` (or any other table), with the row times as its first column |

A timetable's row times come back as `Dates.DateTime` or `Dates.Millisecond`, in a
column named after the row dimension (`Time` unless renamed). A regular timetable, which
MATLAB stores as a start time and a sample rate or time step, has its row times
generated; one stepped in calendar units (months, for example) has no fixed rate and is
left as a `MatlabOpaque`, with a warning. Row times are rounded to whole milliseconds,
the resolution of `Dates.DateTime`.

Note that single element arrays are typically converted to scalars in Julia, because MATLAB cannot distinguish between scalars and `1x1` sized arrays.