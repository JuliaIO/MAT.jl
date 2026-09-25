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
| `datetime`    | `Dates.DateTime` (`TimeZones.ZonedDateTime` if zoned and TimeZones.jl is loaded)    |
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

Note that single element arrays are typically converted to scalars in Julia, because MATLAB cannot distinguish between scalars and `1x1` sized arrays.