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

Note that single element arrays are typically converted to scalars in Julia, because MATLAB cannot distinguish between scalars and `1x1` sized arrays.