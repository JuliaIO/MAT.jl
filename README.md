# MAT.jl
[![Build Status](https://github.com/JuliaIO/MAT.jl/workflows/CI/badge.svg?branch=master)](https://github.com/JuliaIO/MAT.jl/actions)

[![docs-stable][docs-stable-img]][docs-stable-url]
[![docs-dev][docs-dev-img]][docs-dev-url]

### Read and write MATLAB files in Julia

This library can read MATLAB `.mat` files, both in the older v4/v5/v6/v7 format, as well as the newer v7.3 format.

## Installation

This is installed using the standard tools of the [package manager](https://julialang.github.io/Pkg.jl/v1/getting-started/):

```julia
pkg> add MAT
```
where you get the `pkg>` prompt by hitting `]` as the first character of the line. (Exit `pkg` mode by hitting backspace or Ctrl-C as the first character of the line.)

See also the requirements for the [HDF5](https://github.com/timholy/HDF5.jl/) module, used for "v7.3" files and for writing \*.mat files.

## Usage

To load the module:

```julia
using MAT
```

To read a single variable from a MAT file (compressed files are detected and handled automatically):

```julia
file = matopen("matfile.mat")
read(file, "varname") # note that this does NOT introduce a variable ``varname`` into scope
close(file)
```

To write a variable to a MAT file:

```julia
file = matopen("matfile.mat", "w")
write(file, "varname", variable)
close(file)
```

To read all variables from a MAT file as a Dict:

```julia
vars = matread("matfile.mat")
```

To write a Dict to a MAT file, using its keys as variable names.
The `compress` argument is optional, and compression is off by default:

```julia
matwrite("matfile.mat", Dict(
	"myvar1" => 0,
	"myvar2" => 1
); compress = true)
```

To write in MATLAB v4 format:

```julia
matwrite("matfile.mat", Dict(
	"myvar1" => 0,
	"myvar2" => 1
);version="v4")
```

To get a list of variable names in a MAT file:

```julia
file = matopen("matfile.mat")
varnames = keys(file)
close(file)
```

To check for the presence of a variable name in a MAT file:

```julia
file = matopen("matfile.mat")
if haskey(file, "variable")
    # something
end
close(file)
```

## Cell versus struct array writing

Cell arrays are written for Any arrays `Array{Any}` or any other unsupported element type:

```julia
sarr = Any[
    Dict("x"=>1.0, "y"=>2.0),
    Dict("x"=>3.0, "y"=>4.0)
]
matwrite("matfile.mat", Dict("cell" => sarr))

```

Inside MATLAB you will find:

```matlab
>> load('matfile.mat')
>> cell

cell =

  2×1 cell array

    {1×1 struct}
    {1×1 struct}
```

Read and write behavior for struct arrays is different. For struct arrays (in .mat v7.3 files) we use the `MAT.MatlabStructArray` type. You can also write with MAT.jl using Dict arrays `AbstractArray{<:AbstractDict}` if all the Dicts have equal keys, which will automatically convert internally to `MAT.MatlabStructArray`.

```julia
sarr = Dict{String, Any}[
    Dict("x"=>1.0, "y"=>2.0),
    Dict("x"=>3.0, "y"=>4.0)
]
matwrite("matfile.mat", Dict("struct_array" => sarr))

# which is the same as:
matwrite("matfile.mat", Dict("struct_array" => MAT.MatlabStructArray(sarr)))
```

Now you'll find the following inside MATLAB:

```matlab
>> load('matfile.mat')
>> struct_array

struct_array =

[2x1 struct, 576 bytes]
x: 1
y: 2
```

Note that when you read the file again, you'll find the `MAT.MatlabStructArray`, which you can convert back to the Dict array with `Array`:

```julia
julia> sarr = matread("matfile.mat")["struct_array"]
MAT.MAT_types.MatlabStructArray{1}
 "x": Any[1.0, 3.0]
 "y": Any[2.0, 4.0]

julia> sarr["x"]
2-element Vector{Any}:
 1.0
 3.0

julia> Array(sarr)
2-element Vector{Dict{String, Any}}:
 Dict("x" => 1.0, "y" => 2.0)
 Dict("x" => 3.0, "y" => 4.0)

```

Note that in MAT.jl v0.10 and older, or .mat versions before v7.3, will read struct arrays as a Dict with concatenated arrays in the fields/keys, which is equal to `Dict(sarr)`.

## Caveats

* All files are written in MATLAB v7.3 format by default.
* Writing in MATLAB v4 format is provided by the matwrite function with keyword argument.

## Credits

The MAT_HDF5 module, which provides read/write support for MATLAB v7.3 files, was written primarily by [Tim Holy](https://github.com/timholy/). The MAT_v5 module, which provides read support for MATLAB v5/v6/v7 files, was written primarily by [Simon Kornblith](https://github.com/simonster/). The MAT_v4 module, which provides read and write support for MATLAB v4 files, was written primarily by [Victor Saase](https://github.com/vsaase/).


<!-- URLs -->
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaIO.github.io/MAT.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaIO.github.io/MAT.jl/dev
