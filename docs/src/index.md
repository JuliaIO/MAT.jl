```@meta
CurrentModule = MAT
```

# MAT.jl Documentation

## Overview

This Julia package
[(MAT.jl)](https://github.com/JuliaIO/MAT.jl)
provides tools
for reading and writing MATLAB format data files in Julia.

## Basic Usage

To read a single variable from a MAT file (compressed files are detected and handled automatically):

```julia
using MAT
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

## Reading and writing Julia structs

Julia structs are written as MATLAB structs, and can be read back as the original type by passing it to `matread` or `read`:

```julia
struct Person
    name::String
    age::Int
end

matwrite("people.mat", Dict("person" => Person("Alice", 30)))

matread("people.mat", "person", Person)   # Person("Alice", 30)

file = matopen("people.mat")
read(file, "person", Person)
close(file)
```

Typed reads are built on [StructUtils.jl](https://github.com/JuliaServices/StructUtils.jl), whose macros are re-exported by MAT.jl: `@defaults` for default field values, `@tags` for renaming or ignoring fields, `@choosetype` for picking a concrete type from the stored data, plus `@noarg` and `@kwarg`. Field names are matched against the MATLAB field names, and extra fields in the file are ignored.

```julia
@defaults struct Measurement
    timestamp::Float64 = 0.0 &(mat=(name="t",),)     # stored as `t` in the .mat file
    values::Vector{Float64} = Float64[]
    scratch::Float64 = 0.0 &(mat=(ignore=true,),)    # not written
end
```

An ignored field is not written, so a type that has one needs `@defaults` (as above) or `@noarg` to supply a value when the type is read back.

A `@choosetype` function is given the stored MATLAB struct as a `Dict`. The exception is a struct field of the abstract type that holds a MATLAB struct array, which reaches the function as a `MatlabStructArray` even when it has a single element. Indexing it by field name returns that field's whole column, so a function meant to handle this case should test for fields with `haskey` rather than compare their values. Struct arrays stored as the values of a MATLAB struct cannot be read into a `Dict` whose value type is such an abstract type, as in `Dict{String,MyAbstractType}`; read them into a struct with typed fields instead.

Conversion between MATLAB's and Julia's data models happens automatically: element types are converted (a MATLAB `single`, `int32` or `logical` array can be read as a `Matrix{Float64}`), trailing singleton dimensions are added or dropped to match the requested dimensionality, and a MATLAB array with a single non-singleton dimension (a row or column vector) can be read as a `Vector`. A shape that cannot be matched is reported rather than being silently flattened. An empty array reads as `nothing` or `missing` for a `Union{Nothing,T}` or `Union{Missing,T}` field, matching MATLAB's use of `[]` for a missing value, and such fields are written back as `[]`. Empty character data is the exception: it reads as `""`, because that is how MATLAB stores an empty string, and treating it as a missing value would discard something the file does hold. (MATLAB v4 files store `''` and `[]` identically, so there an empty string does read as a missing value.) A value that cannot be constructed as the requested type raises a `MATConstructionError`, whose `cause` field holds the underlying error.

Types that MATLAB has no representation for can be converted on the way out by extending `StructUtils.lower`, which is applied to every value MAT.jl writes:

```julia
struct Celsius
    degrees::Float64
end

StructUtils.lower(::MATWriteStyle, x::Celsius) = x.degrees
```

The way back is `StructUtils.lift`, which builds a value of the requested field type from what was stored. Declaring the type `structlike = false` tells StructUtils it is a single value rather than a struct to be rebuilt field by field:

```julia
StructUtils.structlike(::MATStyle, ::Type{Celsius}) = false
StructUtils.lift(::Type{Celsius}, x::Float64) = Celsius(x)
```
