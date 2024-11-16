# MAT.jl
# Tools for reading MATLAB v5 files in Julia
#
# Copyright (C) 2012   Timothy E. Holy and Simon Kornblith
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

module MAT

using HDF5, SparseArrays

include("MAT_HDF5.jl")
include("MAT_v5.jl")
include("MAT_v4.jl")

using .MAT_HDF5, .MAT_v5, .MAT_v4

export matopen, matread, matwrite, @read, @write

# Open a MATLAB file
const HDF5_HEADER = UInt8[0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]
function matopen(filename::AbstractString, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool, compress::Bool)
    # When creating new files, create as HDF5 by default
    fs = filesize(filename)
    if cr && (tr || fs == 0)
        return MAT_HDF5.matopen(filename, rd, wr, cr, tr, ff, compress)
    elseif fs == 0
        error("File \"$filename\" does not exist and create was not specified")
    end

    rawfid = open(filename, "r")

    # Check for MAT v4 file
    (isv4, swap_bytes) = MAT_v4.checkv4(rawfid)
    if isv4
        return MAT_v4.matopen(rawfid, swap_bytes)
    end

    # Test whether this is a MAT file
    if fs < 128
        close(rawfid)
        error("File \"$filename\" is too small to be a supported MAT file")
    end

    # Check for MAT v5 file
    seek(rawfid, 124)
    version = read(rawfid, UInt16)
    endian_indicator = read(rawfid, UInt16)
    if (version == 0x0100 && endian_indicator == 0x4D49) ||
       (version == 0x0001 && endian_indicator == 0x494D)
        if wr || cr || tr || ff
            error("creating or appending to MATLAB v5 files is not supported")
        end
        return MAT_v5.matopen(rawfid, endian_indicator)
    end

    # Check for HDF5 file
    for offset = 512:512:fs-8
        seek(rawfid, offset)
        if read!(rawfid, Vector{UInt8}(undef, 8)) == HDF5_HEADER
            close(rawfid)
            return MAT_HDF5.matopen(filename, rd, wr, cr, tr, ff, compress)
        end
    end

    close(rawfid)
    error("\"$filename\" is not a MAT file")
end

function matopen(fname::AbstractString, mode::AbstractString; compress::Bool = false)
    mode == "r"  ? matopen(fname, true , false, false, false, false, false)    :
    mode == "r+" ? matopen(fname, true , true , false, false, false, compress) :
    mode == "w"  ? matopen(fname, false, true , true , true , false, compress) :
    # mode == "w+" ? matopen(fname, true , true , true , true , false, compress) :
    # mode == "a"  ? matopen(fname, false, true , true , false, true, compress)  :
    # mode == "a+" ? matopen(fname, true , true , true , false, true, compress)  :
    throw(ArgumentError("invalid open mode: $mode"))
end

matopen(fname::AbstractString; kwargs...) = matopen(fname, "r"; kwargs...)

function matopen(f::Function, args...; kwargs...)
    fid = matopen(args...; kwargs...)
    try
        f(fid)
    finally
        close(fid)
    end
end

"""
    matopen(filename [, mode]; compress = false) -> handle
    matopen(f::Function, filename [, mode]; compress = false) -> f(handle)

Mode defaults to `"r"` for read.
It can also be `"w"` for write,
or `"r+"` for read or write without creation or truncation.

Compression on reading is detected/handled automatically; the `compress`
keyword argument only affects write operations.

Use with `read`, `write`, `close`, `keys`, and `haskey`.

## Examples

### Usage with `read`

```julia-repl
julia> matwrite("test.mat", Dict("foo" => 1, "bar" => 2)) # create file
julia> file = matopen("test.mat")
julia> read(file) # get a dictionary of `key => value` pairs
julia> read(file, "foo") # read a specific variable
julia> close(file)
```

### Usage with `write`

```julia-repl
julia> matwrite("test.mat", Dict("foo" => 1, "bar" => 2)) # create file
julia> file = matopen("test.mat", "w+")
julia> write(file, "baz", 3) # test.mat now contains the variable "baz"
julia> close(file)
```

```julia-repl
julia> matwrite("test.mat", Dict("foo" => 1, "bar" => 2)) # create file
julia> file = matopen("test.mat", "w")
julia> write(file, "bop", 3) # test.mat now ONLY contains the variable "bop"
julia> close(file)
```

### Usage with `keys` and `haskey`

```julia-repl
julia> matwrite("test.mat", Dict("foo" => 1, "bar" => 2)) # create file
julia> file = matopen("test.mat")
julia> keys(file) # ["bar", "foo"]
julia> haskey(file, "bar") # true
julia> close(file)
```
"""
matopen

# Read all variables from a MATLAB file
"""
    matread(filename) -> Dict

Return a dictionary of all the variables and values in a Matlab file,
opening and closing it automatically.

### Example

```julia-repl
julia> matwrite("test.mat", Dict("foo" => 1, "bar" => 2)) # create file
julia> data = matread("test.mat")
julia> keys(data) # KeySet containing "foo" and "bar"
julia> data["foo"] # 1
```
"""
function matread(filename::AbstractString)
    file = matopen(filename)
    local vars
    try
        vars = read(file)
    finally
        close(file)
    end
    vars
end

# Write a dict to a MATLAB file
"""
    matwrite(filename, d::Dict; compress::Bool = false, version::String = "")

Write a dictionary containing variable names as keys and values as values
to a Matlab file, opening and closing it automatically.

### Arguments

- `filename`: The name of the output file. If `filename` already exists, \
its contents will be completely replaced.
- `d`: A dictionary of `key => value` pairs to write. 
  - `key` should be a `String` that starts with a letter and contain only alphanumeric characters and underscores.
  - `value` should be `<:Number`, `Array{<:Number}`, `String`, or `Bool`. `Tuple`, \
  `Nothing` and custom `struct`s will not work.

### Optional Arguments

- `compress`: If `true`, the output data are compressed. Default `false`.
- `version`: A `String` giving the target version. If equal to `"v4"`, \
write to Matlab v4, otherwise write to Matlab v7.3. Default `""`.

### Examples

Write a dictionary of various types:

```julia-repl
julia> d = Dict(
    "var1" => 1,
    "var2" => 4.0 + 5.0*im,
    "var3" => [1.0, 2.0, 3.0],
    "var4" => rand(UInt8, 4, 4),
    "var5" => true,
    "var6" => "foo")
julia> matwrite("test.mat", d)
```

Add a variable to an existing file:

```julia-repl
julia> d = Dict{String, Any}("var7" => "bar") # note use of `Any`
julia> merge!(d, matread("test.mat"))
julia> matwrite("test.mat", d) # "var7" is added to test.mat
```

Completely replace an existing file:

```julia-repl
julia> d = Dict{String, Any}("foo" => "bar")
julia> matwrite("test.mat", d) # "foo" is now the only variable in test.mat
```
"""
function matwrite(filename::AbstractString, dict::AbstractDict{S, T}; compress::Bool = false, version::String ="") where {S, T}

    if version == "v4"
        file = open(filename, "w")
        m = MAT_v4.Matlabv4File(file, false)
        try
            for (k, v) in dict
                local kstring
                try
                    kstring = ascii(convert(String, k))
                catch x
                    error("matwrite requires a Dict with ASCII keys")
                end
                write(m, kstring, v)
            end
        finally
            close(file)
        end

    else
        
        file = matopen(filename, "w"; compress = compress)
        try
            for (k, v) in dict
                local kstring
                try
                    kstring = ascii(convert(String, k))
                catch x
                    error("matwrite requires a Dict with ASCII keys")
                end
                write(file, kstring, v)
            end
        finally
            close(file)
        end

    end
end

###
### v0.10.0 deprecations
###

export exists
@noinline function exists(matfile::Union{MAT_v4.Matlabv4File,MAT_v5.Matlabv5File,MAT_HDF5.MatlabHDF5File}, varname::String)
    Base.depwarn("`exists(matfile, varname)` is deprecated, use `haskey(matfile, varname)` instead.", :exists)
    return haskey(matfile, varname)
end
@noinline function Base.names(matfile::Union{MAT_v4.Matlabv4File,MAT_v5.Matlabv5File,MAT_HDF5.MatlabHDF5File})
    Base.depwarn("`names(matfile)` is deprecated, use `keys(matfile)` instead.", :names)
    return keys(matfile)
end

end
