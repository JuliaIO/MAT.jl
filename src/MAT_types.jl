# MAT_types.jl
# Internal types used by MAT.jl
#
# Copyright (C) 2012   Matthijs Cox
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

###########################################
## Reading and writing MATLAB .mat files ##
###########################################

module MAT_types

using StringEncodings: StringEncodings
import StringEncodings: Encoding
import Dates
import Dates: DateTime, Second, Millisecond
import PooledArrays: PooledArray, RefArray
using Tables: Tables
import OrderedCollections: OrderedDict
using StructUtils

export MatlabStructArray, StructArrayField, convert_struct_array
export MatlabClassObject
export MatlabOpaque, convert_opaque
export MatlabTable
export ScalarOrArray
export FunctionHandle
export MATStyle, MATReadStyle, MATWriteStyle
export MATConstructionError

"""
    MATStyle <: StructUtils.StructStyle

The [StructUtils.jl](https://github.com/JuliaServices/StructUtils.jl) style MAT.jl
uses to convert between Julia values and MAT data. It has two concrete subtypes,
[`MATReadStyle`](@ref) and [`MATWriteStyle`](@ref).

You only need the style types to customize how your own types are converted, by
adding a method for a StructUtils function that dispatches on the style. A method
on `MATStyle` applies to both reading and writing; one on a subtype applies to
that direction only. Methods on a MAT style do not affect other packages built on
StructUtils, such as JSON.jl.

Field tags are namespaced under `mat`, so they can sit alongside tags for other
packages:

```julia
@tags struct Sample
    x::Float64 &(mat=(name="sample_x",), json=(name="x",))
end
```

# Example

A type MATLAB has no representation for is stored as a plain number by extending
`StructUtils.lower` for [`MATWriteStyle`](@ref), and read back by extending
`StructUtils.lift`. Declaring it not `structlike` tells StructUtils it is a single
value rather than a struct to be rebuilt field by field. This two-argument `lift`
takes no style, so unlike the other two methods it applies to every package built
on StructUtils:

```julia
using MAT

struct Celsius
    degrees::Float64
end

StructUtils.lower(::MATWriteStyle, x::Celsius) = x.degrees
StructUtils.structlike(::MATStyle, ::Type{Celsius}) = false
StructUtils.lift(::Type{Celsius}, x::Float64) = Celsius(x)
```
"""
abstract type MATStyle <: StructUtils.StructStyle end

"""
    MATReadStyle <: MATStyle

The StructUtils style used by the typed reads, `matread(filename, varname, T)` and
`read(matfile_handle, varname, T)`. Extend a StructUtils function on this style to
change how a type is built from MAT data without affecting how it is written, for
example `@choosetype MATReadStyle MyAbstractType x -> ...` to limit a type choice
to MAT reads. See [`MATStyle`](@ref).
"""
struct MATReadStyle <: MATStyle end

"""
    MATWriteStyle <: MATStyle

The StructUtils style used by `matwrite` and `write`. Extend `StructUtils.lower`
on this style to convert a value to something MATLAB can store before it is
written; the overload applies wherever the value appears, whether as a variable,
a struct field, a `Dict` value or a cell array element:

```julia
StructUtils.lower(::MATWriteStyle, x::Celsius) = x.degrees
```

See [`MATStyle`](@ref).
"""
struct MATWriteStyle <: MATStyle end

StructUtils.fieldtagkey(::MATStyle) = :mat

# MATLAB spells a missing value as the empty matrix `[]`, which is how `nothing`
# and `missing` are written. Note the asymmetry: these read back as an empty
# array unless the target field type is a `Union` with `Nothing`/`Missing`
# (see `StructUtils.nulllike` below).
const matlab_empty_matrix = zeros(Float64, 0, 0)

# MATLAB stores complex numbers natively, so a `Complex` is a leaf value here
# rather than a two-field struct to be rebuilt from its `re`/`im` components.
StructUtils.structlike(::MATStyle, ::Type{<:Complex}) = false

# MATLAB has no scalar character type distinct from a 1x1 char array, so MAT.jl
# reads a one-character string back as a `Char`. Lift it for string fields.
StructUtils.lift(st::MATStyle, ::Type{T}, x::AbstractChar) where {T<:AbstractString} =
    (convert(T, string(x)), StructUtils.defaultstate(st))

# The other end of the same mismatch: an empty MATLAB char array reads back as an
# empty array rather than as `""`, an empty `Vector{String}` in v5 and v7.3. An
# empty array whose file records no element type (every empty v4 array, which is
# how v4 stores both `''` and `[]`, and empty cells in any version) has element
# type `Union{}`. The signature admits both, since `Union{}` is a subtype of every
# type, and rejects empty numeric data, which is not a string.
function StructUtils.lift(st::MATStyle, ::Type{T},
                          x::AbstractArray{<:Union{AbstractString,AbstractChar}}) where {T<:AbstractString}
    isempty(x) || throw(ArgumentError(
        "cannot read a $(join(size(x), "x")) MATLAB array as `$T`"))
    return (convert(T, ""), StructUtils.defaultstate(st))
end

# MATLAB has no scalars, and v4 files read a scalar back as a 1x1 matrix (v5 and
# v7.3 collapse it on read). A one-element array therefore fills a numeric
# target, the counterpart of a scalar filling a one-element array target.
function StructUtils.lift(st::MATStyle, ::Type{T}, x::AbstractArray{<:Number}) where {T<:Number}
    length(x) == 1 || throw(ArgumentError(
        "cannot read a $(join(size(x), "x")) MATLAB array as `$T`"))
    return (convert(T, only(x)), StructUtils.defaultstate(st))
end

# MATLAB spells "no value" as the empty matrix `[]`, so an empty array narrows a
# `Union{Nothing,T}` or `Union{Missing,T}` field to `nothing`/`missing`. The
# ambiguity is inherent to the format: a genuinely empty numeric array stored in
# such a field is indistinguishable from a missing value.
#
# Empty *string* arrays are excluded: they are how `""` comes back, and reading
# `""` as `nothing` would silently discard a value the file does hold. This
# distinction survives only where the file records the element type; a v4 file
# stores `""` and `[]` identically, as an array of element type `Union{}`, so
# there `""` still narrows to `nothing`. `Union{}` is tested for explicitly
# because, being a subtype of every type, it would otherwise pass as string data.
StructUtils.nulllike(::MATStyle, x::AbstractArray) =
    isempty(x) && (eltype(x) === Union{} || !(eltype(x) <: Union{AbstractString,AbstractChar}))

const ScalarOrArray{T} = Union{T, AbstractArray{T}}

# struct arrays are stored as columns per field name
"""
    MatlabStructArray{N}(
        names::Vector{String},
        values::Vector{Array{Any,N}},
        class::String = "",
    )

Data structure to store matlab struct arrays, which stores the field names separate from the field values.
The field values are stored as columns of `Array{Any,N}` per Matlab field, which is how MAT files store these structures.

These are distinct from cell arrays of structs,
which are handled as in MAT.jl as `Array{Any,N}` with `Dict{String,Any}` inside,
for example `Any[Dict("x"=>1), Dict("x"=>2)]`.

Old class object arrays can be handled by providing a non-empty class name.

# Example

```julia
using MAT

s_arr = MatlabStructArray(["a", "b"], [[1, 2],["foo", 5]])

# write-read
matwrite("matfile.mat", Dict{String,Any}("struct_array" => s_arr))
read_s_arr = matread("matfile.mat")["struct_array"]

# convert to Dict Array
dict_array = Array{Dict{String,Any}}(s_arr)

# convert to Dict (with arrays as fields)
dict = Dict{String,Any}(s_arr)
```
"""
struct MatlabStructArray{N}
    names::Vector{String}
    values::Vector{Array{Any,N}}
    class::String
    function MatlabStructArray(
        names::Vector{String},
        values::Vector{Array{Any,N}},
        class::String="";
        check::Bool=true,
    ) where {N}
        check && check_struct_array(names, values)
        return new{N}(names, values, class)
    end
    function MatlabStructArray{N}(
        names::Vector{String}, values::Vector{Array{Any,N}}, class::String=""
    ) where {N}
        return new{N}(names, values, class)
    end
end

function check_struct_array(names::Vector{String}, values::Vector{Array{Any,N}}) where {N}
    if length(names) != length(values)
        error("MatlabStructArray requires equal number of names and values")
    end
    first_value, rest_values = Iterators.peel(values)
    first_len = length(first_value)
    if !all(x -> length(x) == first_len, rest_values)
        error("MatlabStructArray requires all value columns to be of equal length")
    end
end

function MatlabStructArray(
    names::AbstractVector{<:AbstractString},
    values::AbstractArray{A},
    class="";
    check::Bool=true,
) where {N,A<:AbstractArray{T,N} where {T}}
    return MatlabStructArray(
        string.(names), Vector{Array{Any,N}}(values), string(class); check=check
    )
end
function MatlabStructArray(
    names::Vector{String}, values::AbstractArray{A}, class=""; check::Bool=true
) where {N,A<:AbstractArray{T,N} where {T}}
    return MatlabStructArray(
        names, Vector{Array{Any,N}}(values), string(class); check=check
    )
end

# empty array
function MatlabStructArray(names::AbstractVector{<:AbstractString}, dims::Tuple)
    N = length(dims)
    return MatlabStructArray{N}(names, [Array{Any,N}(undef, dims...) for n in names])
end
function MatlabStructArray(names::AbstractVector{<:AbstractString})
    return MatlabStructArray(names, (0, 0))
end

Base.eltype(::Type{MatlabStructArray{N}}) where {N} = Pair{String,Array{Any,N}}
Base.length(arr::MatlabStructArray) = length(arr.names)
Base.keys(arr::MatlabStructArray) = arr.names
Base.values(arr::MatlabStructArray) = arr.values
Base.haskey(arr::MatlabStructArray, k::AbstractString) = k in keys(arr)
function Base.copy(arr::MatlabStructArray{N}) where {N}
    return MatlabStructArray{N}(copy(arr.names), copy(arr.values))
end
class(arr::MatlabStructArray) = arr.class
class(d::AbstractDict) = ""

function Base.iterate(arr::T, i=next_state(arr)) where {T<:MatlabStructArray}
    if i == 0
        return nothing
    else
        return (eltype(T)(arr.names[i], arr.values[i]), next_state(arr, i))
    end
end
next_state(arr, i=0) = length(arr) == i ? 0 : i + 1

function Base.show(io::IO, ::MIME"text/plain", arr::MatlabStructArray)
    summary(io, arr)
    ncol = length(arr.values)
    print(io, " with $(ncol) ")
    col_word = ncol == 1 ? "column" : "columns"
    print(io, col_word, ":")
    for (k, v) in arr
        print(io, "\n \"$k\": ")
        summary(io, v)
    end
end

function Base.:(==)(m1::MatlabStructArray{N}, m2::MatlabStructArray{N}) where {N}
    return isequal(m1.names, m2.names) &&
           isequal(m1.values, m2.values) &&
           isequal(m1.class, m2.class)
end

function Base.isapprox(m1::MatlabStructArray, m2::MatlabStructArray; kwargs...)
    return isequal(m1.class, m2.class) && issetequal(m1.names, m2.names) &&
    key_based_isapprox(m1.names, m1, m2; kwargs...)
end

function find_index(m::MatlabStructArray, s::AbstractString)
    idx = findfirst(isequal(s), m.names)
    if isnothing(idx)
        error("field \"$s\" not found in MatlabStructArray")
    end
    return idx
end

function Base.getindex(m::MatlabStructArray, s::AbstractString)
    idx = find_index(m, s)
    return getindex(m.values, idx)
end

function Base.setindex!(m::MatlabStructArray{N1}, input::AbstractArray{T,N2}, s::AbstractString) where {N1,N2,T}
    error("size of input is not equal to MatlabStructArray column size")
end

function Base.setindex!(m::MatlabStructArray{N}, input::AbstractArray{T,N}, s::AbstractString) where {N,T}
    idx = find_index(m, s)
    v = Array{Any,N}(input)
    if size(v) != size(m.values[idx])
        error("size of input is not equal to MatlabStructArray column size")
    end
    return setindex!(m.values, v, idx)
end

function Base.get(m::MatlabStructArray, s::AbstractString, default)
    idx = findfirst(isequal(s), m.names)
    if isnothing(idx)
        return default
    else
        return getindex(m.values, idx)
    end
end

# convert Dict array to MatlabStructArray
function MatlabStructArray(arr::AbstractArray{<:AbstractDict,N}, class::String="") where {N}
    first_dict, remaining_dicts = Iterators.peel(arr)
    first_keys = keys(first_dict)
    field_names = string.(first_keys)
    # Ensure same field set for all elements
    for d in remaining_dicts
        if !issetequal(keys(d), first_keys)
            error(
                "Cannot convert Dict array to MatlabStructArray. All elements must share identical field names",
            )
        end
    end
    field_values = Vector{Array{Any,N}}(undef, length(field_names))
    for (idx, k) in enumerate(first_keys)
        this_field_values = Array{Any,N}(undef, size(arr))
        for (idx, d) in enumerate(arr)
            this_field_values[idx] = d[k]
        end
        field_values[idx] = this_field_values
    end
    return MatlabStructArray{N}(field_names, field_values, class)
end

function Base.Dict(arr::MatlabStructArray)
    return Base.Dict{String,Any}(arr)
end
function Base.Dict{String,Any}(arr::MatlabStructArray)
    return Base.Dict{String,Any}(arr.names .=> arr.values)
end

Base.Array{D}(arr::MatlabStructArray{N}) where {D<:AbstractDict,N} = Array{D,N}(arr)

function Base.Array{D,N}(arr::MatlabStructArray{N}) where {D<:AbstractDict,N}
    first_field = first(arr.values)
    sz = size(first_field)
    result = Array{D,N}(undef, sz)
    for idx in eachindex(first_field)
        element_values = (v[idx] for v in arr.values)
        result[idx] = create_struct(D, arr.names, element_values, arr.class)
    end
    return result
end

function Base.getindex(arr::MatlabStructArray, s::Integer)
    row_values = [getindex(v, s) for v in arr.values]
    if isempty(arr.class)
        D = Dict{String,Any}
    else
        D = OrderedDict{String,Any}
    end
    return create_struct(D, arr.names, row_values, arr.class)
end

function create_struct(::Type{D}, keys, values, class::String) where {T,D<:AbstractDict{T}}
    return D(T.(keys) .=> values)
end

# 1D MatlabStructArray also counts as table (mostly for testing purposes)
Tables.istable(::Type{MatlabStructArray{1}}) = true
Tables.columns(t::MatlabStructArray{1}) = Symbol.(t.values)
Tables.columnnames(t::MatlabStructArray{1}) = t.names
Tables.getcolumn(t::MatlabStructArray{1}, nm::String) = t[nm]
Tables.getcolumn(t::MatlabStructArray{1}, nm::Symbol) = Tables.getcolumn(t, string(nm))
function MatlabStructArray{1}(t::Tables.CopiedColumns)
    col_names = Tables.columnnames(t)
    return MatlabStructArray{1}(
        string.(col_names), [Vector{Any}(Tables.getcolumn(t, nm)) for nm in col_names]
    )
end
MatlabStructArray(t::Tables.CopiedColumns) = MatlabStructArray{1}(t)

struct StructArrayField{N}
    values::Array{Any,N}
end
dimension(::StructArrayField{N}) where {N} = N

"""
Internal Marker for Empty Structs with dimensions like 1x0 or 0x0
"""
struct EmptyStruct
    dims::Vector{UInt64}
end
class(m::EmptyStruct) = ""

"""
    FunctionHandle(d::Dict{String, Any}) <: AbstractDict{String, Any}

Type to store function handles which are stored as structs in MATLAB.
"""
struct FunctionHandle{D<:AbstractDict{String,Any}} <: AbstractDict{String,Any}
    d::D
end

Base.eltype(::Type{FunctionHandle}) = Pair{String,Any}
Base.length(m::FunctionHandle) = length(m.d)
Base.keys(m::FunctionHandle) = keys(m.d)
Base.values(m::FunctionHandle) = values(m.d)
Base.getindex(m::FunctionHandle, i) = getindex(m.d, i)
Base.setindex!(m::FunctionHandle, v, k) = setindex!(m.d, v, k)
Base.iterate(m::FunctionHandle, i) = iterate(m.d, i)
Base.iterate(m::FunctionHandle) = iterate(m.d)
Base.haskey(m::FunctionHandle, k) = haskey(m.d, k)
Base.get(m::FunctionHandle, k, default) = get(m.d, k, default)

function Base.:(==)(m1::FunctionHandle, m2::FunctionHandle)
    return m1.d == m2.d
end

function Base.isapprox(m1::FunctionHandle, m2::FunctionHandle; kwargs...)
    return dict_isapprox(m1.d, m2.d; kwargs...)
end

"""
    MatlabClassObject(
        d::Dict{String, Any},
        class::String,
    ) <: AbstractDict{String, Any}

Type to store old class objects. Inside MATLAB a class named \"TestClassOld\" would be defined within `@TestClassOld` folders.

If you want to write these objects you have to make sure the keys in the Dict match the class defined properties/fields.
"""
struct MatlabClassObject{D<:AbstractDict{String,Any}} <: AbstractDict{String,Any}
    d::D
    class::String
end

class(m::MatlabClassObject) = m.class

Base.eltype(::Type{MatlabClassObject}) = Pair{String,Any}
Base.length(m::MatlabClassObject) = length(m.d)
Base.keys(m::MatlabClassObject) = keys(m.d)
Base.values(m::MatlabClassObject) = values(m.d)
Base.getindex(m::MatlabClassObject, i) = getindex(m.d, i)
Base.setindex!(m::MatlabClassObject, v, k) = setindex!(m.d, v, k)
Base.iterate(m::MatlabClassObject, i) = iterate(m.d, i)
Base.iterate(m::MatlabClassObject) = iterate(m.d)
Base.haskey(m::MatlabClassObject, k) = haskey(m.d, k)
Base.get(m::MatlabClassObject, k, default) = get(m.d, k, default)

function Base.:(==)(m1::MatlabClassObject, m2::MatlabClassObject)
    return m1.class == m2.class && m1.d == m2.d
end

function Base.isapprox(m1::MatlabClassObject, m2::MatlabClassObject; kwargs...)
    return m1.class == m2.class && dict_isapprox(m1.d, m2.d; kwargs...)
end

function MatlabStructArray(arr::AbstractArray{<:MatlabClassObject})
    first_obj, remaining_obj = Iterators.peel(arr)
    class = first_obj.class
    if !all(x -> isequal(class, x.class), remaining_obj)
        error(
            "to write a MatlabClassObject array all classes must be equal. Use `Array{Any}` to write a cell array",
        )
    end
    return MatlabStructArray(arr, class)
end

function convert_struct_array(d::AbstractDict{String,Any}, class::String="")
    # there is no possibility of having cell arrays mixed with struct arrays (afaik)
    field_values = first(values(d))
    if field_values isa StructArrayField
        return MatlabStructArray{dimension(field_values)}(
            collect(keys(d)), [arr.values for arr in values(d)], class
        )
    else
        if isempty(class)
            return d
        elseif class == "function_handle"
            return FunctionHandle(d)
        else
            return MatlabClassObject(d, class)
        end
    end
end

function Base.Array(arr::MatlabStructArray{N}) where {N}
    if isempty(arr.class)
        return Array{Dict{String,Any},N}(arr)
    else
        # ordered dict by default, to preserve field order
        D = OrderedDict{String,Any}
        return Array{MatlabClassObject{D},N}(arr)
    end
end

function create_struct(::Type{MatlabClassObject{D}}, keys, values, class::String) where {D<:AbstractDict}
    d = D(string.(keys) .=> values)
    return MatlabClassObject{D}(d, class)
end

"""
    MatlabOpaque(
        d::Dict{String, Any},
        class::String,
    ) <: AbstractDict{String, Any}

Type to store opaque class objects.
These are the 'modern' Matlab classes, different from the old `MatlabClassObject` types.

"""
struct MatlabOpaque{D<:AbstractDict{String,Any}} <: AbstractDict{String,Any}
    d::D
    class::String
end
class(m::MatlabOpaque) = m.class

Base.eltype(::Type{MatlabOpaque}) = Pair{String,Any}
Base.length(m::MatlabOpaque) = length(m.d)
Base.keys(m::MatlabOpaque) = keys(m.d)
Base.values(m::MatlabOpaque) = values(m.d)
Base.getindex(m::MatlabOpaque, i) = getindex(m.d, i)
Base.setindex!(m::MatlabOpaque, v, k) = setindex!(m.d, v, k)
Base.iterate(m::MatlabOpaque, i) = iterate(m.d, i)
Base.iterate(m::MatlabOpaque) = iterate(m.d)
Base.haskey(m::MatlabOpaque, k) = haskey(m.d, k)
Base.get(m::MatlabOpaque, k, default) = get(m.d, k, default)

function Base.:(==)(m1::MatlabOpaque, m2::MatlabOpaque)
    return m1.class == m2.class && m1.d == m2.d
end

function Base.isapprox(m1::MatlabOpaque, m2::MatlabOpaque; kwargs...)
    return m1.class == m2.class && dict_isapprox(m1.d, m2.d; kwargs...)
end

function dict_isapprox(d1::AbstractDict{T}, d2::AbstractDict{T}; kwargs...) where T
    issetequal(keys(d1), keys(d2)) || return false
    return key_based_isapprox(keys(d1), d1, d2; kwargs...)
end
dict_isapprox(d1::AbstractDict{T1}, d2::AbstractDict{T2}; kwargs...) where {T1,T2} = false

function key_based_isapprox(keys, collection1, collection2; kwargs...)
    for k in keys
        v1, v2 = collection1[k], collection2[k]
        value_isapprox(v1, v2; kwargs...) || return false
    end
    return true
end

value_isapprox(x1::AbstractString, x2::AbstractString; kwargs...) = isequal(x1, x2)
value_isapprox(x1, x2; kwargs...) = isapprox(x1, x2; kwargs...)
value_isapprox(d1::AbstractDict, d2::AbstractDict; kwargs...) = dict_isapprox(x1, x2; kwargs...)

function convert_opaque(obj::MatlabOpaque; table::Type=Nothing)
    if obj.class == "string"
        return from_string(obj)
    elseif obj.class == "datetime"
        return from_datetime(obj)
    elseif obj.class == "duration"
        return from_duration(obj)
    elseif obj.class == "categorical"
        return from_categorical(obj)
    elseif obj.class == "table"
        return from_table(obj, table)
    else
        return obj
    end
end

# for reference: https://github.com/foreverallama/matio/blob/main/matio/utils/converters/matstring.py
function from_string(obj::MatlabOpaque, encoding::Encoding=Encoding(Symbol("UTF-16LE")))
    data = obj["any"]
    if isnothing(data) || isempty(data)
        return String[]
    end
    if data[1, 1] != 1
        @warn "String saved from a different MAT-file version. Returning empty string"
        return ""
    end
    ndims = data[1, 2]
    shape = Int.(data[1, 3:(2 + ndims)])
    num_strings = prod(shape)
    char_counts = data[1, (3 + ndims):(2 + ndims + num_strings)]
    byte_data = data[1, (3 + ndims + num_strings):end]
    bytes = reinterpret(UInt8, byte_data)

    strings = String[]
    pos = 1

    for char_count in char_counts
        byte_length = char_count * 2  # UTF-16 encoding
        extracted_bytes = bytes[pos:(pos + byte_length - 1)]
        str = StringEncodings.decode(extracted_bytes, encoding)
        push!(strings, str)
        pos += byte_length
    end

    if num_strings == 1
        return first(strings)
    else
        return reshape(strings, shape...)
    end
end

function from_datetime(obj::MatlabOpaque)
    dat = obj["data"]
    if isnothing(dat) || isempty(dat)
        return DateTime[]
    end
    if haskey(obj, "tz") && !isempty(obj["tz"])
        tz = obj["tz"]
        @warn "no timezone conversion yet for datetime objects. timezone of \"$tz\" ignored"
    end
    #isdate = obj["isDateOnly"] # optional: convert to Date instead of DateTime?
    return map_or_not(ms_to_datetime, dat)
end

# is the complex part the submilliseconds?
ms_to_datetime(ms::Complex) = ms_to_datetime(real(ms))
function ms_to_datetime(ms::Real)
    s, ms_rem = fldmod(ms, 1_000)  # whole seconds and remainder milliseconds
    return DateTime(1970, 1, 1) + Second(s) + Millisecond(round(Int, ms_rem))
end

function to_matlab_data(d::DateTime)
    ms = Dates.value(d)

    # matlab counts w.r.t. 1970
    matlab_offset = Dates.value(DateTime(1970, 1, 1))

    # note: can be negative, that's fine
    matlab_ms = ms - matlab_offset
    return Complex(Float64(matlab_ms))
end

function MatlabOpaque(d::ScalarOrArray{DateTime})
    return MatlabOpaque(
        Dict{String,Any}(
            "tz"   => "",
            "data" => map_or_not(to_matlab_data, d),
            "fmt"  => "",
        ),
        "datetime"
    )
end

function from_duration(obj::MatlabOpaque)
    dat = obj["millis"]
    #fmt = obj["fmt"] # TODO: format, e.g. 'd' to Day
    if isnothing(dat) || isempty(dat)
        return Millisecond[]
    end
    return map_or_not(Millisecond, dat)
end

to_matlab_millis(d::Millisecond) = Float64(Dates.value(d))

function MatlabOpaque(d::ScalarOrArray{Millisecond})
    return MatlabOpaque(
        Dict{String,Any}(
            "millis" => map_or_not(to_matlab_millis, d),
        ),
        "duration"
    )
end

function from_categorical(obj::MatlabOpaque)
    category_names = obj["categoryNames"]
    codes = obj["codes"]
    pool = vec(Array{promoted_eltype(category_names)}(category_names))
    code_type = eltype(codes)
    pool_type = eltype(pool)
    invpool = Dict{pool_type,code_type}(pool .=> code_type.(1:length(pool)))
    RA = typeof(codes)
    N = ndims(codes)
    return PooledArray{pool_type,code_type,N,RA}(RefArray(codes), invpool, pool)
end

function promoted_eltype(v::AbstractArray{Any})
    isempty(v) && return T
    first_el, remaining = Iterators.peel(v)
    T_out = typeof(first_el)
    for el in remaining
        T_out = promote_type(T_out, typeof(el))
    end
    return T_out
end
promoted_eltype(::AbstractArray{T}) where {T} = T

map_or_not(f, dat::AbstractArray) = map(f, dat)
map_or_not(f, dat) = f(dat)

struct MatlabTable
    names::Vector{Symbol}
    columns::Vector
end
Tables.istable(::Type{MatlabTable}) = true
Tables.columns(t::MatlabTable) = t.columns
Tables.columnnames(t::MatlabTable) = t.names
Tables.getcolumn(t::MatlabTable, nm::Symbol) = t[nm]
function find_index(m::MatlabTable, s::Symbol)
    idx = findfirst(isequal(s), m.names)
    if isnothing(idx)
        error("column :$s not found in MatlabTable")
    end
    return idx
end
function Base.getindex(m::MatlabTable, s::Symbol)
    idx = find_index(m, s)
    return getindex(m.columns, idx)
end
Base.getindex(m::MatlabTable, s::AbstractString) = getindex(m, Symbol(s))
MatlabTable(t::Tables.CopiedColumns{MatlabTable}) = Tables.source(t)

function from_table(obj::MatlabOpaque, ::Type{T}=MatlabTable) where {T}
    names = vec(Symbol.(obj["varnames"]))
    cols = vec([try_vec(c) for c in obj["data"]])
    t = MatlabTable(names, cols)
    return T(Tables.CopiedColumns(t))
end
# option to not convert and get the MatlabOpaque as table
from_table(obj::MatlabOpaque, ::Type{Nothing}) = obj

try_vec(c::Vector) = c
try_vec(c) = [c]
function try_vec(c::AbstractArray)
    return (size(c, 2) == 1) ? vec(c) : c
end

# MATLAB char array storage types:
# - miUINT8: ASCII or UTF-8. Can be directly converted to String.
# - miUINT16: UTF-8 encoding packed in 16-bit values. Reinterpret as uint8 and remove 0x00 MSB bytes
# - miUTF8/16/32: UTF-8/16/32 encoding. Transcode to UTF-8 if needed, then convert to String
_decode_row(row::AbstractVector{UInt8},  codec)  = String(row)
_decode_row(row::AbstractVector{UInt32}, codec) = String(transcode(UInt8, row))
function _decode_row(row::AbstractVector{UInt16}, codec)
    if codec == :utf16
        return String(transcode(UInt8, row))
    end

    n = length(row)
    out = Vector{UInt8}(undef, 2n)
    j = 1

    @inbounds for x in row
        msb = UInt8(x >> 8)
        lsb = UInt8(x & 0xff)

        if msb == 0x00
            out[j] = lsb
            j += 1
        else
            out[j] = msb
            out[j+1] = lsb
            j += 2
        end
    end

    resize!(out, j - 1)
    return String(out)
end

function decode_char_array(arr::AbstractArray{T}, codec::Symbol) where T <: Union{UInt8, UInt16, UInt32}
    char_len = size(arr, 2)
    other_dims = (size(arr, 1), size(arr)[3:end]...)

    if char_len == 0
        return fill("", other_dims)
    end

    if ndims(arr) == 2 && size(arr, 1) == 1
        # Return as single string instead of 1-element array
        row = view(arr, 1, :)
        return _decode_row(row, codec)
    end

    out = Array{String}(undef, other_dims...)
    for I in CartesianIndices(out)
        idx = Tuple(I)
        row = view(arr, idx[1], :, idx[2:end]...)
        out[I] = _decode_row(row, codec)
    end

    return out
end

##############################################
## Constructing typed values from MAT data  ##
##############################################

# A value read from a MAT file is already a fully-typed Julia object (an `Array`,
# `Dict`, `String`, `MatlabStructArray`, ...) rather than a generic token stream,
# so the generic StructUtils source-iteration machinery does not fit it:
#
#   * MATLAB has no scalars. MAT.jl collapses a 1x1 array to a scalar on read, so
#     an array-typed field very often sees a scalar source.
#   * MAT arrays are flat and carry their own shape, while StructUtils infers
#     dimensions by descending into nested sources, which would flatten an NxM
#     matrix source into an N*M element vector.
#   * The element type comes from the file (`Float32`, `Int32`, or `Any` for cell
#     arrays) and rarely matches the target element type exactly.
#
# Array targets are therefore constructed here directly from the source array,
# converting element type and shape explicitly.

"""
    MATConstructionError(T, cause)

Raised by the typed read methods when a value read from a MAT file cannot be
constructed as the requested type `T`. The underlying `cause` is reported as
part of the error message.
"""
struct MATConstructionError <: Exception
    T::Type
    cause::Any
end

function Base.showerror(io::IO, e::MATConstructionError)
    print(io, "MATConstructionError: could not construct `", e.T, "` from MAT data.")
    if e.cause isa TypeError
        # how a field with neither a stored value nor a default surfaces
        print(io, "\nA required field may be missing from the file, or hold an incompatible type.")
    end
    print(io, "\nCaused by: ")
    if e.cause isa Exception
        showerror(io, e.cause)
    else
        print(io, e.cause)
    end
end

"""
    construct_from_raw(raw, ::Type{T}) -> T

Convert a raw MAT value (`Dict`, `MatlabStructArray`, `Array`, ...) into type `T`
via StructUtils. Shared by the `MAT_HDF5`, `MAT_v5` and `MAT_v4` typed-read
methods. Failures are wrapped in a [`MATConstructionError`](@ref).

This is internal to MAT.jl; use `read(matfile_handle, varname, T)` or
`matread(filename, varname, T)` instead.
"""
function construct_from_raw(raw, ::Type{T}) where {T}
    try
        return StructUtils.make(T, prepare_source(raw, T), MATReadStyle())::T
    catch e
        if e isa InterruptException || e isa StackOverflowError || e isa OutOfMemoryError
            rethrow()
        end
        throw(MATConstructionError(T, e))
    end
end

"""
    prepare_source(raw, ::Type{T})

Unwrap a source before it reaches `StructUtils.make`: a whole variable, or an
element of an array being constructed.

A one-element struct array is a valid source for a single `T`, and the `make`
method below implements that. It cannot be the only implementation, though:
being general in the target type, it is ambiguous with the per-type `make`
methods `@choosetype` defines, so a `@choosetype` target built from a struct
array would fail on dispatch. Unwrapping here happens before dispatch, which
keeps that case working.
"""
prepare_source(raw, ::Type{T}) where {T} = raw

prepare_source(raw::MatlabStructArray, ::Type{T}) where {T} =
    single_struct_target(MATReadStyle(), T) ? only_struct(raw, T) : raw

function only_struct(x::MatlabStructArray, ::Type{T}) where {T}
    n = struct_array_length(x)
    n == 1 || error(
        "cannot construct a single `$T` from a $n-element MATLAB struct array; " *
        "read it as an array type such as `Vector{$T}`",
    )
    return x[1]
end

# Fast path: the source is already exactly the requested array type.
StructUtils.make(st::MATReadStyle, ::Type{T}, x::T) where {T<:AbstractArray} =
    (x, StructUtils.defaultstate(st))

StructUtils.make(st::MATReadStyle, ::Type{A}, x) where {A<:AbstractArray} =
    make_mat_array(st, A, x)

# A `MatlabStructArray` is the columnar encoding MAT files use for struct arrays;
# expand it into the equivalent array of `Dict`s before constructing the target.
StructUtils.make(st::MATReadStyle, ::Type{A}, x::MatlabStructArray) where {A<:AbstractArray} =
    make_mat_array(st, A, struct_array_dicts(x))

# A struct array is also a valid source for a single struct, if it holds exactly
# one element. Any other target is a collection of its elements: expand the
# struct array and let the generic machinery build it, which also lets `Union`
# targets be narrowed on the expanded source.
StructUtils.make(st::MATReadStyle, ::Type{T}, x::MatlabStructArray) where {T} =
    StructUtils.make(st, T, single_struct_target(st, T) ? only_struct(x, T) : struct_array_dicts(x))

# Whether `make` would build `T` as a single struct rather than as a collection.
# The order mirrors the dispatch in `StructUtils.make` itself, where a type that
# is both (an `AbstractSet` is `arraylike` and a struct type) counts as a collection.
#
# An abstract target counts as a single struct: it is not `structlike`, but the
# concrete type it resolves to (via `@choosetype`, say) will be. `Any` is the
# exception, and keeps the struct array itself, matching the untyped read.
single_struct_target(st::MATReadStyle, ::Type{T}) where {T} =
    T !== Any && !(T <: Tuple) && !StructUtils.arraylike(st, T) &&
    (StructUtils.dictlike(st, T) || StructUtils.noarg(st, T) ||
     StructUtils.structlike(st, T) || isabstracttype(T))

# An `Any` target keeps the raw representation, matching the untyped read.
StructUtils.make(st::MATReadStyle, ::Type{Any}, x::MatlabStructArray) =
    (x, StructUtils.defaultstate(st))

struct_array_length(x::MatlabStructArray) =
    isempty(x.values) ? 0 : length(first(x.values))

function struct_array_dicts(x::MatlabStructArray{N}) where {N}
    isempty(x.values) && return Array{Dict{String,Any},N}(undef, ntuple(_ -> 0, N))
    return Array{Dict{String,Any}}(x)
end

function make_mat_array(st::MATReadStyle, ::Type{A}, x) where {A<:AbstractArray}
    if A isa Union
        error("cannot read a MAT variable as the union of array types `$A`; request a single array type")
    end
    x isa AbstractArray && return build_mat_array(st, A, x)
    # MATLAB has no scalars: a scalar source is a 1x1x...x1 array that MAT.jl
    # collapsed on read, so wrap it back up before converting.
    return build_mat_array(st, A, fill(x, ntuple(_ -> 1, target_ndims(A, 1))))
end

function build_mat_array(st::MATReadStyle, ::Type{A}, src::AbstractArray) where {A<:AbstractArray}
    ET = eltype(A)
    N = target_ndims(A, ndims(src))
    if N == 0 && length(src) != 1
        # a 0-dimensional target holds exactly one element
        throw(DimensionMismatch(shape_message(A, size(src), ndims(src))))
    end
    out = Array{ET,N}(undef, target_dims(A, N, src))
    for (i, j) in zip(eachindex(out), eachindex(src))
        val, _ = StructUtils.make(st, ET, prepare_source(src[j], ET))
        out[i] = val
    end
    return (out isa A ? out : convert(A, out)), StructUtils.defaultstate(st)
end

# `Array{Float64}` and `AbstractArray` leave the dimensionality open, in which
# case the source decides it. `Matrix` and friends are `UnionAll`s too, but their
# dimensionality is known, so ask whether `ndims` applies rather than testing for
# a `UnionAll`.
target_ndims(::Type{A}, unsized::Int) where {A} =
    applicable(ndims, A) ? ndims(A) : unsized

function target_dims(::Type{A}, N::Int, src::AbstractArray) where {A}
    sz = size(src)
    length(sz) == N && return sz
    if N == 1
        # MATLAB stores vectors as 1xN or Nx1, so a source with a single
        # non-singleton dimension flattens unambiguously. Anything else is a
        # shape mismatch that should not be silently flattened.
        n = length(src)
        if n > 1 && count(!isone, sz) > 1
            throw(DimensionMismatch(shape_message(A, sz, ndims(src))))
        end
        return (n,)
    elseif length(sz) < N
        return (sz..., ntuple(_ -> 1, N - length(sz))...)  # pad trailing singleton dimensions
    elseif all(isone, sz[(N + 1):end])
        return sz[1:N]  # drop trailing singleton dimensions
    end
    throw(DimensionMismatch(shape_message(A, sz, ndims(src))))
end

function shape_message(::Type{A}, sz::Tuple, M::Int) where {A}
    return "cannot read a $(join(sz, "x")) MATLAB array as `$A`; " *
           "use a matching array type such as `Array{$(eltype(A)),$M}`"
end

end