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

export MatlabStructArray, StructArrayField, convert_struct_array
export MatlabClassObject
export MatlabOpaque, convert_opaque
export MatlabTable
export ScalarOrArray

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
    return isequal(m1.class, m2.class) && isequal(m1.names, m2.names) && isapprox(m1.values, m2.values; kwargs...)
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
    keys(d1) == keys(d2) || return false
    for k in keys(d1)
        v1, v2 = d1[k], d2[k]
        value_isapprox(v1, v2) || return false
    end
    return true
end
dict_isapprox(d1::AbstractDict{T1}, d2::AbstractDict{T2}; kwargs...) where {T1,T2} = false

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
    return DateTime(1970, 1, 1) + Second(s) + Millisecond(ms_rem)
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

end