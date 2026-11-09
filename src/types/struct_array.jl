# struct_array.jl
# MatlabStructArray: the columnar form MAT files use for struct arrays

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
        throw(KeyError(s))
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
Tables.columnaccess(::Type{MatlabStructArray{1}}) = true
Tables.columns(t::MatlabStructArray{1}) = t
Tables.columnnames(t::MatlabStructArray{1}) = Symbol.(t.names)
Tables.getcolumn(t::MatlabStructArray{1}, i::Int) = t.values[i]
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

# A struct array with no fields, which MATLAB stores as its dimensions alone. A
# `MatlabStructArray` keeps its size in its field columns, so it cannot hold one;
# read it as an array of empty Dicts, or a single empty Dict when it is 1x1.
function fieldless_struct_array(dims::Tuple)
    dims == (1, 1) && return Dict{String,Any}()
    return [Dict{String,Any}() for _ in CartesianIndices(dims)]
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
