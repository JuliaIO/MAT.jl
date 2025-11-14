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

    export MatlabStructArray, StructArrayField, convert_struct_array
    export MatlabClassObject

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
    matwrite("matfile.mat", Dict("struct_array" => s_arr))
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
        function MatlabStructArray(names::Vector{String}, values::Vector{Array{Any,N}}, class::String="") where N
            # call MatlabStructArray{N}() to avoid the check
            check_struct_array(names, values)
            return new{N}(names, values, class)
        end
        function MatlabStructArray{N}(names::Vector{String}, values::Vector{Array{Any,N}}, class::String="") where N
            return new{N}(names, values, class)
        end
    end

    function check_struct_array(names::Vector{String}, values::Vector{Array{Any,N}}) where N
        if length(names) != length(values)
            error("MatlabStructArray requires equal number of names and values")
        end
        first_value, rest_values = Iterators.peel(values)
        first_len = length(first_value)
        if !all(x->length(x)==first_len, rest_values)
            error("MatlabStructArray requires all value columns to be of equal length")
        end
    end

    function MatlabStructArray(names::AbstractVector{<:AbstractString}, values::AbstractArray{<:AbstractArray{T,N}}, class="") where {T,N}
        MatlabStructArray{N}(string.(names), Vector{Array{Any,N}}(values), string(class))
    end

    # empty array
    function MatlabStructArray(names::AbstractVector{<:AbstractString}, dims::Tuple)
        N = length(dims)
        return MatlabStructArray{N}(names, [Array{Any, N}(undef, dims...) for n in names])
    end
    MatlabStructArray(names::AbstractVector{<:AbstractString}) = MatlabStructArray(names, (0,0))

    Base.eltype(::Type{MatlabStructArray{N}}) where N = Pair{String, Array{Any,N}}
    Base.length(arr::MatlabStructArray) = length(arr.names)

    function Base.iterate(arr::T, i=next_state(arr)) where T<:MatlabStructArray
        if i == 0 
            return nothing
        else
            return (eltype(T)(arr.names[i], arr.values[i]), next_state(arr,i))
        end
    end
    next_state(arr, i=0) = length(arr)==i ? 0 : i+1

    function Base.show(io::IO, ::MIME"text/plain", arr::MatlabStructArray)
        summary(io, arr)
        ncol = length(arr.values)
        print(io, " with $(ncol) ")
        col_word = ncol==1 ? "column" : "columns"
        print(io, col_word, ":")
        for (k,v) in arr
            print(io, "\n \"$k\": $v")
        end
    end

    function Base.:(==)(m1::MatlabStructArray{N},m2::MatlabStructArray{N}) where N
        return isequal(m1.names, m2.names) && isequal(m1.values, m2.values) && isequal(m1.class, m2.class)
    end

    function Base.isapprox(m1::MatlabStructArray,m2::MatlabStructArray; kwargs...)
        return isequal(m1.names, m2.names) && isapprox(m1.values, m2.values; kwargs...)
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

    # convert Dict array to MatlabStructArray
    function MatlabStructArray(arr::AbstractArray{<:AbstractDict, N}) where N
        first_keys = keys(first(arr))
        field_names = string.(first_keys)
        # Ensure same field set for all elements
        for d in arr
            if !issetequal(keys(d), first_keys)
                error("Cannot convert Dict array to MatlabStructArray. All elements must share identical field names")
            end
        end
        field_values = Vector{Array{Any,N}}(undef, length(field_names))
        for (idx,k) in enumerate(first_keys)
            this_field_values = Array{Any, N}(undef, size(arr))
            for (idx, d) in enumerate(arr)
                this_field_values[idx] = d[k]
            end
            field_values[idx] = this_field_values
        end
        return MatlabStructArray{N}(field_names, field_values)
    end

    function Base.Dict(arr::MatlabStructArray)
        return Base.Dict{String, Any}(arr)
    end
    function Base.Dict{String, Any}(arr::MatlabStructArray)
        Base.Dict{String, Any}(arr.names .=> arr.values)
    end

    Base.Array(arr::MatlabStructArray) = Array{Dict{String,Any}}(arr)
    function Base.Array{D}(arr::MatlabStructArray{N}) where {T,D<:AbstractDict{T},N}
        first_field = first(arr.values)
        sz = size(first_field)
        result = Array{D, N}(undef, sz)
        for idx in eachindex(first_field)
            element_values = (v[idx] for v in arr.values)
            result[idx] = D(T.(arr.names) .=> element_values)
        end
        return result
    end

    struct StructArrayField{N}
        values::Array{Any,N}
    end
    dimension(::StructArrayField{N}) where N = N

    """
        MatlabClassObject(
            d::Dict{String, Any},
            class::String,
        ) <: AbstractDict{String, Any}

    Type to store old class objects. Inside MATLAB a class named \"TestClassOld\" would be defined within `@TestClassOld` folders.

    If you want to write these objects you have to make sure the keys in the Dict match the class defined properties/fields.
    """
    struct MatlabClassObject <: AbstractDict{String, Any}
        d::Dict{String, Any}
        class::String
    end

    Base.eltype(::Type{MatlabClassObject}) = Pair{String, Any}
    Base.length(m::MatlabClassObject) = length(m.d)
    Base.keys(m::MatlabClassObject) = keys(m.d)
    Base.values(m::MatlabClassObject) = values(m.d)
    Base.getindex(m::MatlabClassObject, i) = getindex(m.d, i)
    Base.setindex!(m::MatlabClassObject, v, k) = setindex!(m.d, v, k)
    Base.iterate(m::MatlabClassObject, i) = iterate(m.d, i)
    Base.iterate(m::MatlabClassObject) = iterate(m.d)
    Base.haskey(m::MatlabClassObject, k) = haskey(m.d, k)
    Base.get(m::MatlabClassObject, k, default) = get(m.d, k, default)

    function convert_struct_array(d::Dict{String, Any}, class::String="")
        # there is no possibility of having cell arrays mixed with struct arrays (afaik)
        field_values = first(values(d))
        if field_values isa StructArrayField
            return MatlabStructArray{dimension(field_values)}(
                collect(keys(d)),
                [arr.values for arr in values(d)],
                class,
            )
        else
            if isempty(class)
                return d
            else
                return MatlabClassObject(d, class)
            end
        end
    end 
end