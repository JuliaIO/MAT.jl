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

    # struct arrays are stored as columns per field name
    struct MatlabStructArray{N}
        names::Vector{String}
        values::Vector{Array{Any,N}}
    end

    function Base.isequal(m1::MatlabStructArray{N},m2::MatlabStructArray{N}) where N
        return isequal(m1.names, m2.names) && isequal(m1.values, m2.values)
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
    function MatlabStructArray(arr::Array{<:AbstractDict{T}, N}) where {T<:AbstractString, N}
        field_names = string.(keys(first(arr)))
        # Ensure same field set for all elements
        for d in arr
            if !issetequal(keys(d), field_names)
                error("Cannot convert Dict array to MatlabStructArray. All elements must share identical field names")
            end
        end
        field_values = Vector{Array{Any,N}}(undef, length(field_names))
        for (idx,f) in enumerate(field_names)
            this_field_values = Array{Any, N}(undef, size(arr))
            for (idx, d) in enumerate(arr)
                this_field_values[idx] = d[f]
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

    function Base.Array(arr::MatlabStructArray{N}) where N
        first_field = first(arr.values)
        sz = size(first_field)
        result = Array{Dict{String,Any}, N}(undef, sz)
        for idx in eachindex(first_field)
            element_values = [v[idx] for v in arr.values]
            result[idx] = Dict{String, Any}(arr.names .=> element_values)
        end
        return result
    end

    struct StructArrayField{N}
        values::Array{Any,N}
    end

    function convert_struct_array(d::Dict{String, Any})
        # there is no possibility of having cell arrays mixed with struct arrays (afaik)
        field_values = first(values(d))
        if field_values isa StructArrayField
            return MatlabStructArray(
                collect(keys(d)),
                [arr.values for arr in values(d)],
            )
        else
            return d
        end
    end

end