# opaque.jl
# MatlabOpaque: modern MATLAB class objects, and their conversion

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
    elseif obj.class == "timetable"
        return from_timetable(obj, table)
    else
        return obj
    end
end
