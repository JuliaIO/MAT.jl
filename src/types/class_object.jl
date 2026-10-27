# class_object.jl
# MatlabClassObject: old-style MATLAB class objects

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
