# function_handle.jl
# FunctionHandle: MATLAB function handles, stored as structs

"""
    FunctionHandle(d::Dict{String, Any}) <: AbstractDict{String, Any}

Type to store function handles which are stored as structs in MATLAB.
"""
struct FunctionHandle{D<:AbstractDict{String,Any}} <: AbstractDict{String,Any}
    d::D
end

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
