# utils.jl
# Small shared helpers used across the MAT.jl types

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
value_isapprox(d1::AbstractDict, d2::AbstractDict; kwargs...) = dict_isapprox(d1, d2; kwargs...)

map_or_not(f, dat::AbstractArray) = map(f, dat)
map_or_not(f, dat) = f(dat)
