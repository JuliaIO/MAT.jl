# opaque_categorical.jl
# Conversion of the MATLAB `categorical` class

function from_categorical(obj::MatlabOpaque)
    category_names = obj["categoryNames"]
    codes = obj["codes"]
    pool = vec(Array{promoted_eltype(category_names)}(category_names))
    if any(iszero, codes)
        # an <undefined> element has the code 0, which a PooledArray cannot hold; it
        # keeps `missing` as a pool entry instead, which may need a wider code type
        pool = Union{Missing,eltype(pool)}[pool; missing]
        n = length(pool)
        C = n <= typemax(eltype(codes)) ? eltype(codes) : widen(eltype(codes))
        codes = map(c -> iszero(c) ? C(n) : C(c), codes)
    end
    code_type = eltype(codes)
    pool_type = eltype(pool)
    invpool = Dict{pool_type,code_type}(pool .=> code_type.(1:length(pool)))
    RA = typeof(codes)
    N = ndims(codes)
    return PooledArray{pool_type,code_type,N,RA}(RefArray(codes), invpool, pool)
end

function promoted_eltype(v::AbstractArray{Any})
    # a categorical with no categories, such as `categorical({})`
    isempty(v) && return Any
    first_el, remaining = Iterators.peel(v)
    T_out = typeof(first_el)
    for el in remaining
        T_out = promote_type(T_out, typeof(el))
    end
    return T_out
end
promoted_eltype(::AbstractArray{T}) where {T} = T
