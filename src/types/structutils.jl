# structutils.jl
# StructUtils styles, conversion rules and typed construction from MAT data

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

