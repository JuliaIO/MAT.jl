using MAT, Test, SparseArrays
using MAT: MATWriteStyle, MATConstructionError, MatlabStructArray

# Define test structs outside of @testset to avoid world-age issues

struct SUPoint
    x::Float64
    y::Float64
end

struct SUInner
    a::Float64
end
struct SUOuter
    inner::SUInner
    b::Float64
end

@defaults struct SUConfig
    name::String = "default"
    value::Float64 = 0.0
    count::Int = 0
end

# See StructUtils.jl docs for @tags syntax: &(mat=(name="...",),)
# The `mat` key namespaces tags for MAT.jl (via fieldtagkey)
@tags struct SURenamed
    julia_name::Float64 &(mat=(name="matlab_name",),)
    normal::Float64
end

@tags struct SUWithIgnored
    keep::Float64
    skip::Float64 &(mat=(ignore=true,),)
end

abstract type SUShape end
struct SUCircle <: SUShape
    radius::Float64
end
struct SUSquare <: SUShape
    side::Float64
end
@choosetype SUShape x -> haskey(x, "radius") ? SUCircle : SUSquare

struct SUOptional
    name::String
    value::Union{Nothing, Float64}
end

struct SUIntFields
    x::Int
    y::Int
end

# --- Additional test structs ---

struct SUWithArrays
    values::Vector{Float64}
    matrix::Matrix{Float64}
    label::String
end

struct SUWithBool
    flag::Bool
    name::String
end

struct SUWithDict
    metadata::Dict{String, Any}
    id::Float64
end

struct SUWithComplex
    z::ComplexF64
    r::Float64
end

struct SUSymbolField
    name::Symbol
    value::Float64
end
# Symbol needs to be lowered to String for MAT write (MAT can't write raw Symbols).
# Reading back works via the default lift: convert(Symbol, "hello") = :hello
StructUtils.lower(::MATWriteStyle, x::Symbol) = String(x)

@tags struct SUBadRename
    x::Float64 &(mat=(name="1invalid",),)
end

@noarg mutable struct SUNoarg
    x::Float64 = 1.0
    y::Float64 = 2.0
    name::String = "default"
end

@kwarg struct SUKwarg
    x::Float64 = 0.0
    y::Float64 = 0.0
    label::String = "origin"
end

@defaults struct SUInnerDefaults
    a::Float64 = 0.0
    b::String = "inner_default"
end
struct SUOuterWithDefaults
    inner::SUInnerDefaults
    c::Float64
end

abstract type SUAnimal end
struct SUDog <: SUAnimal
    name::String
    breed::String
end
struct SUCat <: SUAnimal
    name::String
    indoor::Bool
end
@choosetype SUAnimal x -> haskey(x, "breed") ? SUDog : SUCat

struct SUHousehold
    owner::String
    pet::SUAnimal
end

struct SUNested
    point::SUPoint
    value::Float64
end

# Matches the struct "s" in test/v7/struct.mat: a=1.0, b=[1.0 2.0], c=[1.0 2.0 3.0]
struct SUMatlabStruct
    a::Float64
    b::Matrix{Float64}
    c::Matrix{Float64}
end

# --- Structs covering the MATLAB/Julia data model mismatches ---

struct SUOneElement
    values::Vector{Float64}
    names::Vector{String}
    zs::Vector{ComplexF64}
end

struct SUMatrixField
    m::Matrix{Float64}
end

struct SUVectorField
    v::Vector{Float64}
end

# A struct array nested inside a struct, the most common MATLAB shape
struct SUStructArrayField
    items::Vector{SUPoint}
    label::String
end

struct SUOptionalMissing
    name::String
    value::Union{Missing, Float64}
end

# Lowered to a plain number on write, to check that `lower` is honored wherever
# a value can appear, not only as a struct field
struct SUTagged
    v::Float64
end
StructUtils.lower(::MATWriteStyle, x::SUTagged) = x.v

struct SUNamed
    name::String
    age::Float64
end

struct SUOptionalString
    name::String
    note::Union{Nothing, String}
end

struct SUOptionalVector
    values::Union{Nothing, Vector{Float64}}
end

# The extension pattern the README documents: a wrapper stored as a plain number
struct SUCelsius
    degrees::Float64
end
StructUtils.lower(::MATWriteStyle, x::SUCelsius) = x.degrees
StructUtils.structlike(::MATStyle, ::Type{SUCelsius}) = false
StructUtils.lift(::Type{SUCelsius}, x::Float64) = SUCelsius(x)

struct SURoom
    temperature::SUCelsius
end

@testset "StructUtils" begin
    mktempdir() do dir
        tmpfile = joinpath(dir, "test.mat")

        @testset "basic struct round-trip" begin
            p = SUPoint(1.0, 2.0)
            matwrite(tmpfile, Dict("p" => p))
            p2 = matread(tmpfile, "p", SUPoint)
            @test p2 == p
            @test p2 isa SUPoint
        end

        @testset "nested structs" begin
            o = SUOuter(SUInner(1.0), 2.0)
            matwrite(tmpfile, Dict("o" => o))
            o2 = matread(tmpfile, "o", SUOuter)
            @test o2 == o
            @test o2 isa SUOuter
            @test o2.inner isa SUInner
        end

        @testset "@defaults" begin
            # Write only some fields via a Dict
            matwrite(tmpfile, Dict("c" => Dict("name" => "test", "value" => 42.0)))
            c = matread(tmpfile, "c", SUConfig)
            @test c.name == "test"
            @test c.value == 42.0
            @test c.count == 0  # default value

            # Full round-trip
            c2 = SUConfig("hello", 3.14, 7)
            matwrite(tmpfile, Dict("c" => c2))
            c3 = matread(tmpfile, "c", SUConfig)
            @test c3 == c2
        end

        @testset "@tags field renaming" begin
            r = SURenamed(1.5, 2.5)
            matwrite(tmpfile, Dict("r" => r))

            # Verify the .mat file has the renamed key
            raw = matread(tmpfile)
            @test haskey(raw["r"], "matlab_name")
            @test !haskey(raw["r"], "julia_name")
            @test haskey(raw["r"], "normal")

            # Read back with type
            r2 = matread(tmpfile, "r", SURenamed)
            @test r2 == r
        end

        @testset "@tags field ignoring" begin
            w = SUWithIgnored(1.0, 2.0)
            matwrite(tmpfile, Dict("w" => w))

            raw = matread(tmpfile)
            @test haskey(raw["w"], "keep")
            @test !haskey(raw["w"], "skip")
            @test raw["w"]["keep"] == 1.0

            # Reading back into SUWithIgnored fails because `skip` has no
            # default and was not written. Types with ignored fields should
            # use @defaults or @noarg if they need to be read back.
            @test_throws MATConstructionError matread(tmpfile, "w", SUWithIgnored)
        end

        @testset "@choosetype abstract dispatch" begin
            matwrite(tmpfile, Dict("c" => SUCircle(5.0)))
            s = matread(tmpfile, "c", SUShape)
            @test s isa SUCircle
            @test s.radius == 5.0

            matwrite(tmpfile, Dict("s" => SUSquare(3.0)))
            s2 = matread(tmpfile, "s", SUShape)
            @test s2 isa SUSquare
            @test s2.side == 3.0
        end

        @testset "Union{T, Nothing}" begin
            o1 = SUOptional("hello", 1.5)
            matwrite(tmpfile, Dict("o" => o1))
            o1b = matread(tmpfile, "o", SUOptional)
            @test o1b.name == "hello"
            @test o1b.value == 1.5
        end

        @testset "numeric type conversion (Float64 to Int)" begin
            # Write a Dict with Float64 values (simulating MATLAB doubles)
            matwrite(tmpfile, Dict("s" => Dict("x" => 3.0, "y" => 7.0)))

            # Raw read gives Float64
            raw = matread(tmpfile)
            @test raw["s"]["x"] isa Float64

            # Typed read converts to Int
            s2 = matread(tmpfile, "s", SUIntFields)
            @test s2 == SUIntFields(3, 7)
            @test s2.x isa Int
            @test s2.y isa Int
        end

        @testset "MatlabStructArray typed read" begin
            sa = MatlabStructArray(["x", "y"], [Any[1.0, 2.0, 3.0], Any[4.0, 5.0, 6.0]])
            matwrite(tmpfile, Dict("points" => sa))

            points = matread(tmpfile, "points", Vector{SUPoint})
            @test length(points) == 3
            @test points[1] == SUPoint(1.0, 4.0)
            @test points[2] == SUPoint(2.0, 5.0)
            @test points[3] == SUPoint(3.0, 6.0)
        end

        @testset "backward compatibility" begin
            p = SUPoint(1.0, 2.0)
            matwrite(tmpfile, Dict("p" => p))

            # matread without type still returns Dict
            result = matread(tmpfile)
            @test result isa Dict{String, Any}
            @test result["p"] isa Dict{String, Any}
            @test result["p"]["x"] == 1.0
        end

        @testset "read(handle, varname, T)" begin
            p = SUPoint(3.0, 4.0)
            matwrite(tmpfile, Dict("p" => p))

            fid = matopen(tmpfile, "r")
            p2 = read(fid, "p", SUPoint)
            close(fid)
            @test p2 == p
        end

        # ---- Struct field types ----

        @testset "struct with array fields" begin
            s = SUWithArrays([1.0, 2.0, 3.0], [1.0 2.0; 3.0 4.0], "test")
            matwrite(tmpfile, Dict("s" => s))
            s2 = matread(tmpfile, "s", SUWithArrays)
            @test s2.values == [1.0, 2.0, 3.0]
            @test s2.matrix == [1.0 2.0; 3.0 4.0]
            @test s2.matrix isa Matrix{Float64}  # verifies array make short-circuit
            @test s2.label == "test"
        end

        @testset "struct with Bool fields" begin
            s = SUWithBool(true, "yes")
            matwrite(tmpfile, Dict("s" => s))
            s2 = matread(tmpfile, "s", SUWithBool)
            @test s2.flag === true
            @test s2.name == "yes"

            s3 = SUWithBool(false, "no")
            matwrite(tmpfile, Dict("s" => s3))
            s4 = matread(tmpfile, "s", SUWithBool)
            @test s4.flag === false
        end

        @testset "struct with Dict field" begin
            inner = Dict{String, Any}("a" => 1.0, "b" => "hello")
            s = SUWithDict(inner, 42.0)
            matwrite(tmpfile, Dict("s" => s))
            s2 = matread(tmpfile, "s", SUWithDict)
            @test s2.metadata isa Dict{String, Any}
            @test s2.metadata["a"] == 1.0
            @test s2.metadata["b"] == "hello"
            @test s2.id == 42.0
        end

        @testset "struct with complex number fields" begin
            s = SUWithComplex(1.0 + 2.0im, 3.0)
            matwrite(tmpfile, Dict("s" => s))
            s2 = matread(tmpfile, "s", SUWithComplex)
            @test s2.z == 1.0 + 2.0im
            @test s2.r == 3.0
        end

        # ---- Write path edge cases ----

        @testset "lower on write (Symbol field)" begin
            s = SUSymbolField(:hello, 1.0)
            matwrite(tmpfile, Dict("s" => s))

            # Verify written as string
            raw = matread(tmpfile)
            @test raw["s"]["name"] == "hello"
            @test raw["s"]["name"] isa String

            # Read back: default lift uses convert(Symbol, "hello") = :hello
            s2 = matread(tmpfile, "s", SUSymbolField)
            @test s2.name === :hello
            @test s2.value == 1.0
        end

        @testset "@tags rename to invalid MATLAB name errors" begin
            @test_throws ErrorException matwrite(tmpfile, Dict("s" => SUBadRename(1.0)))
        end

        @testset "@noarg mutable struct" begin
            s = SUNoarg()
            s.x = 10.0
            s.y = 20.0
            s.name = "test"
            matwrite(tmpfile, Dict("s" => s))
            s2 = matread(tmpfile, "s", SUNoarg)
            @test s2.x == 10.0
            @test s2.y == 20.0
            @test s2.name == "test"

            # With defaults (partial data)
            matwrite(tmpfile, Dict("s" => Dict("x" => 5.0)))
            s3 = matread(tmpfile, "s", SUNoarg)
            @test s3.x == 5.0
            @test s3.y == 2.0   # default
            @test s3.name == "default"  # default
        end

        @testset "@kwarg struct" begin
            s = SUKwarg(x=3.0, y=4.0, label="point")
            matwrite(tmpfile, Dict("s" => s))
            s2 = matread(tmpfile, "s", SUKwarg)
            @test s2.x == 3.0
            @test s2.y == 4.0
            @test s2.label == "point"

            # With defaults (partial data)
            matwrite(tmpfile, Dict("s" => Dict("x" => 1.0)))
            s3 = matread(tmpfile, "s", SUKwarg)
            @test s3.x == 1.0
            @test s3.y == 0.0      # default
            @test s3.label == "origin"  # default
        end

        # ---- Read path edge cases ----

        @testset "extra keys in source ignored" begin
            matwrite(tmpfile, Dict("s" => Dict("x" => 1.0, "y" => 2.0, "z" => 999.0, "extra" => "ignored")))
            p = matread(tmpfile, "s", SUPoint)
            @test p == SUPoint(1.0, 2.0)
        end

        @testset "missing required field errors" begin
            matwrite(tmpfile, Dict("s" => Dict("x" => 1.0)))  # missing "y"
            @test_throws MATConstructionError matread(tmpfile, "s", SUPoint)

            # The error names the type that could not be built and reports the cause
            err = try
                matread(tmpfile, "s", SUPoint)
            catch e
                e
            end
            msg = sprint(showerror, err)
            @test occursin("SUPoint", msg)
            @test occursin("Caused by:", msg)
            @test err.cause isa Exception

            # a missing variable is a file-level error, not a construction failure
            @test_throws KeyError matread(tmpfile, "no_such_variable", SUPoint)
        end

        @testset "multiple typed reads from same file" begin
            matwrite(tmpfile, Dict(
                "point" => SUPoint(1.0, 2.0),
                "config" => SUConfig("test", 3.14, 5)
            ))

            p = matread(tmpfile, "point", SUPoint)
            @test p == SUPoint(1.0, 2.0)

            c = matread(tmpfile, "config", SUConfig)
            @test c == SUConfig("test", 3.14, 5)

            # Also via matopen handle
            fid = matopen(tmpfile, "r")
            p2 = read(fid, "point", SUPoint)
            c2 = read(fid, "config", SUConfig)
            close(fid)
            @test p2 == SUPoint(1.0, 2.0)
            @test c2 == SUConfig("test", 3.14, 5)
        end

        @testset "nested struct with @defaults on inner" begin
            # Write with missing inner fields
            matwrite(tmpfile, Dict("s" => Dict("inner" => Dict("a" => 5.0), "c" => 10.0)))
            s = matread(tmpfile, "s", SUOuterWithDefaults)
            @test s.inner.a == 5.0
            @test s.inner.b == "inner_default"  # default from SUInnerDefaults
            @test s.c == 10.0
        end

        @testset "@choosetype inside parent struct" begin
            h = SUHousehold("Alice", SUDog("Rex", "Labrador"))
            matwrite(tmpfile, Dict("h" => h))
            h2 = matread(tmpfile, "h", SUHousehold)
            @test h2.owner == "Alice"
            @test h2.pet isa SUDog
            @test h2.pet.name == "Rex"
            @test h2.pet.breed == "Labrador"

            h3 = SUHousehold("Bob", SUCat("Whiskers", true))
            matwrite(tmpfile, Dict("h" => h3))
            h4 = matread(tmpfile, "h", SUHousehold)
            @test h4.owner == "Bob"
            @test h4.pet isa SUCat
            @test h4.pet.name == "Whiskers"
            @test h4.pet.indoor == true
        end

        @testset "MatlabStructArray with nested struct data" begin
            sa = MatlabStructArray(
                ["point", "value"],
                [
                    Any[Dict("x" => 1.0, "y" => 2.0), Dict("x" => 3.0, "y" => 4.0)],
                    Any[10.0, 20.0]
                ]
            )
            matwrite(tmpfile, Dict("arr" => sa))
            arr = matread(tmpfile, "arr", Vector{SUNested})
            @test length(arr) == 2
            @test arr[1].point == SUPoint(1.0, 2.0)
            @test arr[1].value == 10.0
            @test arr[2].point == SUPoint(3.0, 4.0)
            @test arr[2].value == 20.0
        end

        # ---- MATLAB/Julia data model conversions ----

        @testset "one-element array fields" begin
            # MATLAB has no scalars, so a one-element array is written as 1x1 and
            # read back as a scalar; it must still fill an array-typed field.
            s = SUOneElement([1.0], ["a"], [1.0 + 2.0im])
            matwrite(tmpfile, Dict("s" => s))
            s2 = matread(tmpfile, "s", SUOneElement)
            @test s2.values == [1.0]
            @test s2.names == ["a"]
            @test s2.zs == [1.0 + 2.0im]

            s3 = SUOneElement([1.0, 2.0], ["ab", "cd"], [1.0im, 2.0im])
            matwrite(tmpfile, Dict("s" => s3))
            s4 = matread(tmpfile, "s", SUOneElement)
            @test s4.values == [1.0, 2.0]
            @test s4.names == ["ab", "cd"]
            @test s4.zs == [1.0im, 2.0im]
        end

        @testset "1x1 matrix field" begin
            m = reshape([1.0], 1, 1)
            matwrite(tmpfile, Dict("s" => SUMatrixField(m)))
            s = matread(tmpfile, "s", SUMatrixField)
            @test s.m == m
            @test s.m isa Matrix{Float64}
        end

        @testset "array element type conversion" begin
            for src in (Float32[1 2; 3 4], Int32[1 2; 3 4], Any[1.0 2.0; 3.0 4.0],
                        [true false; false true])
                matwrite(tmpfile, Dict("s" => Dict("m" => src)))
                s = matread(tmpfile, "s", SUMatrixField)
                @test s.m isa Matrix{Float64}
                @test s.m == Float64.(src)
            end
        end

        @testset "array shape adaptation" begin
            # a MATLAB row or column vector reads as a Vector
            for src in ([1.0 2.0 3.0], reshape([1.0, 2.0, 3.0], 3, 1))
                matwrite(tmpfile, Dict("s" => Dict("v" => src)))
                @test matread(tmpfile, "s", SUVectorField).v == [1.0, 2.0, 3.0]
            end

            # a vector source fills a matrix field by padding a singleton dimension
            matwrite(tmpfile, Dict("s" => Dict("m" => [1.0, 2.0, 3.0])))
            @test matread(tmpfile, "s", SUMatrixField).m == reshape([1.0, 2.0, 3.0], 3, 1)

            # a genuine shape mismatch is reported rather than silently flattened
            matwrite(tmpfile, Dict("s" => Dict("v" => [1.0 2.0; 3.0 4.0])))
            @test_throws MATConstructionError matread(tmpfile, "s", SUVectorField)
        end

        @testset "empty array fields" begin
            matwrite(tmpfile, Dict("s" => SUVectorField(Float64[])))
            v = matread(tmpfile, "s", SUVectorField).v
            @test v isa Vector{Float64}
            @test isempty(v)
        end

        @testset "array target types" begin
            matwrite(tmpfile, Dict("v" => [1.0, 2.0, 3.0], "m" => [1.0 2.0; 3.0 4.0]))

            # dimensionality known but element type open
            @test matread(tmpfile, "m", Matrix) == [1.0 2.0; 3.0 4.0]
            @test matread(tmpfile, "v", Vector) == [1.0, 2.0, 3.0]
            @test matread(tmpfile, "v", Matrix) == reshape([1.0, 2.0, 3.0], 3, 1)
            @test_throws MATConstructionError matread(tmpfile, "m", Vector)

            # element type known but dimensionality open
            @test matread(tmpfile, "m", Array{Float64}) == [1.0 2.0; 3.0 4.0]
            @test matread(tmpfile, "v", AbstractVector{Float64}) == [1.0, 2.0, 3.0]

            # 0-dimensional targets hold exactly one element
            matwrite(tmpfile, Dict("s" => 5.0, "m" => [1.0 2.0; 3.0 4.0]))
            zerodim = matread(tmpfile, "s", Array{Float64,0})
            @test zerodim isa Array{Float64,0}
            @test zerodim[] == 5.0
            @test_throws MATConstructionError matread(tmpfile, "m", Array{Float64,0})

            # array types that are not `Array`
            matwrite(tmpfile, Dict("m" => [1.0 0.0; 0.0 4.0], "b" => [true false]))
            @test matread(tmpfile, "m", SparseMatrixCSC{Float64,Int}) ==
                  sparse([1.0 0.0; 0.0 4.0])
            @test matread(tmpfile, "b", BitMatrix) == BitMatrix([true false])

            # a stored value that already satisfies a union of array types is kept,
            # but one that would have to be converted cannot be disambiguated
            matwrite(tmpfile, Dict("m" => [1.0 2.0; 3.0 4.0], "f" => Float32[1 2; 3 4]))
            U = Union{Vector{Float64},Matrix{Float64}}
            @test matread(tmpfile, "m", U) == [1.0 2.0; 3.0 4.0]
            @test_throws MATConstructionError matread(tmpfile, "f", U)
        end

        @testset "collection targets from a struct array" begin
            sa = MatlabStructArray(["x", "y"], [Any[1.0, 2.0], Any[3.0, 4.0]])
            matwrite(tmpfile, Dict("a" => sa))
            pts = [SUPoint(1.0, 3.0), SUPoint(2.0, 4.0)]

            @test matread(tmpfile, "a", Set{SUPoint}) == Set(pts)
            @test matread(tmpfile, "a", Tuple{SUPoint,SUPoint}) == (pts[1], pts[2])
            @test matread(tmpfile, "a", Union{Nothing,Vector{SUPoint}}) == pts
            @test matread(tmpfile, "a", Union{Missing,Vector{SUPoint}}) == pts

            # a struct target still requires a single element
            @test_throws MATConstructionError matread(tmpfile, "a", SUPoint)
        end

        @testset "struct array as a struct field" begin
            sa = MatlabStructArray(["x", "y"], [Any[1.0, 2.0], Any[3.0, 4.0]])
            matwrite(tmpfile, Dict("s" => Dict("items" => sa, "label" => "L")))
            s = matread(tmpfile, "s", SUStructArrayField)
            @test s.label == "L"
            @test s.items == [SUPoint(1.0, 3.0), SUPoint(2.0, 4.0)]
        end

        @testset "typed read returns the requested type" begin
            sa = MatlabStructArray(["x", "y"], [Any[1.0, 2.0], Any[3.0, 4.0]])
            matwrite(tmpfile, Dict("a" => sa))
            @test matread(tmpfile, "a", Vector{SUPoint}) isa Vector{SUPoint}
            @test matread(tmpfile, "a", Matrix{SUPoint}) isa Matrix{SUPoint}

            # a one-element struct array is stored as a plain 1x1 struct, and can
            # be read back either as a single struct or as a one-element array
            one = MatlabStructArray(["x", "y"], [Any[1.0], Any[2.0]])
            matwrite(tmpfile, Dict("a" => one))
            @test matread(tmpfile, "a", SUPoint) == SUPoint(1.0, 2.0)
            @test matread(tmpfile, "a", Vector{SUPoint}) == [SUPoint(1.0, 2.0)]
        end

        @testset "nothing and missing round-trip via the empty matrix" begin
            matwrite(tmpfile, Dict("o" => SUOptional("hi", nothing)))
            @test size(matread(tmpfile)["o"]["value"]) == (0, 0)  # MATLAB `[]`
            @test matread(tmpfile, "o", SUOptional) == SUOptional("hi", nothing)

            matwrite(tmpfile, Dict("o" => SUOptionalMissing("hi", missing)))
            o = matread(tmpfile, "o", SUOptionalMissing)
            @test o.name == "hi"
            @test o.value === missing

            matwrite(tmpfile, Dict("o" => SUOptionalMissing("hi", 2.0)))
            @test matread(tmpfile, "o", SUOptionalMissing).value == 2.0
        end

        @testset "lower applies wherever a value appears" begin
            matwrite(tmpfile, Dict(
                "top" => SUTagged(1.0),
                "dict" => Dict("k" => SUTagged(2.0)),
                "cell" => Any[SUTagged(3.0), SUTagged(4.0)],
            ))
            raw = matread(tmpfile)
            @test raw["top"] == 1.0
            @test raw["dict"]["k"] == 2.0
            @test raw["cell"] == Any[3.0, 4.0]
        end

        @testset "write path type handling" begin
            # NamedTuples are written as structs, including the columnar shape
            # that Tables.jl would report as a table
            matwrite(tmpfile, Dict("nt" => (a=1.0, b="hi", c=[1.0, 2.0])))
            nt = matread(tmpfile)["nt"]
            @test nt["a"] == 1.0
            @test nt["b"] == "hi"
            @test nt["c"] == [1.0, 2.0]

            matwrite(tmpfile, Dict("nt" => (a=[1.0, 2.0], b=[3.0, 4.0])))
            nt2 = matread(tmpfile)["nt"]
            @test nt2["a"] == [1.0, 2.0]
            @test nt2["b"] == [3.0, 4.0]

            # a type with no MATLAB representation is reported by name
            @test_throws ErrorException("cannot write a value of type `$(typeof(sin))` to a MAT file") matwrite(tmpfile, Dict("f" => sin))
        end

        @testset "empty strings round-trip" begin
            matwrite(tmpfile, Dict("p" => SUNamed("", 30.0)))
            @test matread(tmpfile)["p"]["name"] == String[]   # how MATLAB stores `''`
            @test matread(tmpfile, "p", SUNamed) == SUNamed("", 30.0)

            # an empty string is a value, so it must not narrow to `nothing`,
            # while an absent value written as `[]` still must
            matwrite(tmpfile, Dict("o" => SUOptionalString("n", "")))
            @test matread(tmpfile, "o", SUOptionalString) == SUOptionalString("n", "")
            matwrite(tmpfile, Dict("o" => SUOptionalString("n", nothing)))
            @test matread(tmpfile, "o", SUOptionalString) == SUOptionalString("n", nothing)
            matwrite(tmpfile, Dict("o" => SUOptionalString("n", "hi")))
            @test matread(tmpfile, "o", SUOptionalString).note == "hi"

            # a non-empty array is still not a string
            matwrite(tmpfile, Dict("p" => Dict("name" => [1.0, 2.0], "age" => 30.0)))
            @test_throws MATConstructionError matread(tmpfile, "p", SUNamed)

            # nor is an empty numeric one: only empty character data reads as `""`
            matwrite(tmpfile, Dict("p" => Dict("name" => zeros(0, 0), "age" => 30.0)))
            @test_throws MATConstructionError matread(tmpfile, "p", SUNamed)
        end

        @testset "an empty numeric array narrows to nothing" begin
            # the inherent ambiguity the docs call out: an empty `Vector{Float64}`
            # is written as `[]` and cannot be told apart from a missing value
            matwrite(tmpfile, Dict("o" => SUOptionalVector(Float64[])))
            @test matread(tmpfile, "o", SUOptionalVector).values === nothing
            matwrite(tmpfile, Dict("o" => SUOptionalVector([1.0, 2.0])))
            @test matread(tmpfile, "o", SUOptionalVector).values == [1.0, 2.0]
        end

        @testset "@choosetype against a struct array source" begin
            one = MatlabStructArray(["radius"], [Any[1.0]])
            matwrite(tmpfile, Dict("s" => one))
            @test matread(tmpfile, "s", SUShape) == SUCircle(1.0)

            two = MatlabStructArray(["side"], [Any[1.0, 2.0]])
            matwrite(tmpfile, Dict("s" => two))
            @test matread(tmpfile, "s", Vector{SUShape}) == [SUSquare(1.0), SUSquare(2.0)]

            # a single target still needs a single element, and says so
            err = try
                matread(tmpfile, "s", SUShape)
            catch e
                e
            end
            @test err isa MATConstructionError
            @test occursin("2-element MATLAB struct array", sprint(showerror, err))
        end

        @testset "@choosetype field holding a struct array" begin
            # a one-element struct array stored in a field reaches the chooser as a
            # `MatlabStructArray` rather than a `Dict`; this checks that a chooser
            # testing with `haskey`, as the docs advise, still picks the right type
            pet = MatlabStructArray(["name", "breed"], [Any["Rex"], Any["Labrador"]])
            matwrite(tmpfile, Dict("h" => Dict("owner" => "Alice", "pet" => pet)))
            @test matread(tmpfile, "h", SUHousehold) ==
                  SUHousehold("Alice", SUDog("Rex", "Labrador"))

            # cell array elements are unwrapped before dispatch, like a whole variable
            dog(name) = MatlabStructArray(["name", "breed"], [Any[name], Any["Labrador"]])
            matwrite(tmpfile, Dict("c" => Any[dog("Rex"), dog("Max")]))
            @test matread(tmpfile, "c", Vector{SUAnimal}) ==
                  [SUDog("Rex", "Labrador"), SUDog("Max", "Labrador")]
        end

        @testset "lower/lift extension pattern" begin
            matwrite(tmpfile, Dict("r" => SURoom(SUCelsius(21.5))))
            @test matread(tmpfile)["r"]["temperature"] == 21.5   # stored as a number
            @test matread(tmpfile, "r", SURoom) == SURoom(SUCelsius(21.5))
        end

        # ---- Format coverage ----

        @testset "compressed round-trip" begin
            p = SUPoint(1.0, 2.0)
            matwrite(tmpfile, Dict("p" => p); compress=true)
            p2 = matread(tmpfile, "p", SUPoint)
            @test p2 == p
        end

        @testset "v4 format typed read" begin
            # v4 files hold only numeric matrices, so there are no structs to build
            v4file = joinpath(dir, "v4.mat")
            matwrite(v4file, Dict("m" => [1.0 2.0; 3.0 4.0], "v" => [1.0 2.0 3.0]);
                     version="v4")
            @test matread(v4file, "m", Matrix{Float64}) == [1.0 2.0; 3.0 4.0]
            @test matread(v4file, "v", Vector{Float64}) == [1.0, 2.0, 3.0]
            @test matread(v4file, "m", Matrix{Float32}) isa Matrix{Float32}

            # the v4 writer lowers values and writes `nothing` as `[]`, as v7.3 does
            matwrite(v4file, Dict("t" => SUTagged(1.0), "n" => nothing); version="v4")
            raw = matread(v4file)
            @test only(raw["t"]) == 1.0   # v4 reads scalars back as 1x1 matrices
            @test isempty(raw["n"])

            # a 1x1 matrix still fills a scalar target, as a collapsed scalar does
            # in v5 and v7.3, but a larger array does not
            @test matread(v4file, "t", Float64) === 1.0
            @test matread(v4file, "t", Int) === 1
            matwrite(v4file, Dict("v" => [1.0 2.0]); version="v4")
            @test_throws MATConstructionError matread(v4file, "v", Float64)
            # a one-element cell array is not a number
            matwrite(tmpfile, Dict("c" => Any[1.0]))
            @test_throws MATConstructionError matread(tmpfile, "c", Float64)

            # v4 stores `''` and `[]` identically, so an empty string narrows to
            # `nothing` like `[]` does, yet still reads as `""` for a string target
            matwrite(v4file, Dict("e" => ""); version="v4")
            @test matread(v4file, "e", Union{Nothing,String}) === nothing
            @test matread(v4file, "e", String) == ""

            # the typed read is also reachable from a v4 file handle
            matwrite(v4file, Dict("m" => [1.0 2.0; 3.0 4.0]); version="v4")
            fid = matopen(v4file)
            try
                @test read(fid, "m", Matrix{Float64}) == [1.0 2.0; 3.0 4.0]
            finally
                close(fid)
            end
        end
    end

    # ---- v5 format coverage (outside mktempdir since we read existing files) ----

    @testset "v5 format typed read" begin
        # test/v7/struct.mat is a v5-format file with struct "s" containing
        # a=1.0, b=[1.0 2.0], c=[1.0 2.0 3.0]
        v7dir = joinpath(dirname(@__FILE__), "v7")
        s = matread(joinpath(v7dir, "struct.mat"), "s", SUMatlabStruct)
        @test s isa SUMatlabStruct
        @test s.a == 1.0
        @test s.b == [1.0 2.0]
        @test s.c == [1.0 2.0 3.0]

        # and from a v5 file handle
        fid = matopen(joinpath(v7dir, "struct.mat"))
        try
            s2 = read(fid, "s", SUMatlabStruct)
            @test s2 isa SUMatlabStruct
            @test (s2.a, s2.b, s2.c) == (s.a, s.b, s.c)
        finally
            close(fid)
        end
    end
end
