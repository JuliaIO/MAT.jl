using MAT, Test, Dates
using SparseArrays, LinearAlgebra

tmpfile = string(tempname(), ".mat")

function test_write(data; kwargs...)
    matwrite(tmpfile, data; kwargs...)

    fid = matopen(tmpfile, "r")
    local result
    try
        result = read(fid)
    finally
        close(fid)
    end

    @test isequal(result, data)
end

function test_write(data)
    test_write(data; compress = false)
    test_write(data; compress = true)
end

function test_compression_effective(data)
    test_write(data; compress = false)
    sizeUncompressed = stat(tmpfile).size
    test_write(data; compress = true)
    sizeCompressed = stat(tmpfile).size

    if sizeCompressed >= sizeUncompressed
        error("Compression was not effective")
    end
end

@testset "write error messages" begin
    msg = "writing for \"v7\" is not supported"
    @test_throws ErrorException(msg) matwrite(tmpfile, Dict("s" => 1); version="v7")

    msg = "matwrite requires a Dict with ASCII keys"
    @test_throws ErrorException(msg) matwrite(tmpfile, Dict(:s => 1))
    @test_throws ErrorException(msg) matwrite(tmpfile, Dict(:s => 1); version="v4")
end

test_write(Dict(
    "int8" => Int8(1),
    "uint8" => UInt8(1),
    "int16" => Int16(1),
    "uint16" => UInt16(1),
    "int32" => Int32(1),
    "uint32" => UInt32(1),
    "int64" => Int64(1),
    "uint64" => UInt64(1),
    "single" => Float32(1),
    "double" => Float64(1),
    "logical" => true
))

test_write(Dict(
    "ComplexInt" => Complex{Int}[1 -1 1+1im 1-1im -1+1im -1-1im 1im],
    "ComplexF32" => ComplexF32[1.0 -1.0 1.0+1.0im 1.0-1.0im -1.0+1.0im -1.0-1.0im 1.0im],
    "ComplexF64" => [1.0 -1.0 1.0+1.0im 1.0-1.0im -1.0+1.0im -1.0-1.0im 1.0im],
    "ComplexPair" => [1 2-3im 4+5im]
))
test_write(Dict("ComplexF64" => 1.0im, "ComplexPair" => 2-3im))

test_write(Dict(
    "simple_string" => "the quick brown fox",
    "accented_string" => "thé qüîck browñ fòx",
    "concatenated_strings" => ["this is a string", "this is another string"],
    "cell_strings" => ["this is a string" "this is another string"],
    "empty_string" => ""
))

test_write(Dict(
    "a1x2" => [1.0 2.0],
    "a2x1" => zeros(2, 1)+[1.0, 2.0],
    "a2x2" => [1.0 3.0; 4.0 2.0],
    "a2x2x2" => cat([1.0 3.0; 4.0 2.0], [1.0 2.0; 3.0 4.0], dims=3),
    "empty" => zeros(0, 0),
    "bitarray" => trues(3, 3),
    "string" => "string"
))

# cannot distinguish char from single element string
test_write(Dict("char" => 'a'))
# inconsistent behavior in v4
matwrite(tmpfile, Dict("char" => 'a'), version="v4")
@test matread(tmpfile)["char"] == "a"

test_write(Dict(
    "cell" => Any[1 2.01 "string" Any["string1" "string2"]]
))

s3 = Dict()
for i in 1:526
    s3["field$i"] = i
end

test_write(Dict(
    "s" => Dict(
        "a" => 1.0,
        "b" => [1.0 2.0],
        "c" => [1.0 2.0 3.0]
    ),
    "s2" => Dict("a" => [1.0 2.0]),
    "s3" => s3
))

test_write(Dict(
    "sparse_empty" => sparse(Matrix{Float64}(undef, 0, 0)),
    "sparse_eye" => sparse(1.0I, 20, 20),
    "sparse_logical" => SparseMatrixCSC{Bool,Int}(5, 5, [1:6;], [1:5;], fill(true, 5)),
    "sparse_random" => sparse([0 6. 0; 8. 0 1.; 0 0 9.]),
    "sparse_complex" => sparse([0 6. 0; 8. 0 1.; 0 0 9.]*(1. + 1.0im)),
    "sparse_zeros" => SparseMatrixCSC(20, 20, ones(Int, 21), Int[], Float64[])
))

@test_throws ErrorException test_write(Dict("1invalidkey" => "starts with a number"))
@test_throws ErrorException test_write(Dict("another invalid key" => "invalid characters"))
@test_throws ErrorException test_write(Dict("yetanotherinvalidkeyyetanotherinvalidkeyyetanotherinvalidkeyyetanotherinvalidkey" => "too long"))

struct TestCompositeKind
    field1::AbstractString
end
fid = matopen(tmpfile, "w")
write(fid, "test", TestCompositeKind("test value"))
close(fid)
fid = matopen(tmpfile, "r")
result = read(fid, "test")
close(fid)
@assert result == Dict("field1" => "test value")


fid = matopen(tmpfile, "w")
@test_throws ErrorException write(fid, "1invalidvarname", "1invalidvarvalue")
close(fid)

using DataStructures
sd = SortedDict(Dict(
    "uint16" => UInt16(1),
    "ComplexF64" => [1.0 -1.0 1.0+1.0im 1.0-1.0im -1.0+1.0im -1.0-1.0im 1.0im],
    "simple_string" => "the quick brown fox",
    "a1x2" => [1.0 2.0],
    "sparse_empty" => sparse(Matrix{Float64}(undef, 0, 0))
))
test_write(sd)

# note: compression is NOT effective when the dict contains many duplicate entries
# which are not compressible by themselves!
test_compression_effective(Dict("data" => fill(1.0, 1000)))

# test adjoint/reshape array
test_write(Dict("adjoint_arr"=>[1 2 3;4 5 6;7 8 9]'))
test_write(Dict("reshape_arr"=>reshape([1 2 3;4 5 6;7 8 9]',1,9)))

test_write(Dict("adjoint_arr"=>Any[1 2 3;4 5 6;7 8 9]'))
test_write(Dict("reshape_arr"=>reshape(Any[1 2 3;4 5 6;7 8 9]',1,9)))

@testset "named tuple" begin
    nt = (x = 5, y = Any[6, "string"])
    matwrite(tmpfile, Dict("nt" => nt))
    nt_read = matread(tmpfile)["nt"]
    @test nt_read["x"] == 5
    @test nt_read["y"] == nt.y
end

@testset "tuple" begin
    # NTuple{T}
    t = (5, 6)
    matwrite(tmpfile, Dict("t" => (5, 6)))
    r = matread(tmpfile)["t"]
    @test r == [x for x in t]

    # otherwise cell array
    t = (5, "string")
    matwrite(tmpfile, Dict("t" => t))
    r = matread(tmpfile)["t"]
    @test r == [x for x in t]
end

@testset "symbol" begin
    matwrite(tmpfile, Dict("s" => :symbol))
    r = matread(tmpfile)["s"]
    @test r == "symbol"
end

# test nested struct array - interface via Dict array
@testset "MatlabStructArray writing" begin
    sarr = Dict{String, Any}[
        Dict("x"=>[1.0,2.0], SubString("y")=>3.0),
        Dict("x"=>[5.0,6.0], "y"=>[Dict("a"=>7), Dict("a"=>8)])
    ]
    # we have to test Array size is maintained inside mat files
    sarr = reshape(sarr, 1, 2)
    matwrite(tmpfile, Dict("s_array" => sarr))
    read_sarr = matread(tmpfile)["s_array"]
    @test read_sarr isa MatlabStructArray
    @test read_sarr["y"][2] isa MatlabStructArray

    sarr = Dict{String, Any}[
        Dict("x"=>[1.0,2.0], SubString("y")=>3.0),
        Dict("x"=>[5.0,6.0], "y"=>[])
    ]
    test_write(Dict("s_array" => MatlabStructArray(sarr)))

    empty_sarr = MatlabStructArray(["a", "b", "c"])
    test_write(Dict("s_array" => empty_sarr))

    # old matlab class object array
    carr = MatlabStructArray(["foo"], [[5, "bar"]], "TestClassOld")
    test_write(Dict("class_array" => carr))

    d = Dict{String,Any}("foo" => 5)
    obj = MatlabClassObject(d, "TestClassOld")
    test_write(Dict("tc_old" => obj))

    carr = [MatlabClassObject(d, "TestClassOld"), MatlabClassObject(d, "TestClassOld")]
    matwrite(tmpfile, Dict("class_array" => carr))
    carr_read = matread(tmpfile)["class_array"]
    @test carr_read == MatlabStructArray(carr)

    s_large1 = Dict()
    s_large2 = Dict()
    for i in 1:600
        s_large1["field$i"] = i
        s_large2["field$i"] = i + 1000
    end
    sarr = Dict{String, Any}[s_large1, s_large2]
    test_write(Dict("s_array_large" => MatlabStructArray(sarr)))
end

@testset "MatlabOpaque simple" begin
    d = Dict{String,Any}("a" => 1, "b" => Any[1.0, 2.0])
    obj = MatlabOpaque(d, "TestClass")
    var_dict = Dict("var" => obj)

    mktempdir() do tmpdir
        tmpfile = joinpath(tmpdir, "test.mat")
        matwrite(tmpfile, var_dict)
        read_var = matread(tmpfile)

        @test haskey(read_var, "var")
        @test isa(read_var["var"], MatlabOpaque)
        @test read_var["var"].class == obj.class

        @test haskey(read_var["var"], "a")
        @test haskey(read_var["var"], "b")
        @test read_var["var"]["a"] == obj["a"]
        @test isequal(read_var["var"]["b"], obj["b"])
    end
end

@testset "Empty Struct 1x1" begin
    var_dict = Dict("empty_struct" => Dict{String,Any}())
    mktempdir() do tmpdir
        tmpfile = joinpath(tmpdir, "test.mat")
        matwrite(tmpfile, var_dict)
        read_var = matread(tmpfile)

        @test haskey(read_var, "empty_struct")
        @test isa(read_var["empty_struct"], Dict{String,Any})
        @test length(keys(read_var["empty_struct"])) == 0
    end
end

@testset "MatlabOpaque handle" begin
    d = Dict{String,Any}("a" => 1, "b" => Any[1.0, 2.0])
    obj = MatlabOpaque(d, "TestClassHandle")
    var_dict = Dict("var1" => obj, "var2" => obj)

    mktempdir() do tmpdir
        tmpfile = joinpath(tmpdir, "test.mat")
        matwrite(tmpfile, var_dict)
        read_var = matread(tmpfile)

        @test haskey(read_var, "var1")
        @test haskey(read_var, "var2")
        @test isa(read_var["var1"], MatlabOpaque)
        @test isa(read_var["var2"], MatlabOpaque)
        @test read_var["var1"] === read_var["var2"]  # same object
    end
end

@testset "MatlabOpaque Array" begin
    d1 = Dict{String, Any}("a" => 1)
    d2 = Dict{String, Any}("a" => 2)
    obj1 = MatlabOpaque(d1, "TestClassArray")
    obj2 = MatlabOpaque(d2, "TestClassArray")
    obj_arr = Array{MatlabOpaque}(undef, 2, 1)
    obj_arr[1, 1] = obj1
    obj_arr[2, 1] = obj2
    var_dict = Dict("obj_array" => obj_arr)

    mktempdir() do tmpdir
        tmpfile = joinpath(tmpdir, "test.mat")
        matwrite(tmpfile, var_dict)
        read_var = matread(tmpfile)

        @test haskey(read_var, "obj_array")
        @test isa(read_var["obj_array"], Array{MatlabOpaque})
        @test size(read_var["obj_array"]) == (2, 1)
        @test read_var["obj_array"][1, 1]["a"] == 1
        @test read_var["obj_array"][2, 1]["a"] == 2
    end
end

@testset "MatlabOpaque Nested object" begin
    inner_dict = Dict{String,Any}("a" => 1, "b" => 2)
    inner_obj = MatlabOpaque(inner_dict, "InnerClass")
    outer_dict = Dict{String,Any}("inner" => inner_obj, "c" => 3)
    outer_obj = MatlabOpaque(outer_dict, "OuterClass")
    var_dict = Dict("outer_obj" => outer_obj)

    mktempdir() do tmpdir
        tmpfile = joinpath(tmpdir, "test.mat")
        matwrite(tmpfile, var_dict)
        read_var = matread(tmpfile)

        @test haskey(read_var, "outer_obj")
        @test isa(read_var["outer_obj"], MatlabOpaque)
        @test read_var["outer_obj"].class == outer_obj.class

        @test haskey(read_var["outer_obj"], "inner")
        @test haskey(read_var["outer_obj"], "c")
        @test read_var["outer_obj"]["c"] == outer_obj["c"]

        inner_read = read_var["outer_obj"]["inner"]
        @test isa(inner_read, MatlabOpaque)
        @test inner_read.class == inner_obj.class

        @test haskey(inner_read, "a")
        @test haskey(inner_read, "b")
        @test inner_read["a"] == inner_obj["a"]
        @test inner_read["b"] == inner_obj["b"]
    end

end

@testset "Dates" begin
    dt = DateTime[
        DateTime(2016, 12, 20) # 20-Dec-2016
        DateTime(2016, 12, 21) # 21-Dec-2016
        DateTime(2016, 12, 22) # 22-Dec-2016
    ]
    test_write(Dict{String,Any}("dt" => dt))
    test_write(Dict{String,Any}("dt" => dt[1]))

    ms = Millisecond(500)
    test_write(Dict{String,Any}("ms" => ms))
    test_write(Dict{String,Any}("ms" => [ms, Millisecond(50000)]))
end
