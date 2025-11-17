using MAT, Test
using Dates

@testset "MatlabStructArray" begin
    d_arr = Dict{String, Any}[
        Dict("x"=>[1.0,2.0], SubString("y")=>3.0),
        Dict("x"=>[5.0,6.0], "y"=>[])
    ]
    s_arr = MatlabStructArray(d_arr)
    @test s_arr["y"][2] == d_arr[2]["y"]
    @test s_arr["x"][1] == d_arr[1]["x"]

    # constructor errors to protect the user
    @test_throws ErrorException MatlabStructArray(["a", "b"], [[]])
    @test_throws ErrorException MatlabStructArray(["a", "b"], [[],[0.1, 0.2]])

    # equality checks
    @test isequal(MatlabStructArray(["a"], [[0.1, 0.2]]), MatlabStructArray(["a"], [[0.1, 0.2]]))
    @test !isequal(MatlabStructArray(["a"], [[0.1, 0.2]]), MatlabStructArray(["a"], [[0.1, 0.2]], "TestClass"))
    @test !isequal(MatlabStructArray(["a"], [[0.1, 0.2]]), MatlabStructArray(["b"], [[0.1, 0.2]]))
    @test isapprox(MatlabStructArray(["a"], [[0.1, 0.2]]), MatlabStructArray(["a"], [[0.1+eps(0.1), 0.2]]))
    @test !isapprox(MatlabStructArray(["a"], [[0.1, 0.2]]), MatlabStructArray(["b"], [[0.1, 0.2]]))
    @test !isapprox(MatlabStructArray(["a"], [[0.1, 0.2]]), MatlabStructArray(["a"], [[0.11, 0.2]]))

    # empty struct array constructor
    s_arr = MatlabStructArray(["x", "y"], (0,1))
    @test s_arr["x"] == Matrix{Any}(undef, 0, 1)
    @test s_arr["y"] == Matrix{Any}(undef, 0, 1)
    @test MatlabStructArray(["a"]) == MatlabStructArray(["a"], (0,0))

    # convert to Dict to support easy conversion to legacy read behavior
    s_arr = MatlabStructArray(d_arr)
    d = Dict(s_arr)
    @test d isa Dict{String, Any}
    @test collect(keys(d)) == keys(s_arr)
    @test collect(values(d)) == values(s_arr)
    # Dict like interfaces
    @test length(s_arr) == 2
    @test collect(Dict(s_arr)) == collect(s_arr)
    @test haskey(s_arr, "x")
    @test get(s_arr, "x", nothing) == s_arr["x"]
    @test !haskey(s_arr, "wrong")
    @test get(s_arr, "wrong", nothing) === nothing

    # possibility to convert back to dict array via `Array`
    s_arr = MatlabStructArray(d_arr)
    @test Array(s_arr) == d_arr
    d_arr_reshape = reshape(d_arr, 1, 2)
    @test Array(MatlabStructArray(d_arr_reshape)) == d_arr_reshape
    d_symbol = Array{Dict{Symbol,Any}}(MatlabStructArray(d_arr))
    @test d_symbol[2][:x] == d_arr[2]["x"]
    @test Array(MatlabStructArray(d_symbol)) == d_arr

    # class object array conversion
    s_arr = MatlabStructArray(d_arr, "TestClass")
    c_arr = Array(s_arr)
    @test c_arr isa Array{MatlabClassObject}
    @test all(c->c.class=="TestClass", c_arr)
    @test MatlabStructArray(c_arr) == s_arr

    # test error of unequal structs
    wrong_sarr = Dict{String, Any}[
        Dict("x"=>[1.0,2.0], "y"=>[3.0,4.0]),
        Dict("x"=>[5.0,6.0])
    ]
    msg = "Cannot convert Dict array to MatlabStructArray. All elements must share identical field names"
    @test_throws ErrorException(msg) MatlabStructArray(wrong_sarr)
end

@testset "MatlabClassObject" begin
    d = Dict{String,Any}("a" => 5)
    obj = MatlabClassObject(d, "TestClassOld")
    @test keys(obj) == keys(d)
    @test values(obj) == values(d)
    @test collect(obj) == collect(d)
    @test obj["a"] == d["a"]
    @test haskey(obj, "a")
    @test get(obj, "b", "default") == "default"

    obj["b"] = 7
    @test obj["b"] == 7

    c_arr = [MatlabClassObject(d, "TestClassOld"), MatlabClassObject(d, "TestClassOld")]
    s_arr = MatlabStructArray(c_arr)
    @test s_arr.class == "TestClassOld"

    wrong_arr = [MatlabClassObject(d, "TestClassOld"), MatlabClassObject(d, "Bah")]
    @test_throws ErrorException MatlabStructArray(wrong_arr)
end

@testset "MatlabOpaque to_string" begin
    dat = UInt64[
        0x0000000000000001
        0x0000000000000002
        0x0000000000000003
        0x0000000000000001
        0x0000000000000005
        0x0000000000000005
        0x0000000000000005
        0x0065006e006f004a
        0x006f007200420073
        0x006d0053006e0077
        0x0000006800740069
    ]
    dat = reshape(dat, 1, length(str))
    obj = MatlabOpaque(Dict{String, Any}("any" => str), "string")
    str = MAT.MAT_types.to_string(obj)
    @test size(str) == (3,1)
    @test vec(str) == ["Jones", "Brown", "Smith"]

    dat = [
        0x0000000000000001
        0x0000000000000002
        0x0000000000000001
        0x0000000000000001
        0x0000000000000005
        0x0065006e006f004a
        0x0000000000000073
    ]
    dat = reshape(dat, 1, length(str))
    obj = MatlabOpaque(Dict{String, Any}("any" => str), "string")
    str = MAT.MAT_types.to_string(obj)
    @test str == "Jones"
end

@testset "MatlabOpaque to_datetime" begin
    d = Dict{String, Any}(
        "tz"         => "",
        "data"       => ComplexF64[
            1482192000000.0+0.0im;
            1482278400000.0+0.0im;
            1482364800000.0+0.0im;;
            ],
        "fmt"        => "",
        "isDateOnly" => true,
    )
    obj = MatlabOpaque(d, "datetime")
    expected_dates = [
        DateTime(2016, 12, 20) # 20-Dec-2016
        DateTime(2016, 12, 21) # 21-Dec-2016
        DateTime(2016, 12, 22) # 22-Dec-2016
    ]
    @test all(MAT.MAT_types.to_datetime(obj) .== expected_dates)

    d = Dict{String, Any}(
        "tz"         => "",
        "data"       => 1575304969634.0+0.0im,
        "fmt"        => "",
        "isDateOnly" => false,
    )
    obj = MatlabOpaque(d, "datetime")
    # "02-Dec-2019 16:42:49"
    expected_dt = DateTime(2019, 12, 2, 16, 42, 49)
    @test MAT.MAT_types.to_datetime(obj) - expected_dt < Second(1)
end