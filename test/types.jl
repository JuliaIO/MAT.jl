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
    # name order shouldn't matter in isapprox
    @test isapprox(MatlabStructArray(["a", "b"], [[0.1], [0.2]]), MatlabStructArray(["b", "a"], [[0.2], [0.1]]))

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

    # setindexing
    msg = "field \"z\" not found in MatlabStructArray"
    @test_throws ErrorException(msg) s_arr["z"] = [1.0, 2.0]
    s_arr["x"] = [1.0, 2.0]
    @test s_arr["x"] == Any[1.0, 2.0]
    msg = "size of input is not equal to MatlabStructArray column size"
    @test_throws ErrorException(msg) s_arr["x"] = Any[1.0]
    @test_throws ErrorException(msg) s_arr["x"] = reshape(Any[1.0,2.0], 1, 2)

    # possibility to convert back to dict array via `Array`
    s_arr = MatlabStructArray(d_arr)
    @test Array(s_arr) == d_arr
    d_arr_reshape = reshape(d_arr, 1, 2)
    @test Array(MatlabStructArray(d_arr_reshape)) == d_arr_reshape
    d_symbol = Array{Dict{Symbol,Any}}(MatlabStructArray(d_arr))
    @test d_symbol[2][:x] == d_arr[2]["x"]
    @test Array(MatlabStructArray(d_symbol)) == d_arr
    @test s_arr == MatlabStructArray(d_arr)

    # class object array conversion
    s_arr_class = MatlabStructArray(d_arr, "TestClass")
    c_arr = Array(s_arr_class)
    @test c_arr isa Array{<:MatlabClassObject}
    @test all(c->c.class=="TestClass", c_arr)
    @test MatlabStructArray(c_arr) == s_arr_class
    @test s_arr_class != s_arr

    # test error of unequal structs
    wrong_sarr = Dict{String, Any}[
        Dict("x"=>[1.0,2.0], "y"=>[3.0,4.0]),
        Dict("x"=>[5.0,6.0])
    ]
    msg = "Cannot convert Dict array to MatlabStructArray. All elements must share identical field names"
    @test_throws ErrorException(msg) MatlabStructArray(wrong_sarr)
end

@testset "MatlabStructArray integer indexing" begin
    d_arr = Dict{String, Any}[
        Dict("x"=>[1.0,2.0], SubString("y")=>3.0),
        Dict("x"=>[5.0,6.0], "y"=>[])
    ]
    s_arr = MatlabStructArray(d_arr)
    @test s_arr[1] == Array(s_arr)[1]
    @test s_arr[2] == Array(s_arr)[2]

    s_arr = MatlabStructArray(reshape(d_arr, 1, 2))
    @test s_arr[1] == Array(s_arr)[1]
    @test s_arr[2] == Array(s_arr)[2]

    s_arr = MatlabStructArray(reshape(d_arr, 2, 1))
    @test s_arr[1] == Array(s_arr)[1]
    @test s_arr[2] == Array(s_arr)[2]

    s_arr = MatlabStructArray(reshape(d_arr, 2, 1), "TestClass")
    @test s_arr[1] == Array(s_arr)[1]
    @test s_arr[2] == Array(s_arr)[2]
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

    obj["b"] = "str"
    @test obj["b"] == "str"

    c_arr = [MatlabClassObject(d, "TestClassOld"), MatlabClassObject(d, "TestClassOld")]
    s_arr = MatlabStructArray(c_arr)
    @test s_arr.class == "TestClassOld"

    wrong_arr = [MatlabClassObject(d, "TestClassOld"), MatlabClassObject(d, "Bah")]
    @test_throws ErrorException MatlabStructArray(wrong_arr)

    d2 = deepcopy(obj.d)
    @test obj == MatlabClassObject(d2, "TestClassOld")
    @test obj != MatlabClassObject(obj.d, "Banana")
    d2["a"] = 5.0 + 1e-9
    @test obj != MatlabClassObject(d2, "TestClassOld")
    @test obj ≈ MatlabClassObject(d2, "TestClassOld")
    @test !(obj ≈ MatlabClassObject(d2, "Banana"))
end

@testset "MatlabOpaque string" begin
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
    dat = reshape(dat, 1, length(dat))
    obj = MatlabOpaque(Dict{String, Any}("any" => dat), "string")
    str = MAT.convert_opaque(obj)
    @test size(str) == (3,1)
    @test vec(str) == ["Jones", "Brown", "Smith"]

    # single element string array is a single string in matlab, ofcourse
    dat = [
        0x0000000000000001
        0x0000000000000002
        0x0000000000000001
        0x0000000000000001
        0x0000000000000005
        0x0065006e006f004a
        0x0000000000000073
    ]
    dat = reshape(dat, 1, length(dat))
    obj = MatlabOpaque(Dict{String, Any}("any" => dat), "string")
    str = MAT.convert_opaque(obj)
    @test str == "Jones"
end

@testset "MatlabOpaque datetime" begin
    @testset "datetime array" begin
        d = Dict{String, Any}(
            "tz"         => "",
            "data"       => ComplexF64[
                1482192000000.0+0.0im;
                1482278400000.0+0.0im;
                1482364800000.0+0.0im;;
                ],
            "fmt"        => "",
            #"isDateOnly" => true, # Note: "isDateOnly" not available in all versions
        )
        obj = MatlabOpaque(d, "datetime")
        expected_dates = [
            DateTime(2016, 12, 20) # 20-Dec-2016
            DateTime(2016, 12, 21) # 21-Dec-2016
            DateTime(2016, 12, 22) # 22-Dec-2016
        ]
        dt = MAT.convert_opaque(obj)
        @test all(dt .== expected_dates)

        obj_conv = MatlabOpaque(dt)
        @test obj_conv == obj
    end

    @testset "datetime element" begin
        d = Dict{String, Any}(
            "tz"         => "",
            "data"       => 1575304969634.0+0.0im,
            "fmt"        => "",
            #"isDateOnly" => false,
        )
        obj = MatlabOpaque(d, "datetime")
        # "02-Dec-2019 16:42:49"
        expected_dt = DateTime(2019, 12, 2, 16, 42, 49)
        # still have some millisecond rounding issue?
        dt = MAT.convert_opaque(obj)
        @test dt - expected_dt < Second(1)

        # convert back to MatlabOpaque
        obj_conv = MatlabOpaque(dt)
        @test obj_conv == obj
    end
end

@testset "MatlabOpaque duration" begin
    d = Dict{String,Any}(
        "millis" => [3.6e6 7.2e6],
        # "fmt"    => 'h' # optional format
    )
    obj = MatlabOpaque(d, "duration")
    ms = MAT.convert_opaque(obj)
    @test ms == map(Millisecond, d["millis"])
    @test MatlabOpaque(ms) == obj

    d = Dict{String,Any}(
        "millis" => 12000.0,
        # "fmt"    => 'h',
    )
    obj = MatlabOpaque(d, "duration")
    ms = MAT.convert_opaque(obj)
    @test ms == Millisecond(d["millis"])
    @test MatlabOpaque(ms) == obj
end

@testset "MatlabOpaque categorical" begin
    d = Dict{String,Any}(
        "isProtected"   => false,
        "codes"         => reshape(UInt8[0x02, 0x03, 0x01, 0x01, 0x01, 0x02], 3, 2),
        "categoryNames" => Any["Fair"; "Good"; "Poor";;],
        "isOrdinal"     => false,
    )
    obj = MatlabOpaque(d, "categorical")

    c = MAT.convert_opaque(obj)
    @test c == [ 
        "Good"  "Fair"
        "Poor"  "Fair"
        "Fair"  "Good"
    ]
    
end

@testset "MatlabOpaque table" begin
    # simplified table struct; there's some other properties as well
    d = Dict{String,Any}(
        "varnames" => Any["FlightNum" "Customer"],
        "nrows" => 3.0,
        "data"  => reshape(Any[[1261.0; 547.0; 3489.0;;], ["Jones"; "Brown"; "Smith";;]], 1, 2),
        "ndims" => 2.0,
        "nvars" => 2.0,
    )
    obj = MatlabOpaque(d, "table")

    # Note: this should work with DataFrames.DataFrame, but that's a big dependency to add for testing
    t = MAT.convert_opaque(obj; table = MatlabTable)
    @test t.names == [:FlightNum, :Customer]
    @test t[:FlightNum] isa Vector{Float64}
    @test t[:FlightNum] == [1261.0, 547.0, 3489.0]
    @test t[:Customer] isa Vector{String}
    @test t["Customer"] == ["Jones", "Brown", "Smith"]

    t = MAT.convert_opaque(obj; table = MatlabStructArray{1})
    @test t isa MatlabStructArray{1}
    @test t["FlightNum"] == [1261.0, 547.0, 3489.0]
    @test t["Customer"] == ["Jones", "Brown", "Smith"]

    t = MAT.convert_opaque(obj; table = Nothing)
    @test t === obj

    nd_array = reshape(1:12, 2, 3, 2)    

    # ND-arrays as columns
    # Note: does not convert to DataFrame
    d = Dict{String,Any}(
        "varnames" => Any["Floats" "NDArray"],
        "nrows" => 2.0,
        "data"  => reshape(Any[[1261.0; 547.0;;], nd_array], 1, 2),
        "ndims" => 2.0,
        "nvars" => 2.0,
    )
    obj = MatlabOpaque(d, "table")
    t = MAT.convert_opaque(obj; table = MatlabTable)
    @test size(t[:Floats]) == (2,)
    @test size(t[:NDArray]) == (2,3,2)

    # single row table
    d = Dict{String,Any}(
        "varnames" => Any["Age" "Name" "Matrix"],
        "nrows" => 1.0,
        "data"  => reshape([25.0, "Smith", [1.0 2.0]], 1, 3),
        "ndims" => 2.0,
        "nvars" => 2.0,
    )
    obj = MatlabOpaque(d, "table")
    t = MAT.convert_opaque(obj; table = MatlabTable)
    @test t[:Age] == [25.0]
    @test t[:Name] == ["Smith"]
    @test t[:Matrix] == [1.0 2.0]
end