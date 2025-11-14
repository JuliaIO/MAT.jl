using MAT, Test

# MatlabStructArray construction from dict array
d_arr = Dict{String, Any}[
    Dict("x"=>[1.0,2.0], SubString("y")=>[3.0,4.0]),
    Dict("x"=>[5.0,6.0], "y"=>[])
]
s_arr = MAT.MatlabStructArray(d_arr)
@test s_arr["y"][2] == d_arr[2]["y"]
@test s_arr["x"][1] == d_arr[1]["x"]

# constructor errors to protect the user
@test_throws ErrorException MAT.MatlabStructArray(["a", "b"], [[]])
@test_throws ErrorException MAT.MatlabStructArray(["a", "b"], [[],[0.1, 0.2]])

# equality checks
@test isequal(MAT.MatlabStructArray(["a"], [[0.1, 0.2]]), MAT.MatlabStructArray(["a"], [[0.1, 0.2]]))
@test !isequal(MAT.MatlabStructArray(["a"], [[0.1, 0.2]]), MAT.MatlabStructArray(["b"], [[0.1, 0.2]]))
@test isapprox(MAT.MatlabStructArray(["a"], [[0.1, 0.2]]), MAT.MatlabStructArray(["a"], [[0.1+eps(0.1), 0.2]]))
@test !isapprox(MAT.MatlabStructArray(["a"], [[0.1, 0.2]]), MAT.MatlabStructArray(["b"], [[0.1, 0.2]]))
@test !isapprox(MAT.MatlabStructArray(["a"], [[0.1, 0.2]]), MAT.MatlabStructArray(["a"], [[0.11, 0.2]]))

# empty struct array constructor
s_arr = MAT.MatlabStructArray(["x", "y"], (0,1))
@test s_arr["x"] == Matrix{Any}(undef, 0, 1)
@test s_arr["y"] == Matrix{Any}(undef, 0, 1)
@test MAT.MatlabStructArray(["a"]) == MAT.MatlabStructArray(["a"], (0,0))

# convert to Dict to support easy conversion to legacy read behavior
s_arr = MAT.MatlabStructArray(d_arr)
d = Dict(s_arr)
@test d isa Dict{String, Any}
@test collect(keys(d)) == s_arr.names
@test collect(values(d)) == s_arr.values
# iteration similar to Dict
@test length(s_arr) == 2
@test collect(Dict(s_arr)) == collect(s_arr)

# possibility to convert back to dict array via `Array`
s_arr = MAT.MatlabStructArray(d_arr)
@test Array(s_arr) == d_arr
d_symbol = Array{Dict{Symbol,Any}}(MAT.MatlabStructArray(d_arr))
@test d_symbol[2][:x] == d_arr[2]["x"]
@test Array(MAT.MatlabStructArray(d_symbol)) == d_arr

# test error of unequal structs
wrong_sarr = Dict{String, Any}[
    Dict("x"=>[1.0,2.0], "y"=>[3.0,4.0]),
    Dict("x"=>[5.0,6.0])
]
msg = "Cannot convert Dict array to MatlabStructArray. All elements must share identical field names"
@test_throws ErrorException(msg) MAT.MatlabStructArray(wrong_sarr)