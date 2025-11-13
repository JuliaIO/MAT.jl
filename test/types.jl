# MatlabStructArray construction from dict array
d_arr = Dict{String, Any}[
    Dict("x"=>[1.0,2.0], SubString("y")=>[3.0,4.0]),
    Dict("x"=>[5.0,6.0], "y"=>[])
]
s_arr = MAT.MatlabStructArray(d_arr)
@test s_arr["y"][2] == d_arr[2]["y"]
@test s_arr["x"][1] == d_arr[1]["y"]

# convert to Dict (legacy read behavior)
d = Dict(s_arr)
@test d isa Dict{String, Any}
@test collect(keys(d)) == s_arr.names
@test collect(values(d)) == s_arr.values

# convert back to dict array
@test Array(s_arr) == d_arr