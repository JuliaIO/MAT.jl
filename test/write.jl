using MAT

tmpfile = string(tempname, ".mat")

function test_write(data)
	matwrite(tmpfile, data)

	fid = matopen(tmpfile, "r")
	local result
	try
		result = read(fid)
	finally
		close(fid)
	end

	if !isequal(result, data)
		error("Data mismatch")
	end
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
	"Complex128" => [1.0 -1.0 1.0+1.0im 1.0-1.0im -1.0+1.0im -1.0-1.0im 1.0im],
	"ComplexPair" => [1 2-3im 4+5im]
))
test_write(Dict("Complex128" => 1.0im, "ComplexPair" => 2-3im))

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
	"a2x2x2" => cat(3, [1.0 3.0; 4.0 2.0], [1.0 2.0; 3.0 4.0]),
	"empty" => zeros(0, 0),
	"string" => "string"
))

test_write(Dict(
	"cell" => Any[1 2.01 "string" Any["string1" "string2"]]
))

test_write(Dict(
	"s" => Dict(
		"a" => 1.0,
		"b" => [1.0 2.0],
		"c" => [1.0 2.0 3.0]
	),
	"s2" => Dict("a" => [1.0 2.0])
))

test_write(Dict(
	"sparse_empty" => sparse(Matrix{Float64}(0, 0)),
	"sparse_eye" => speye(20),
	"sparse_logical" => SparseMatrixCSC{Bool,Int64}(5, 5, [1:6;], [1:5;], fill(true, 5)),
	"sparse_random" => sparse([0 6. 0; 8. 0 1.; 0 0 9.]),
	"sparse_complex" => sparse([0 6. 0; 8. 0 1.; 0 0 9.]*(1. + 1.0im)),
	"sparse_zeros" => SparseMatrixCSC(20, 20, ones(Int, 21), Int[], Float64[])
))

@test_throws ErrorException test_write(Dict("1invalidkey" => "starts with a number"))
@test_throws ErrorException test_write(Dict("another invalid key" => "invalid characters"))
@test_throws ErrorException test_write(Dict("yetanotherinvalidkeyyetanotherinvalidkeyyetanotherinvalidkeyyetanotherinvalidkey" => "too long"))

type TestCompositeKind
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
	"Complex128" => [1.0 -1.0 1.0+1.0im 1.0-1.0im -1.0+1.0im -1.0-1.0im 1.0im],
	"simple_string" => "the quick brown fox",
	"a1x2" => [1.0 2.0],
	"sparse_empty" => sparse(Matrix{Float64}(0, 0))
))
test_write(sd)
