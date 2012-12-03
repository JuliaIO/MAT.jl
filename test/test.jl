load("src/MAT")
using MAT

function check(filename, result)
	matfile = matopen(filename)
	if read(matfile) != result
		close(matfile)
		return false
	end
	for (k, v) in result
		if read(matfile, k) != v
			close(matfile)
			return false
		end
	end
	close(matfile)
	return true
end

result = {
	"int8" => int8(1),
	"uint8" => uint8(1),
	"int16" => int16(1),
	"uint16" => uint16(1),
	"int32" => int32(1),
	"uint32" => uint32(1),
	"int64" => int64(1),
	"uint64" => uint64(1),
	"single" => float32(1),
	"double" => float64(1)
}
@assert check("test/simple.mat", result)
matfile = matopen("test/simple.mat")
mat = read(matfile)
close(matfile)
for (k, v) in result
	if(typeof(mat[k]) != typeof(v))
		error("Types for $k didn't match (expected $(typeof(v)), got $(typeof(mat[k])))")
	end
end

result = {
	"imaginary" => [1 -1 1+im 1-im -1+im -1-im im]
}
@assert check("test/complex.mat", result)

result = {
	"simple_string" => "the quick brown fox",
	"accented_string" => "thé qüîck browñ fòx",
	"concatenated_strings" => ["this is a string", "this is another string"],
	"cell_strings" => ["this is a string" "this is another string"]
}
@assert check("test/string.mat", result)
@assert check("test/string_unicode.mat", result)

result = {
	"a1x2" => [1.0 2.0],
	"a2x1" => zeros(2, 1)+[1.0, 2.0],
	"a2x2" => [1.0 3.0; 4.0 2.0],
	"a2x2x2" => cat(3, [1.0 3.0; 4.0 2.0], [1.0 2.0; 3.0 4.0]),
	"empty" => zeros(0, 0),
	"string" => "string"
}
@assert check("test/array.mat", result)

result = {
	"cell" => {1 2.01 "string" {"string1" "string2"}}
}
@assert check("test/cell.mat", result)

result = {
	"s" => {
		"a" => 1.0,
		"b" => [1.0 2.0],
		"c" => [1.0 2.0 3.0]
	},
	"s2" => [{ "a" => 1 } { "a" => 2 }]
}
@assert check("test/struct.mat", result)

result = {
	"vector" => float64([1:10]'),
	"struct" => { "my"=>"struct" },
	"cell" => {"my" "cell"}
}
@assert check("test/compressed.mat", result)

matfile = matopen("test/partial.mat")
var1 = read(matfile, "var1")
@assert var1[28, 33] == 5
var2 = read(matfile, "var2")
@assert var2[27, 90] == 10
close(matfile)

matfile = matopen("test/partial_compressed.mat")
var1 = read(matfile, "var1")
@assert var1[28, 33] == 5
var2 = read(matfile, "var2")
@assert var2[27, 90] == 10
close(matfile)