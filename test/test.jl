load("src/MAT")
using MAT

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
@time mat = read(matopen("test/simple.mat"))
@assert result == mat
for (k, v) in result
	if(typeof(mat[k]) != typeof(v))
		error("Types for $k didn't match (expected $(typeof(v)), got $(typeof(mat[k])))")
	end
end

result = {
	"imaginary" => [1 -1 1+im 1-im -1+im -1-im im]
}
@time mat = read(matopen("test/complex.mat"))
@assert result == mat

result = {
	"simple_string" => "the quick brown fox",
	"accented_string" => "thé qüîck browñ fòx",
	"concatenated_strings" => ["this is a string", "this is another string"],
	"cell_strings" => ["this is a string" "this is another string"]
}
@time mat = read(matopen("test/string.mat"))
@assert result == mat
@time mat = read(matopen("test/string_unicode.mat"))
@assert result == mat

result = {
	"a1x2" => [1.0 2.0],
	"a2x1" => zeros(2, 1)+[1.0, 2.0],
	"a2x2" => [1.0 3.0; 4.0 2.0],
	"a2x2x2" => cat(3, [1.0 3.0; 4.0 2.0], [1.0 2.0; 3.0 4.0]),
	"empty" => zeros(0, 0),
	"string" => "string"
}
@time mat = read(matopen("test/array.mat"))
@assert result == mat

result = {
	"cell" => {1 2.01 "string" {"string1" "string2"}}
}
@time mat = read(matopen("test/cell.mat"))
@assert result == mat

result = {
	"s" => {
		"a" => 1.0,
		"b" => [1.0 2.0],
		"c" => [1.0 2.0 3.0]
	},
	"s2" => [{ "a" => 1 } { "a" => 2 }]
}
@time mat = read(matopen("test/struct.mat"))
@assert result == mat

result = {
	"vector" => float64([1:10]'),
	"struct" => { "my"=>"struct" },
	"cell" => {"my" "cell"}
}
@time mat = read(matopen("test/compressed.mat"))
@assert result == mat

matfile = matopen("test/partial.mat")
@time var1 = read(matfile, "var1")
@assert var1[28, 33] == 5
@time var2 = read(matfile, "var2")
@assert var2[27, 90] == 10
close(matfile)

matfile = matopen("test/partial_compressed.mat")
@time var1 = read(matfile, "var1")
@assert var1[28, 33] == 5
@time var2 = read(matfile, "var2")
@assert var2[27, 90] == 10
close(matfile)