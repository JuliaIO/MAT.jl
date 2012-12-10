load("src/MAT")
using MAT

function test_write(data)
	fid = matopen("/tmp/matwrite.mat", "w")
	for (k, v) in data
		write(fid, k, v)
	end
	close(fid)
	fid = matopen("/tmp/matwrite.mat", "r")
	result = read(fid)
	close(fid)
	if result != data
		error("Data mismatch")
	end
end

test_write({
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
})

test_write({
	"simple_string" => "the quick brown fox",
	"accented_string" => "thé qüîck browñ fòx",
	"concatenated_strings" => ["this is a string", "this is another string"],
	"cell_strings" => ["this is a string" "this is another string"]
})

test_write({
	"a1x2" => [1.0 2.0],
	"a2x1" => zeros(2, 1)+[1.0, 2.0],
	"a2x2" => [1.0 3.0; 4.0 2.0],
	"a2x2x2" => cat(3, [1.0 3.0; 4.0 2.0], [1.0 2.0; 3.0 4.0]),
	"empty" => zeros(0, 0),
	"string" => "string"
})

test_write({
	"cell" => {1 2.01 "string" {"string1" "string2"}}
})

test_write({
	"s" => {
		"a" => 1.0,
		"b" => [1.0 2.0],
		"c" => [1.0 2.0 3.0]
	},
	"s2" => { "a" => [1.0 2.0] }
})