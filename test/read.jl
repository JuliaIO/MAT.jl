load("src/MAT")
using MAT

function check(filename, result)
	matfile = matopen(filename)
	mat = read(matfile)
	if mat != result
			error("Data mismatch reading $filename")
		close(matfile)
		return false
	end
	for (k, v) in result
		if read(matfile, k) != v
			close(matfile)
			error("Data mismatch reading $k from $filename")
		end
	end
	close(matfile)
	return true
end

for format in ["v6", "v7", "v7.3"]
	tests = "test/"*format

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
	check("$tests/simple.mat", result)
	matfile = matopen("$tests/simple.mat")
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
	#check("$tests/complex.mat", result)

	result = {
		"simple_string" => "the quick brown fox",
		"accented_string" => "thé qüîck browñ fòx",
		"concatenated_strings" => ["this is a string", "this is another string"],
		"cell_strings" => ["this is a string" "this is another string"]
	}
	check("$tests/string.mat", result)

	result = {
		"a1x2" => [1.0 2.0],
		"a2x1" => zeros(2, 1)+[1.0, 2.0],
		"a2x2" => [1.0 3.0; 4.0 2.0],
		"a2x2x2" => cat(3, [1.0 3.0; 4.0 2.0], [1.0 2.0; 3.0 4.0]),
		"empty" => zeros(0, 0),
		"string" => "string"
	}
	check("$tests/array.mat", result)

	result = {
		"cell" => {1 2.01 "string" {"string1" "string2"}}
	}
	check("$tests/cell.mat", result)

	result = {
		"s" => {
			"a" => 1.0,
			"b" => [1.0 2.0],
			"c" => [1.0 2.0 3.0]
		},
		"s2" => { "a" => [1.0 2.0] }
	}
	check("$tests/struct.mat", result)

	matfile = matopen("$tests/partial.mat")
	var1 = read(matfile, "var1")
	@assert var1[28, 33] == 5
	var2 = read(matfile, "var2")
	@assert var2[27, 90] == 10
	close(matfile)
end