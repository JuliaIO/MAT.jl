using MAT, Base.Test

function check(filename, result)
	matfile = matopen(filename)
	for (k, v) in result
		@test exists(matfile, k)
		got = read(matfile, k)
		if !isequal(got, v) || typeof(got) != typeof(v)
			close(matfile)
			error("""
				Data mismatch reading $k from $filename ($format)

				Got $(typeof(got)):

				$(repr(got))

				Expected $(typeof(v)):

				$(repr(v))
				""")
		end
	end
	@test union!(Set(), names(matfile)) == union!(Set(), keys(result))
	close(matfile)

	mat = matread(filename)
	if !isequal(mat, result)
		error("""
			Data mismatch reading $filename ($format)

			Got:

			$(repr(mat))

			Expected:

			$(repr(result))
			""")
		close(matfile)
		return false
	end

	return true
end

global format
for format in ["v6", "v7", "v7.3"]
	cd(joinpath(dirname(@__FILE__), format))

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
		"double" => float64(1),
		"logical" => true
	}
	check("simple.mat", result)
	matfile = matopen("simple.mat")
	mat = read(matfile)
	close(matfile)
	for (k, v) in result
		if(typeof(mat[k]) != typeof(v))
			error("Types for $k didn't match (expected $(typeof(v)), got $(typeof(mat[k])))")
		end
	end

	result = {
		"imaginary" => Complex128[1 -1 1+im 1-im -1+im -1-im im]
	}
	check("complex.mat", result)

	result = {
		"simple_string" => "the quick brown fox",
		"accented_string" => "thé qüîck browñ fòx",
		"concatenated_strings" => ByteString["this is a string", "this is another string"],
		"cell_strings" => {"this is a string" "this is another string"},
		"empty_string" => ""
	}
	check("string.mat", result)

	result = {
		"a1x2" => [1.0 2.0],
		"a2x1" => zeros(2, 1)+[1.0, 2.0],
		"a2x2" => [1.0 3.0; 4.0 2.0],
		"a2x2x2" => cat(3, [1.0 3.0; 4.0 2.0], [1.0 2.0; 3.0 4.0]),
		"empty" => zeros(0, 0),
		"string" => "string"
	}
	check("array.mat", result)

	result = {
		"cell" => {1.0 2.01 "string" {"string1" "string2"}}
	}
	check("cell.mat", result)

	result = {
		"s" => (ASCIIString=>Any)[
			"a" => 1.0,
			"b" => [1.0 2.0],
			"c" => [1.0 2.0 3.0]
		],
		"s2" => (ASCIIString=>Any)[ "a" => {1.0 2.0} ]
	}
	check("struct.mat", result)

	result = {
		"logical" => false,
		"logical_mat" => [
			true false false
			false true false
			true false false
		]
	}
	check("logical.mat", result)
	
	result = {
		"empty_cells" => {zeros(0, 0) "test" zeros(0, 0)}
	}
	check("empty_cells.mat", result)

	result = {
		"sparse_empty" => sparse(Array(Float64, 0, 0)),
		"sparse_eye" => speye(20),
		"sparse_logical" => SparseMatrixCSC{Bool,Int64}(5, 5, [1:6], [1:5], bitunpack(trues(5))),
		"sparse_random" => sparse([0 6. 0; 8. 0 1.; 0 0 9.]),
		"sparse_complex" => sparse([0 6. 0; 8. 0 1.; 0 0 9.]*(1. + 1.im)),
		"sparse_zeros" => SparseMatrixCSC(20, 20, ones(Int, 21), Int[], Float64[])
	}
	check("sparse.mat", result)

	matfile = matopen("partial.mat")
	var1 = read(matfile, "var1")
	@assert var1[28, 33] == 5
	var2 = read(matfile, "var2")
	@assert var2[27, 90] == 10
	close(matfile)
    
    if format != "v7.3"
        # Class tests, not implemented for v7.3 yet
        result = {
            "s" => (ASCIIString=>Any)[
                "char_field_1" => "char_field",
                "array_field_2" => [0.0 1.0 2.0 3.0 4.0 5.0],
                "cell_field_3" => {1.0  "one"  2.0  "two"},
                "_class" => ("", "SimpleClass")
            ]
        }
        check("simpleclass.mat", result)
    end
    
end
