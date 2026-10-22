using MAT, Test

function check(filename, result)
    matfile = matopen(filename)
    for (k, v) in result
        @test haskey(matfile, k)
        got = read(matfile, k)
        if !isequal(got, v) || (typeof(got) != typeof(v) && (!isa(got, String) || !(isa(v, String))))
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
    @test union!(Set(), keys(matfile)) == union!(Set(), keys(result))
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
cd(dirname(@__FILE__))
for filename in readdir("v4")
    endswith(filename, ".mat") || continue
    d = matread("v4/$filename")
    matwrite("v4/tmp.mat", d; version="v4")
    check("v4/tmp.mat", d)
    rm("v4/tmp.mat")
end

# support read var name without '\0', xref https://github.com/JuliaIO/MAT.jl/pull/202
let tmp = "v4/tmp.mat"
    data = rand(1,9);
    open(tmp, "w") do fid
        M, N = size(data);
        data_type = 0;
        name = "testnamelen";
        L = ncodeunits(name)
        write(fid, Int32(data_type));
        write(fid, Int32(M));
        write(fid, Int32(N));
        write(fid, Int32(0));
        write(fid, Int32(L));
        # name is store without null-terminator '\0'
        write(fid, name)
        write(fid, data)
    end
    d = MAT.matread(tmp)
    @test haskey(d, "testnamelen")
    @test first(values(d)) == data
    rm(tmp)
end

# Read zero-dimensional vars
let tmp = "v4/tmp.mat"
    open(tmp, "w") do fid
        name = "testnamelen";
        L = ncodeunits(name)
        write(fid, Int32(0));
        write(fid, Int32(0));
        write(fid, Int32(0));
        write(fid, Int32(0));
        write(fid, Int32(L));
        # name is store without null-terminator '\0'
        write(fid, name)
    end
    d = MAT.matread(tmp)
    @test haskey(d, "testnamelen")
    rm(tmp)
end

@testset "v4 sparse matrices and text as MATLAB writes them" begin
    # a v4 variable: type (1 text, 2 sparse), rows, columns, complex flag, name, data
    v4var(name, type, mrows, ncols, data; imagf=0) =
        vcat(reinterpret(UInt8, Int32[type, mrows, ncols, imagf, ncodeunits(name) + 1]),
             Vector{UInt8}(name), UInt8[0], reinterpret(UInt8, Float64.(data)))
    mktempdir() do dir
        tmp = joinpath(dir, "tmp.mat")
        # sparse matrices end with a row of (m, n, 0) giving the size
        write(tmp, vcat(
            v4var("s", 2, 2, 3, [1, 3, 1, 3, 5, 0]),                # sparse(1, 1, 5, 3, 3)
            v4var("c", 2, 2, 4, [1, 2, 1, 4, 1, 0, 2, 0]),          # sparse(1, 1, 1+2i, 2, 4)
            v4var("z", 2, 1, 3, [3, 3, 0]),                         # sparse(3, 3)
            v4var("t", 1, 1, 6, [181, 109, 32, 177, 176, 67])))     # 'µm ±°C'
        d = matread(tmp)
        @test d["s"] == sparse([1], [1], [5.0], 3, 3) && nnz(d["s"]) == 1
        @test d["c"] == sparse([1], [1], [1.0 + 2.0im], 2, 4) && nnz(d["c"]) == 1
        @test size(d["z"]) == (3, 3) && nnz(d["z"]) == 0
        @test d["t"] == "µm ±°C"

        # MAT.jl writes no size row: the last row is a stored value
        write(tmp, v4var("s", 2, 1, 3, [1, 2, 5]))
        @test matread(tmp)["s"] == sparse([1], [2], [5.0], 1, 2)
    end
end
