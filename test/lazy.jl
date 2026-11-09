using MAT, Test
import HDF5

@testset "readlazy" begin
    mktempdir() do dir
        path = joinpath(dir, "lazy.mat")
        A = reshape(collect(1.0:24.0), 2, 3, 4)
        big = reshape(collect(1.0:60_000.0), 200, 300)
        matwrite(path, Dict(
            "A" => A, "big" => big, "lg" => [true false; false true],
            "cx" => [1.0+2.0im 3.0-4.0im], "i8" => Int8[1 2 3],
            "c" => Any[1.0, "a"], "s" => Dict("a" => 1.0), "t" => "text", "e" => zeros(0, 3),
        ); compress = true)  # compressed datasets are chunked

        matopen(path) do fid
            if HDF5.Dataset <: AbstractArray  # HDF5.jl 0.18 and later
                lazyA = MAT.readlazy(fid, "A")
                @test size(lazyA) == size(A)
                @test lazyA[2, 1:2, 3] == A[2, 1:2, 3]
                @test collect(lazyA) == A
                @test sum(view(MAT.readlazy(fid, "big"), :, 7)) == sum(big[:, 7])

                lg = MAT.readlazy(fid, "lg")
                @test eltype(lg) == Bool
                @test collect(lg) == [true false; false true]
                cx = MAT.readlazy(fid, "cx")
                @test eltype(cx) == ComplexF64
                @test collect(cx) == [1.0+2.0im 3.0-4.0im]
                @test eltype(MAT.readlazy(fid, "i8")) == Int8

                for name in ("c", "s", "t", "e")
                    @test_throws ArgumentError MAT.readlazy(fid, name)
                end
            else
                @test_throws r"needs HDF5.jl 0.18" MAT.readlazy(fid, "A")
            end
        end
    end
    @test_throws ArgumentError matopen(fid -> MAT.readlazy(fid, "x"), joinpath(@__DIR__, "v7", "simple.mat"))
end
