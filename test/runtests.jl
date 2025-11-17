using SparseArrays, LinearAlgebra
using Test, MAT

@testset "MAT" begin
    include("types.jl")
    include("read.jl")
    include("readwrite4.jl")
    include("write.jl")
end
