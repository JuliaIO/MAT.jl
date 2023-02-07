
cd(joinpath(@__DIR__,".."))
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))


using SparseArrays, LinearAlgebra

include("read.jl")
include("write.jl")



