
cd(joinpath(@__DIR__,".."))
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))


using SparseArrays, LinearAlgebra

include("read.jl")
include("readwrite4.jl")
include("write.jl")
include("runtests_modelica.jl")


