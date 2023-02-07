cd(joinpath(@__DIR__,".."))
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# test only MAT_v4
using Test
include("../src/MAT.jl")
include("../src/MAT_v4_pr132.jl")

mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"

@testset "checkv4" begin
  # function matopen(filename::AbstractString, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool, compress::Bool)
  # function matopen(fname::AbstractString, mode::AbstractString; compress::Bool = false)
  # mat = MAT.matopen(mat1s, "r" )

  rawfid = open(mat1s, "r")
  (isv4, swap_bytes) = MAT_v4.checkv4(rawfid)
  close(rawfid)

  @test isv4 == true
end

rawfid = open(mat1s, "r")
(isv4, swap_bytes) = MAT_v4.checkv4(rawfid)
mat = MAT_v4_pr132.matopen(rawfid, swap_bytes )

@testset "getvarnames" begin
  # function getvarnames(matfile::Matlabv4File)
  @show gvn = MAT_v4_pr132.getvarnames(mat)
  # gvn = MAT_v4.getvarnames(mat) = Dict("name" => 71, "Aclass" => 0, "dataInfo" => 448328, "data_2" => 501296, "data_1" => 488197, "description" => 117126)
# ...and this is not correct
  @test false
end


close(rawfid)