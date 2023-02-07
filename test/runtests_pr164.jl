cd(joinpath(@__DIR__,".."))
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# test only MAT_v4
using Test
include("../src/MAT.jl")
include("../src/MAT_v4_pr164.jl")

mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"



@testset "unpack_header" begin
  rawfid = open(mat1s, "r")
  # (isv4, swap_bytes) = MAT_v4.checkv4(rawfid)

  # From 'Version-4-MAT-File-Format.pdf':
  # Each matrix starts with a fixed-length 20-byte header that contains information describing certain attributes of the Matrix.
  # The 20-byte header consists of five long (4-byte) integers:
  @show type, mrows, ncols, imagf, namlen = read!(rawfid, Vector{Int32}(undef, 5)) # int32=4byte * 5 = 20byte
  #The type flag contains an integer whose decimal digits encode storage information. If the integer is represented as MOPT where M is the thousands digit...
  T, P, O, M = digits(type; pad=4)
  # P indicates which format the data is stored: 5 = 8bit uint
  # T indicates matrix type: 1 = text matrix

  @show eltype = MAT_v4_pr164.unpack_header_type(type)


  close(rawfid)

  # @test isv4 == true
end

# rawfid = open(mat1s, "r")
# (isv4, swap_bytes) = MAT_v4.checkv4(rawfid)
# mat = MAT_v4.matopen(rawfid, swap_bytes )

# @testset "getvarnames" begin
#   # function getvarnames(matfile::Matlabv4File)
#   @show gvn = MAT_v4.getvarnames(mat)
#   # gvn = MAT_v4.getvarnames(mat) = Dict("name" => 71, "Aclass" => 0, "dataInfo" => 448328, "data_2" => 501296, "data_1" => 488197, "description" => 117126)
# # ...and this is not correct
#   @test false
# end


# close(rawfid)