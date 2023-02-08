cd(joinpath(@__DIR__,".."))
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# test only MAT_v4
using Test
# include("../src/MAT.jl")
# include("../src/MAT_v4_Modelica.jl")

mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"

struct MATrix
  name::String
  # data::T where T<:AbstractMatrix
  data::Any
  # type
end

@enum numberFormat little=0 big=1 vaxD=2 vaxG=3 cray=4
@enum dataFormat double=0 single=1 int32=2 int16=3 uint16=4 uint8=5
@enum matrixFormat numericFull=0 text=1 spars=2

"""
Read a single data matrix from the opened MAT file
"""
function readMATMatrix(ios::IOStream)
  @show startP = position(ios)
  # The 20-byte header consists of five long (4-byte) integers:
  type, nrows, ncols, imagf, namelen = read!(ios, Vector{Int32}(undef, 5)) # 32bits=4byte, * 5 qty

  #The type flag contains an integer whose decimal digits encode storage information. If the integer is represented as MOPT where M is the thousands digit...
  T, P, O, M = digits(type; pad=4)

  nameuint = read!(ios, Vector{UInt8}(undef, namelen)) # read the full namelen to make the pointer ready to read the data
  matrixName = replace(String(nameuint), '\0'=>"")

  #  real Real part of the matrix consists of nrows ∗ ncols numbers in the format specified by the P element of the type flag. The data is stored column-wise such that the second column follows the first column, etc.
  realint = []
  @show dataFormat(P)
  if dataFormat(P) == uint8 
    realint = read!(ios, Matrix{UInt8}(undef, nrows,ncols))  
  elseif dataFormat(P) == uint16 
    realint = read!(ios, Matrix{UInt16}(undef, nrows,ncols))  
  elseif dataFormat(P) == int32 
    realint = read!(ios, Matrix{Int32}(undef, nrows,ncols))  
  elseif dataFormat(P) == single 
    realint = read!(ios, Matrix{Float32}(undef, nrows,ncols))  
  elseif dataFormat(P) == double 
    realint = read!(ios, Matrix{Float64}(undef, nrows,ncols))  
  else
    error("Unknown data format for P=$P")
  end

  #  imag Imaginary part of the matrix, if any. If the imaginary flag imagf is nonzero, the imaginary part of a matrix is placed here. It is stored in the same manner as the real data. 
  if imagf==1
    error("imaginary not implemented")
  #   if dataFormat(P)
  #   @show imagint = read!(ios, Vector{UInt8}(undef, nrows*ncols)) # read the full namelen to make the pointer ready to read the data
  end

  @show matrixName
  @show endP = position(ios)

  return MATrix(matrixName, realint)
end

"""
  # The variables stored in the MAT-file are (in the order required by OpenModelica):
  #   Aclass
  # • Aclass(1,:) is always Atrajectory
  # • Aclass(2,:) is 1.1 in OpenModelica
  # • Aclass(3,:) is empty
  # • Aclass(4,:) is either binTrans or binNormal


"""
function parseAclass(mat::MATrix)
  # @show mat
  # Aclass1 = replace(String(mat.data[1,:]), '\0'=>"")
  # Aclass2 = replace(String(data[2,:]), '\0'=>"")
  # Aclass3 = replace(String(data[3,:]), '\0'=>"")
  # Aclass4 = replace(String(data[4,:]), '\0'=>"")
  Aclass = [  replace(String(mat.data[1,:]), '\0'=>""),
              replace(String(mat.data[2,:]), '\0'=>""),
              replace(String(mat.data[3,:]), '\0'=>""),
              replace(String(mat.data[4,:]), '\0'=>"") ]
  return Aclass
end


"""
name
Is an n x m character (int8) matrix, where n is the number of variables stored in the result-file (including time). 
m is the length of the longest variable.
OpenModelica stores the trailing part of the name as NIL bytes (0) whereas other tools use spaces for the trailing part.
"""
function parseVariableNames(mat::MATrix)
  varNames = []
  m,n = size(mat.data)
  for i in 1:n
    push!(varNames, replace(String(mat.data[:,i]), '\0'=>""))
  end
  return varNames
end

"""
Find the index of the given variable, returning -1 if not found.
"""
function getVariableIndex(variableNames, name)
  vecAll = findall( x->x==name, variableNames)
  n = length(vecAll)
  
  if isempty(vecAll) == true
    return -1
  else
    if n>1
      error("Found $n instances of variable [$name], but variables should be unique.")
    end
    return vecAll[1] 
  end
end  


"""
description
Is an n x m character (int8) matrix containing the comment-string corresponding to the variable in the name matrix
"""
function parseVariableDescriptions(mat::MATrix)
  varDesc = []
  m,n = size(mat.data)
  for i in 1:n
    push!(varDesc, replace(String(mat.data[:,i]), '\0'=>""))
  end
  return varDesc 
end

"""
dataInfo provides indicies to access variable data

dataInfo
Is an n x 4 integer matrix containing information for each variable (in the same order as the name and description matrices).
• dataInfo(i,1) is 1 or 2, saying if variable i is stored in the data_1 or data_2 matrix. If it is 0, it is the abscissa (time variable).
• dataInfo(i,2) contains the index in the data_1 or data_2 matrix. The index is 1-based and may contain several variables pointing to the same row (alias variables). A negative value means that the variable is a negated alias variable.
• dataInfo(i,3) is 0 to signify linear interpolation. In other tools the value is the number of times differentiable this variable is, which may improve plotting.
• dataInfo(i,4) is -1 in OpenModelica to signify that the value is not defined outside the time range. 0 keeps the first/last value when going outside the time range and 1 performs linear interpolation on the first/last two points.
"""
function parseDataInfo(mat::MATrix, varIndex::Int)
  # @show mat.data[:,varIndex]
  return Dict("locatedInData"=>mat.data[1,varIndex], "indexInData"=>mat.data[2,varIndex], "isInterpolated"=>mat.data[3,varIndex], "isWithinTimeRange"=>mat.data[4,varIndex] )
end

"""
data_1
If it is an n x 1 matrix it contains the values of parameters. If it is an n x 2 matrix, the first and second column
signify start and stop-values.
"""
function parseData1(mat::MATrix, varIndex::Int)
  @show size(mat.data)
  @show mat.data[varIndex,:]
end


@testset "matrix reader" begin
  matio = open(mat1s, "r")
  seekstart(matio)
  ret = readMATMatrix(matio)
  Aclass = parseAclass(ret)
  @test Aclass[1] == "Atrajectory"
  @test Aclass[2] == "1.1"
  @test isempty(Aclass[3]) == true
  @test Aclass[4] == "binTrans" || Aclass[4] == "binNormal"

  ret = readMATMatrix(matio) # name
  variableNames = parseVariableNames(ret)
  @test variableNames[1] == "time"
  @test variableNames[2490] == "world.z_label.color[3]"
  @test getVariableIndex(variableNames, "time") == 1


  @show position(matio)
  ret = readMATMatrix(matio) # description
  variableDescriptions = parseVariableDescriptions(ret)
  @test variableDescriptions[2490] == "Color of cylinders"

  # @allocated readMATMatrix(matio) # dataInfo
  ret = readMATMatrix(matio) # dataInfo
  @show dataInfo = parseDataInfo(ret, 1)

  ret = readMATMatrix(matio) # data_1
  

  # itime = getVariableIndex(variableNames, "time") 
  # # ditime = parse
  # parseData1(ret, itime) 

  # parseData1(ret, iphi) 
  # # parseData1(ret, ihfaf1) 

  # iphi = getVariableIndex(variableNames, "revolute.phi") 
  # ihfaf1 = getVariableIndex(variableNames, "head.frame_a.f[1]")

  ret = readMATMatrix(matio) # data_2
  close(matio)

  @test true

end
