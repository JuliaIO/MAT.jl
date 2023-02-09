cd(joinpath(@__DIR__,".."))
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# test only MAT_v4
using Test
using JSON
# include("../src/MAT.jl")
include("../src/MAT_v4_Modelica.jl")




# The Modelica MATv4 file takes the basic v4 Matrix format and adds some requirments to the contents and ordering of the matrices
# The first matrix, Aclass is narrowly defined 

@testset "isLittleEndian" begin
  @test MAT_v4_Modelica.isLittleEndian(0) == true
  @test MAT_v4_Modelica.isLittleEndian(1000) == false
  @test MAT_v4_Modelica.isLittleEndian(2000) == false
  @test MAT_v4_Modelica.isLittleEndian(3000) == false
end

@testset "dataFormat" begin
  @test MAT_v4_Modelica.dataFormat(0000) <: Float64
  @test MAT_v4_Modelica.dataFormat(0010) <: Float32
  @test MAT_v4_Modelica.dataFormat(0020) <: Int32
  @test MAT_v4_Modelica.dataFormat(0030) <: Int16
  @test MAT_v4_Modelica.dataFormat(0040) <: UInt16
  @test MAT_v4_Modelica.dataFormat(0050) <: UInt8
end

@testset "typeBytes" begin
  @test MAT_v4_Modelica.typeBytes(Int32) == 4
end


@testset "Aclass" begin
  mat = "W:/sync/mechgits/library/julia/MAT.jl/test/Modelica/BouncingBall/BouncingBall_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat)
  @test ac.positionStart == 0
  @test ac.positionEnd == 71
end


@testset "readVariableNames" begin
  mat = "W:/sync/mechgits/library/julia/MAT.jl/test/Modelica/BouncingBall/BouncingBall_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  # @show vn
  @test length(vn.names) == 11
  @test vn.names[1] == "time"
  @test vn.names[3] == "vel"
  @test vn.names[11] == "grav"
  @test vn.positionStart == 71
  @test vn.positionEnd == 228
end


@testset "getVariableIndex" begin
  mat = "W:/sync/mechgits/library/julia/MAT.jl/test/Modelica/BouncingBall/BouncingBall_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  @test MAT_v4_Modelica.getVariableIndex(vn, vn.names[3]) == 3
  @test MAT_v4_Modelica.getVariableIndex(vn, vn.names[10]) == 10
end


@testset "readVariableDescriptions" begin
  mat = "W:/sync/mechgits/library/julia/MAT.jl/test/Modelica/BouncingBall/BouncingBall_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  @test length(vd.descriptions) == 11  
  @test vd.descriptions[1] == "Simulation time [s]"
  @test vd.descriptions[3] == "velocity of ball"
  @test vd.descriptions[11] == "gravity acceleration"
end

@testset "readDataInfo" begin
  mat = "W:/sync/mechgits/library/julia/MAT.jl/test/Modelica/BouncingBall/BouncingBall_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  di = MAT_v4_Modelica.readDataInfo(ac,vd)
  # @show di.info[3]
  @test di.info[1]["isWithinTimeRange"] == -1
  @test di.info[3]["locatedInData"] == 2 
  @test di.info[4]["isInterpolated"] == 0
  @test di.info[11]["isWithinTimeRange"] == 0
end

@testset "readVariable: BouncingBall" begin
  mat = "W:/sync/mechgits/library/julia/MAT.jl/test/Modelica/BouncingBall/BouncingBall_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  di = MAT_v4_Modelica.readDataInfo(ac,vd)

  # println(JSON.json(di.info, 2)) #get the data 1/2 info

  eff = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "eff") #data1
  @test length(eff) == 2
  @test eff[1] ≈ 0.77
  @test eff[2] ≈ 0.77

  grav = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "grav") #data1
  @test length(grav) == 2
  @test grav[1] ≈ 9.81
  @test grav[2] ≈ 9.81

  time = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "time") # data0
  @test all(isapprox.(time, [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1], rtol=1e-3))

  height = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "height") #data2
  @test isapprox(height[1], 111, rtol=1e-3)
  @test isapprox(height[2], 110.9509, rtol=1e-3)

  vel = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "vel") #data2
  @test isapprox(vel[2], -0.981, rtol=1e-3)
end

@testset "readVariable: FallingBodyBox" begin
  mat = "W:/sync/mechgits/library/julia/MAT.jl/test/Modelica/FallingBodyBox/FallingBodyBox_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  di = MAT_v4_Modelica.readDataInfo(ac,vd)

  # println(JSON.json(di.info, 2)) 

  var = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "time") 
  # display(var)
  ret = true
  for i = 2:length(var)-1 #last time is duplicated
    ret &= isapprox(var[i]-var[i-1], 0.002, rtol=1e-4)
  end
  @test ret == true

  #point-check values read from FallingBodyBox_res.csv
  var = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.frame_a.r_0[1]") 
  @test isapprox(var[16], 0.002923239, rtol=1e-3)

  var = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.frame_a.R.T[1,1]") 
  @test isapprox(var[26], 0.983794001, rtol=1e-3)

  @test_throws ErrorException MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.frame_a.v_0[1]")  # there is no frame_A.v_0

  var = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.v_0[2]") 
  @test isapprox(var[33], -0.58818129, rtol=1e-3)

  var = MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.frame_b.r_0[1]") 
  @test isapprox(var[72], 0.935886479, rtol=1e-3)
end

;

