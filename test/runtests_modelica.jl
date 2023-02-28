
# Test origins:
#  test/v4_Modelica/BouncingBall/BouncingBall_res.mat - simulated with OpenModelica v1.19.0
#  test/v4_Modelica/BouncingBall/BouncingBall_dymola2021.mat - simulated with Dymola v2021
#  test/v4_Modelica/FallingbodyBox/FallingBodyBox_res.mat - simulated with OpenModelica v1.19.0
# These exercise every function in MAT_v4_Modelica.jl...but often use hand-observed values or otherwise require knowledge of the mat's contents

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using JSON
cd(joinpath(@__DIR__,".."))
include("../src/MAT.jl")
# include("../src/MAT_v4_Modelica.jl")

#OpenModelica v1.19.0
matBBO = joinpath(@__DIR__, "v4_Modelica","BouncingBall","BouncingBall_res.mat")
matFBB = joinpath(@__DIR__, "v4_Modelica","FallingBodyBox","FallingBodyBox_res.mat")
#Dymola v2021 -- not implemented
# matBBD = joinpath(@__DIR__, "v4_Modelica","BouncingBall","BouncingBall_dymola2021.mat")


@testset "isLittleEndian" begin
  @test MAT.MAT_v4_Modelica.isLittleEndian(0) == true
  @test MAT.MAT_v4_Modelica.isLittleEndian(1000) == false
  @test MAT.MAT_v4_Modelica.isLittleEndian(2000) == false
  @test MAT.MAT_v4_Modelica.isLittleEndian(3000) == false
end

@testset "dataFormat" begin
  @test MAT.MAT_v4_Modelica.dataFormat(0000) <: Float64
  @test MAT.MAT_v4_Modelica.dataFormat(0020) <: Int32
  @test MAT.MAT_v4_Modelica.dataFormat(0030) <: Int16
  @test MAT.MAT_v4_Modelica.dataFormat(0040) <: UInt16
  @test MAT.MAT_v4_Modelica.dataFormat(0050) <: UInt8
end

@testset "typeBytes" begin
  @test MAT.MAT_v4_Modelica.typeBytes(Int32) == 4
end

@testset "Aclass" begin
  ac = MAT.MAT_v4_Modelica.readAclass(matBBO)
  @test ac.positionStart == 0
  @test ac.positionEnd == 71
end

@testset "readVariableNames" begin
  ac = MAT.MAT_v4_Modelica.readAclass(matBBO)
  vn = MAT.MAT_v4_Modelica.readVariableNames(ac)
  # @show vn
  @test length(vn.names) == 11
  @test vn.names[1] == "time"
  @test vn.names[3] == "vel"
  @test vn.names[11] == "grav"
  @test vn.positionStart == 71
  @test vn.positionEnd == 228
end

@testset "getVariableIndex" begin
  ac = MAT.MAT_v4_Modelica.readAclass(matBBO)
  vn = MAT.MAT_v4_Modelica.readVariableNames(ac)
  @test MAT.MAT_v4_Modelica.getVariableIndex(vn, vn.names[3]) == 3
  @test MAT.MAT_v4_Modelica.getVariableIndex(vn, vn.names[10]) == 10
end

@testset "readVariableDescriptions" begin
  ac = MAT.MAT_v4_Modelica.readAclass(matBBO)
  vn = MAT.MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT.MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  @test length(vd.descriptions) == 11  
  @test vd.descriptions[1] == "Simulation time [s]"
  @test vd.descriptions[3] == "velocity of ball"
  @test vd.descriptions[11] == "gravity acceleration"
end

@testset "readDataInfo" begin
  ac = MAT.MAT_v4_Modelica.readAclass(matBBO)
  vn = MAT.MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT.MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  di = MAT.MAT_v4_Modelica.readDataInfo(ac,vd)
  # @show di.info[3]
  @test di.info[1]["isWithinTimeRange"] == -1
  @test di.info[3]["locatedInData"] == 2 
  @test di.info[4]["isInterpolated"] == 0
  @test di.info[11]["isWithinTimeRange"] == 0
end

@testset "readVariable: BouncingBall" begin
  ac = MAT.MAT_v4_Modelica.readAclass(matBBO)
  vn = MAT.MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT.MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  di = MAT.MAT_v4_Modelica.readDataInfo(ac,vd)

  # println(JSON.json(di.info, 2)) #get the data 1/2 info

  eff = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "eff") #data1
  @test length(eff) == 2
  @test eff[1] ≈ 0.77
  @test eff[2] ≈ 0.77

  grav = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "grav") #data1
  @test length(grav) == 2
  @test grav[1] ≈ 9.81
  @test grav[2] ≈ 9.81

  time = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "time") # data0
  @test all(isapprox.(time, [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1], rtol=1e-3))

  height = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "height") #data2
  @test isapprox(height[1], 111, rtol=1e-3)
  @test isapprox(height[2], 110.9509, rtol=1e-3)

  vel = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "vel") #data2
  @test isapprox(vel[2], -0.981, rtol=1e-3)
end

@testset "readVariable: FallingBodyBox" begin
  ac = MAT.MAT_v4_Modelica.readAclass(matFBB)
  vn = MAT.MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT.MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  di = MAT.MAT_v4_Modelica.readDataInfo(ac,vd)

  # println(JSON.json(di.info, 2)) 

  var = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "time") 
  # display(var)
  ret = true
  for i = 2:length(var)-1 #last time is duplicated
    ret &= isapprox(var[i]-var[i-1], 0.002, rtol=1e-4)
  end
  @test ret == true

  #point-check values read from FallingBodyBox_res.csv
  var = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.frame_a.r_0[1]") 
  @test isapprox(var[16], 0.002923239, rtol=1e-3)

  var = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.frame_a.R.T[1,1]") 
  @test isapprox(var[26], 0.983794001, rtol=1e-3)

  @test_throws ErrorException MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.frame_a.v_0[1]")  # there is no frame_A.v_0

  var = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.v_0[2]") 
  @test isapprox(var[33], -0.58818129, rtol=1e-3)

  var = MAT.MAT_v4_Modelica.readVariable(ac, vn, vd, di, "bodyBox.frame_b.r_0[1]") 
  @test isapprox(var[72], 0.935886479, rtol=1e-3)
end

;

