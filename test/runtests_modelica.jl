cd(joinpath(@__DIR__,".."))
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# test only MAT_v4
using Test
# include("../src/MAT.jl")
include("../src/MAT_v4_Modelica.jl")

const mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"



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
  mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat1s)
  @test ac.positionStart == 0
  @test ac.positionEnd == 71
end


@testset "readVariableNames" begin
  mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat1s)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  # @show vn
  @test length(vn.names) == 2490
  @test vn.names[1] == "time"
  @test vn.names[3] == "revolute.w"
  @test vn.names[30] == "der(alignElastoBacklash.frame_a.r_0[1])"
  @test vn.names[2490] == "world.z_label.color[3]"
  @test vn.positionStart == 71
  @test vn.positionEnd == 117126
end


@testset "getVariableIndex" begin
  mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat1s)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  @test MAT_v4_Modelica.getVariableIndex(vn, vn.names[3]) == 3
  @test MAT_v4_Modelica.getVariableIndex(vn, vn.names[30]) == 30
end


@testset "readVariableDescriptions" begin
  mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat1s)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  @test length(vd.descriptions) == 2490  
  @test vd.descriptions[1] == "Simulation time [s]"
  @test vd.descriptions[3] == "First derivative of angle phi (relative angular velocity) [rad/s]"
  @test vd.descriptions[30] == "Position vector from world frame to the connector frame origin, resolved in world frame"
  @test vd.descriptions[2490] == "Color of cylinders"
end



@testset "readDataInfo" begin
  mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"
  ac = MAT_v4_Modelica.readAclass(mat1s)
  vn = MAT_v4_Modelica.readVariableNames(ac)
  vd = MAT_v4_Modelica.readVariableDescriptions(ac,vn)
  di = MAT_v4_Modelica.readDataInfo(ac,vd)
  # @show di.info[3]
  @test di.info[1]["isWithinTimeRange"] == -1
  @test di.info[3]["locatedInData"] == 2 
  @test di.info[30]["isInterpolated"] == 0
  @test di.info[2490]["isWithinTimeRange"] == 0
end

using JSON
@testset "readVariable" begin
  # mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"
  matbb = "W:/sync/mechgits/library/julia/MAT.jl/test/Modelica/BouncingBall/BouncingBall_res.mat"
  ac = MAT_v4_Modelica.readAclass(matbb)
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


;

