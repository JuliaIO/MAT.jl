cd(joinpath(@__DIR__,".."))
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# test only MAT_v4
using Test
# include("../src/MAT.jl")
# include("../src/MAT_v4_Modelica.jl")

mat1s = "W:/sync/mechgits/library/julia/ConvenientModelica/test/ViseHammer_result_1s/ViseHammer_res.mat"


@testset "manual" begin
  rawfid = open(mat1s, "r")

  # This reader is collected from two references:
  # (isv4, swap_bytes) = MAT_v4.checkv4(rawfid)
  #  Matlab_MATFileFormatV4.pdf from https://www.eiscat.se/wp-content/uploads/2016/03/Version-4-MAT-File-Format.pdf
  #  OpenModelica_TheMATv4ResultFileFormat_OpenModelicaUserGuide_v1.21.0-dev-211-gcd8969a225.pdf from https://www.openmodelica.org/doc/OpenModelicaUsersGuide/latest/technical_details.html
  
  # Start by reading the header, following Matlab:
  # Each matrix starts with a fixed-length 20-byte header that contains information describing certain attributes of the Matrix.
  # The 20-byte header consists of five long (4-byte) integers:
  @show type, nrows, ncols, imagf, namelen = read!(rawfid, Vector{Int32}(undef, 5)) # int32=4byte * 5 = 20byte

  #The type flag contains an integer whose decimal digits encode storage information. If the integer is represented as MOPT where M is the thousands digit...
  T, P, O, M = digits(type; pad=4)

  # M indicates the numeric format of binary numbers on the machine that wrote the file.
  # Use this table to determine the number to use for your machine:
  #   0 IEEE Little Endian (PC, 386, 486, DEC Risc)
  #   1 IEEE Big Endian (Macintosh, SPARC, Apollo,SGI, HP 9000/300, other Motorola)
  #   2 VAX D-float
  #   3 VAX G-float
  #   4 Cray
  @enum numberFormat little=0 big=1 vaxD=2 vaxG=3 cray=4
  @show numberFormat(M)

  # O is always 0 (zero) and is reserved for future use.
  # @show O

  # P indicates which format the data is stored:
  #   0 double-precision (64-bit) floating point numbers
  #   1 single-precision (32-bit) floating point numbers
  #   2 32-bit signed integers
  #   3 16-bit signed integers
  #   4 16-bit unsigned integers
  #   5 8-bit unsigned integers
  @enum dataFormat double=0 single=1 int32=2 int16=3 uint16=4 uint8=5
  @show dataFormat(P)
  @test dataFormat(P) == uint8

  # T indicates matrix type: 1 = text matrix
  @enum matrixFormat numericFull=0 text=1 spars=2
  @show matrixFormat(T)

  # nrows The row dimension contains an integer with the number of rows in the matrix.
  @show nrows
  # ncols The column dimension contains an integer with the number of columns in the matrix.
  @show ncols
  # imagf The imaginary flag is an integer whose value is either 0 or 1. If 1, then the matrix has an imaginary part. If 0, there is only real data.
  @show imagf
  # namelen The name length contains an integer with 1 plus the length of the matrix name.
  @show namelen

  # Immediately following the fixed length header is the data whose length is dependent on the variables in the fixed length header:
  #  name The matrix name consists of namlen ASCII bytes, the last one of which must be a null character (’\0’ ).
  @show nameuint = read!(rawfid, Vector{UInt8}(undef, namelen)) # read the full namelen to make the pointer ready to read the data
  # pop!(nameuint) #trim \0
  # @show name = String(nameuint)
  @show name = replace(String(nameuint), '\0'=>"")

  #  real Real part of the matrix consists of nrows ∗ ncols numbers in the format specified by the P element of the type flag. The data is stored column-wise such that the second column follows the first column, etc.
  # @show realint = read!(rawfid, Vector{UInt8}(undef, 1*nrows*ncols))  # UInt8 from P is 8 bytes long
  # @show realint = read!(rawfid, Vector{UInt16}(undef, 8*nrows*ncols))  # UInt8 from P is 8 bytes long
  # @show realint = read(rawfid, 32)  # UInt8 from P is 8 bytes long
  @show realint = read!(rawfid, Matrix{UInt8}(undef, nrows,ncols))  # UInt8 from P is 8 bytes long

  #  imag Imaginary part of the matrix, if any. If the imaginary flag imagf is nonzero, the imaginary part of a matrix is placed here. It is stored in the same manner as the real data. 
  if imagf==1
    @show imagint = read!(rawfid, Vector{UInt8}(undef, nrows*ncols)) # read the full namelen to make the pointer ready to read the data
  end

  # The variables stored in the MAT-file are (in the order required by OpenModelica):
  #   Aclass
  # • Aclass(1,:) is always Atrajectory
  # • Aclass(2,:) is 1.1 in OpenModelica
  # • Aclass(3,:) is empty
  # • Aclass(4,:) is either binTrans or binNormal

  Aclass1 = replace(String(realint[1,:]), '\0'=>"")
  Aclass2 = replace(String(realint[2,:]), '\0'=>"")
  Aclass3 = replace(String(realint[3,:]), '\0'=>"")
  Aclass4 = replace(String(realint[4,:]), '\0'=>"")
  @test Aclass1 == "Atrajectory"
  @test Aclass2 == "1.1"
  @test isempty(Aclass3)
  @test Aclass4 == "binNormal" || Aclass4 == "binTrans"

  #now to read the next matrix header...
  # name:
  # Is an n x m character (int8) matrix, where n is the number of variables stored in the result-file (including time). m is the length of the longest variable. OpenModelica stores the trailing part of the name as NIL bytes (0) whereas other tools use spaces for the trailing part.
  @show type, nrows, ncols, imagf, namelen = read!(rawfid, Vector{Int32}(undef, 5)) # int32=4byte * 5 = 20byte
  T, P, O, M = digits(type; pad=4)
  @show numberFormat(M)
  @show matrixFormat(T)
  @show dataFormat(P)

  @show nrows
  @show ncols
  @show imagf
  @show namelen
  nameuint = read!(rawfid, Vector{UInt8}(undef, namelen)) # read the full namelen to make the pointer ready to read the data
  @show name = replace(String(nameuint), '\0'=>"")

  realint = read!(rawfid, Matrix{UInt8}(undef, nrows,ncols))  # UInt8 from P is 8 bytes long
  @show varname1 = replace(String(realint[:,1]), '\0'=>"") # note :,1 implicit transpose
  @show varname2 = replace(String(realint[:,2]), '\0'=>"") # note :,1 implicit transpose
  @show varname3 = replace(String(realint[:,3]), '\0'=>"") # note :,1 implicit transpose
  @show varname4 = replace(String(realint[:,4]), '\0'=>"") # note :,1 implicit transpose
  @show varname5 = replace(String(realint[:,5]), '\0'=>"") # note :,1 implicit transpose

  close(rawfid)
end


