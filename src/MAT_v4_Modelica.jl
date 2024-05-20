# A module to read MAT files written by OpenModelica tools

# Copyright (C) 2023   Ben Conrad
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Given the often large size of mat files, each read option is atomic, opening and closing the file for every operation with state managed by the user.

# The OpenModelica MATv4 file takes the basic v4 matrix format and adds some requirments on the contents and ordering of the matrices
# The format is described at https://openmodelica.org/doc/OpenModelicaUsersGuide/latest/technical_details.html#the-matv4-result-file-format, consisting of a series of matrices that first describe the data then store it.
#-  Aclass:
#    Aclass(1,:) is always Atrajectory
#    Aclass(2,:) is 1.1 in OpenModelica
#    Aclass(3,:) is empty
#    Aclass(4,:) is either binTrans or binNormal, which determines if the data is stored striped (rows of values at a time instance) or tansposed (rows being a single variable across time)
#- name:
#    a NxM matrix giving the names of the N variables as int8 characters
#- description:
#    a NxM matrix giving the descriptions of the N variables as int8 characters
#- dataInfo:
#    a Nx4 matrix describing the data of each variable, with
#    dataInfo(i,1) locating the data in data_1 or data_2
#    dataInfo(i,2) providing the start index within the data_ matrix
#    dataInfo(i,3) = 0 to indicate that the variable is interpolated
#    dataInfo(i,4) = -1 to indicate that the variable is undefined outside the time range
#- data_1:
#    is either an Nx1 matrix giving the variable's constant value, or Nx2 giving the start and end values
#- data_2:
#    holds the values of the continuously-varying variables in rows of [time1, var1(@time1), var2(@time1), ...varN(@time1), time2, var1(@time2)...]


module MAT_v4_Modelica

function isLittleEndian(dtype) :: Bool 
  #The type flag contains an integer whose decimal digits encode storage information. If the integer is represented as MOPT where M is the thousands digit...
  T, P, O, M = digits(dtype; pad=4)
  # @enum numberFormat little=0 big=1 vaxD=2 vaxG=3 cray=4
  return M == 0
end

function dataFormat(type) :: DataType 
  #The type flag contains an integer whose decimal digits encode storage information. If the integer is represented as MOPT where M is the thousands digit...
  T, P, O, M = digits(type; pad=4)
  # @enum dataFormat double=0 single=1 int32=2 int16=3 uint16=4 uint8=5
  if P == 0
    return Float64
  end
  if P == 1
    return Float32
  end
  if P == 2
    return Int32
  end
  if P == 3
    return Int16
  end
  if P == 4
    return UInt16
  end
  if P == 5
    return UInt8
  end
end

function typeBytes(type::T)::Int where T<:DataType
  if type == Float64
    return 8
  end
  if type == Float32
    return 4
  end
  if type == Int32
    return 4
  end
  if type == Int16
    return 2
  end
  if type == UInt16
    return 2
  end
  if type == UInt8
    return 1
  end
end


struct Aclass
  filepath::String
  isTranspose::Bool
  positionStart::Int
  positionEnd::Int
end

"""
Reads the Aclass matrix, returing if the data is stored binNormal or binTranspose
"""
function readAclass( filepath::String )
  open(filepath, "r") do matio
    seekstart(matio) # always start from the start, don't assume the position
    startP = position(matio)

    # The 20-byte header consists of five long (4-byte) integers:
    dtype = 0
    nrows = 0
    ncols = 0
    namelen = 0
    includesImaginary = 0
    try
      dtype, nrows, ncols, includesImaginary, namelen = read!(matio, Vector{Int32}(undef, 5)) # 32bit=4byte * 5 qty
    catch e
      error("caught error $e while reading $filepath")
    end

    if !isLittleEndian(dtype) 
      error("Only the little-endian encoding is implemented, cannot read $filepath")
    end

    nameuint = read!(matio, Vector{UInt8}(undef, namelen)) # read the full namelen to make the pointer ready to read the data
    name = strip(replace(String(nameuint), '\0'=>""))
    if name != "Aclass"
      error("First matrix must be named Aclass, is instead [$name]. This likely means that [$filepath] is not a Modelica MAT v4 file.")
    end

    # real: the Real part of the matrix consists of nrows ∗ ncols numbers in the format specified by the P element of the type flag. The data is stored column-wise such that the second column follows the first column, etc.
    fmt = dataFormat(dtype) # read the format type before reading
    realint = read!(matio, Matrix{UInt8}(undef, nrows,ncols))  

    Aclass1 = strip(replace(String(realint[1,:]), '\0'=>""))
    Aclass2 = strip(replace(String(realint[2,:]), '\0'=>""))
    Aclass3 = strip(replace(String(realint[3,:]), '\0'=>""))
    Aclass4 = strip(replace(String(realint[4,:]), '\0'=>""))
    if Aclass1 == "Atrajectory" && Aclass2 == "1.1" && isempty(Aclass3) && Aclass4 == "binNormal" || Aclass4 == "binTrans"
      return Aclass( filepath, Aclass4 == "binTranspose", startP, position(matio) )
    end
  end #open
end

"""
  Struct to record variable names and their start/stop file positions.
"""
struct VariableNames
  # names::Vector{T}(undef,undef) where T<:AbstractString
  names::Vector{String}
  positionStart::Int
  positionEnd::Int
end

function readVariableNames(ac::Aclass)
  open(ac.filepath, "r") do matio
    seek(matio, ac.positionEnd) #skip over Aclass
    startP = position(matio)

    #read the matrix header
    # The 20-byte header consists of five long (4-byte) integers:
    dtype = 0
    nrows = 0
    ncols = 0
    namelen = 0
    includesImaginary = 0
    try
      dtype, nrows, ncols, includesImaginary, namelen = read!(matio, Vector{Int32}(undef, 5)) # int32=4byte * 5 = 20byte
    catch e
      error("caught error $e while reading $ac.filepath")
    end

    #read the matrix name
    nameuint = read!(matio, Vector{UInt8}(undef, namelen)) # read the full namelen to make the pointer ready to read the data
    matrixName = strip(replace(String(nameuint), '\0'=>""))
    if matrixName != "name"
      error("trying to read matrix [name] but read [$matrixName]")
    end

    #read the matrix data
    fmt = dataFormat(dtype) # read the format type before reading
    realint = []
    try
      realint = read!(matio, Matrix{fmt}(undef, nrows,ncols))  # UInt8 from P is 8 bytes long
    catch e
      error("caught error $e while reading $ac.filepath")
    end

    #pull the names out of the matrix
    vnames = []
    for i in 1:ncols
      push!(vnames, strip(replace(String(realint[:,i]), '\0'=>""))) # note :,1 = implicit transpose
    end

    return VariableNames(vnames, startP, position(matio))

  end #open
end


function getVariableIndex(vn::VariableNames, name::String)
  vecAll = findall( x->x==name, vn.names)
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


struct VariableDescriptions
  names::Vector{String}
  descriptions::Vector{String}
  positionStart::Int
  positionEnd::Int
end

function readVariableDescriptions(ac::Aclass, vn::VariableNames)
   open(ac.filepath, "r") do matio
    seek(matio, vn.positionEnd) #this follows the VariableNames matrix
    startP = position(matio)

    # The 20-byte header consists of five long (4-byte) integers:
    dtype = 0
    nrows = 0
    ncols = 0
    namelen = 0
    includesImaginary = 0
    try
      dtype, nrows, ncols, includesImaginary, namelen = read!(matio, Vector{Int32}(undef, 5)) # int32=4byte * 5 = 20byte
    catch e
      error("caught error $e while reading $ac.filepath")
    end

    #read the matrix name
    nameuint = read!(matio, Vector{UInt8}(undef, namelen)) # read the full namelen to make the pointer ready to read the data
    matrixName = strip(replace(String(nameuint), '\0'=>""))
    if matrixName != "description"
      error("trying to read matrix [description] but read [$matrixName]")
    end

    #read the matrix data
    fmt = dataFormat(dtype) # read the format type before reading
    realread = []
    try
      realread = read!(matio, Matrix{fmt}(undef, nrows,ncols))  # UInt8 from P is 8 bytes long
    catch e
      error("caught error $e while reading $ac.filepath")
    end

    vdesc = []
    for i in 1:ncols
      push!(vdesc, strip(replace(String(realread[:,i]), '\0'=>""))) # note :,1 = implicit transpose
    end

    return VariableDescriptions(vn.names, vdesc, startP, position(matio))
  end #open
end


struct DataInfo
  info
  positionStart::Int
  positionEnd::Int
end

"""
dataInfo provides indicies to access variable data

dataInfo
Is an n x 4 integer matrix containing information for each variable (in the same order as the name and description matrices).
  dataInfo(i,1) is 1 or 2, saying if variable i is stored in the data_1 or data_2 matrix. If it is 0, it is the abscissa (time variable).
  dataInfo(i,2) contains the index in the data_1 or data_2 matrix. The index is 1-based and may contain several variables pointing to the same row (alias variables). A negative value means that the variable is a negated alias variable.
  dataInfo(i,3) is 0 to signify linear interpolation. In other tools the value is the number of times differentiable this variable is, which may improve plotting.
  dataInfo(i,4) is -1 in OpenModelica to signify that the value is not defined outside the time range. 0 keeps the first/last value when going outside the time range and 1 performs linear interpolation on the first/last two points.
"""
function readDataInfo(ac::Aclass, vd::VariableDescriptions)
   open(ac.filepath, "r") do matio
    seek(matio, vd.positionEnd) #this follows the VariableNames matrix
    startP = position(matio)

    # The 20-byte header consists of five long (4-byte) integers:
    dtype = 0
    nrows = 0
    ncols = 0
    namelen = 0
    includesImaginary = 0
    try
      dtype, nrows, ncols, includesImaginary, namelen = read!(matio, Vector{Int32}(undef, 5)) # int32=4byte * 5 = 20byte
    catch e
      error("caught error $e while reading $ac.filepath")
    end

    #read the matrix name
    nameuint = read!(matio, Vector{UInt8}(undef, namelen)) # read the full namelen to make the pointer ready to read the data
    matrixName = strip(replace(String(nameuint), '\0'=>""))
    if matrixName != "dataInfo"
      error("trying to read variable [dataInfo] but read [$matrixName]")
    end

    #read the matrix data
    fmt = dataFormat(dtype) # read the format type before reading
    realread = []
    try
      realread = read!(matio, Matrix{fmt}(undef, nrows,ncols))  # UInt8 from P is 8 bytes long
    catch e
      error("caught error $e while reading $ac.filepath")
    end

    dinfo = []
    for i in 1:ncols
      push!(dinfo, Dict("name"=>vd.names[i], "description"=>vd.descriptions[i], "locatedInData"=>realread[1,i], "indexInData"=>realread[2,i], "isInterpolated"=>realread[3,i], "isWithinTimeRange"=>realread[4,i] ) )
    end
    return DataInfo( dinfo, startP, position(matio))
  end #open
end

struct MatrixHeader
  type::Int
  nRows::Int
  nCols::Int
  hasImaginary::Bool
  lName::Int
  name::String
  format::DataType
end

"""
Reads the matix header, assuming matio's position is correct to read the header
"""
function readMatrixHeader!(matio::IOStream) :: MatrixHeader
  dtype = 0
  nrows = 0
  ncols = 0
  namelen = 0
  includesImaginary = 0
  try
    dtype, nrows, ncols, includesImaginary, namelen = read!(matio, Vector{Int32}(undef, 5)) # int32=4byte * 5 = 20byte
  catch e
    error("caught error $e while reading matrix header")
  end

  # data1MatrixName = mark(matio)
  nameuint = read!(matio, Vector{UInt8}(undef, namelen)) # read the full namelen to make the pointer ready to read the data
  matrixName = strip(replace(String(nameuint), '\0'=>""))

  fmt = dataFormat(dtype) # read the format type before reading

  return MatrixHeader(dtype,nrows,ncols,includesImaginary,namelen,matrixName, fmt)
end

"""
read one variable from the file
"""
function readVariable(ac::Aclass, vn::VariableNames, vd::VariableDescriptions, di::DataInfo, name::String)
  open(ac.filepath, "r") do matio
    seek(matio, di.positionEnd) #this follows the VariableNames matrix

    # read data1 header:
    mh1 = readMatrixHeader!(matio)
    @assert mh1.name == "data_1" "trying to read matrix [data_1] but read $matrixName"

    data1MatrixStart = mark(matio)
    try
      toskip = mh1.nRows*mh1.nCols*typeBytes(mh1.format) #817*2*8 = 13072, 488197+20+7+13072 = 501296
      skip(matio, toskip )
    catch e
      throw( ErrorException("Caught error $e while reading $ac.filepath") )
    end

    # read data2 header:
    mh2 = readMatrixHeader!(matio)

    @assert mh2.name == "data_2" "trying to read matrix [data_2] but read $(mh2.name)"
    data2MatrixStart = mark(matio)

    #with the positions marked, read the desired variable
    varInd = getVariableIndex(vn, name)
    if varInd < 1
      throw( ArgumentError("Variable [$name] not found in file [$(ac.filepath)]") )
    end

    if di.info[varInd]["locatedInData"] == 1 #data_1
      #read the matrix data_1
      # dataInfo(i,4) is -1 in OpenModelica to signify that the value is not defined outside the time range. 0 keeps the first/last value when going outside the time range and 1 performs linear interpolation on the first/last two points.
      if di.info[varInd]["isWithinTimeRange"]== 0 #linear interpolation
        #data format is: time(tInitial), var1(tI), ... varN(tI), time(tFinal), var1(tF), ... varN(tF)
        readns = Vector{mh1.format}(undef, mh1.nCols)
        for ind = 1:mh1.nCols
          seek(matio, data1MatrixStart + (di.info[varInd]["indexInData"]-1)*typeBytes(mh1.format) + ((ind-1)*mh1.nRows*typeBytes(mh1.format)) )
          readns[ind] = read(matio, mh1.format)
        end
        return readns
      end

    elseif name == "Time" || name == "time" || di.info[varInd]["locatedInData"] == 2 #data_2
      #read the matrix data_2
      if ac.isTranspose == false
        # data is sequential: time(t0), var1(t0), var2(t0),... varN(t0), time(t1), var1(t1),...
        readns = Vector{mh2.format}(undef, mh2.nCols)
        for ind = 1:mh2.nCols
          seek(matio, data2MatrixStart + (di.info[varInd]["indexInData"]-1)*typeBytes(mh2.format) + ((ind-1)*mh2.nRows*typeBytes(mh2.format)) )
          readns[ind] = read(matio, mh2.format)
        end
        return readns
      else
        throw(ErrorException("reading binTranspose not implemented, lack test data") )
      end
    else
      throw(ErrorException("variable [$name] is located in an unknown location") )
    end
  end #open
end

"""
Determine whether the `filepath` is a mat file in the Modelica format.
"""
function isMatV4Modelica(filepath::String)
  #first check the extension
  ret = splitext(filepath)[2] == ".mat" || splitext(filepath)[2] == ".MAT"

  #now we need to interrogate the contents, by ensuring that all of the internal functions return correctly
  ac = readAclass(filepath)
  vn = readVariableNames(ac)
  vd = readVariableDescriptions(ac,vn)
  di = readDataInfo(ac,vd)

  ret &= typeof(ac) == Aclass
  ret &= typeof(vn) == VariableNames
  ret &= typeof(vd) == VariableDescriptions
  ret &= typeof(di) == DataInfo
  return ret;
end


"""
All-in-one reading of variable `name` from `filepath`, returning a Dict with keys "time" and `name`
"""
function readVariable(filepath::String, name::String) :: Dict
  ac = readAclass(filepath)
  vn = readVariableNames(ac)
  vd = readVariableDescriptions(ac,vn)
  di = readDataInfo(ac,vd)

  time = readVariable(ac, vn, vd, di, "time") 
  varn = readVariable(ac, vn, vd, di, name) 
  return Dict(["time"=>time, name=>varn]);
end

"""
Reads the vector of variable `names` from mat file `filepath`, returning a Dict with columns "time" and `names`
"""
# function readVariables(filepath::String, names::AbstractVector{T}) :: Dict where T<:AbstractString 
function readVariables(filepath::String, names::AbstractVector{T}) where T<:AbstractString 
  ac = readAclass(filepath)
  vn = readVariableNames(ac)
  vd = readVariableDescriptions(ac,vn)
  di = readDataInfo(ac,vd)

  data = Dict(["time"=> readVariable(ac, vn, vd, di, "time")])
  for name in names
    var = readVariable(ac,vn,vd,di, name)
    if length(var) == 1 # a constant value
      # data[name] = var[1]
      data[name] = [var[1]]
    elseif length(var) == 2 && var[1] == var[2] # this is a 2-element constant array: [123, 123]
      # data[name] = var[1]
      data[name] = [var[1]] # Dict infers Vector64
    elseif length(var) == length(data["time"]) # a regular vector
      data[name] = var
    else
      throw(DimensionMismatch("Length of $name [$(length(var))] differs from the dataframe [$(length(data["time"]))], cannot add it"))
    end
  end
  return data
end



end #MAT_v4_Modelica
