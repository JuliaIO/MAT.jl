# Copied from pr132 of https://github.com/JuliaIO/MAT.jl/pull/132

# MAT_v4.jl
# Tools for reading MATLAB v4 files in Julia
#
# Copyright (C) 2012   Simon Kornblith
# Copyright (C) 2019   Victor Saase (modified from MAT_v5.jl)
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

# MATLAB's file format documentation can be found at
# http://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf

module MAT_v4_pr132
using BufferedStreams, HDF5, SparseArrays
import Base: read, write, close
import HDF5: names

round_uint8(data) = round.(UInt8, data)
complex_array(a, b) = complex.(a, b)

mutable struct Matlabv4File <: HDF5.H5DataStore
    ios::IOStream
    swap_bytes::Bool
    varnames::Dict{String, Int64}

    Matlabv4File(ios, swap_bytes) = new(ios, swap_bytes)
end

const mLITTLE_ENDIAN = 0
const mBIG_ENDIAN = 1
const mVAX_DFLOAT = 2
const mVAX_GFLOAT = 3
const mGRAY = 4

const pTYPE = Dict(
    0 => Float64,
    1 => Float32,
    2 => Int32,
    3 => Int16,
    4 => UInt16,
    5 => UInt8
)

const tNUMERIC = 0
const tTEXT = 1
const tSPARSE = 2

const imagfREAL = 0
const imagfCOMPLEX = 1

read_bswap(f::IO, swap_bytes::Bool, ::Type{T}) where T = swap_bytes ? bswap(read(f, T)) : read(f, T)

function read_bswap(f::IO, swap_bytes::Bool, ::Type{T}, dim::Union{Int, Tuple{Vararg{Int}}}) where T
    d = read!(f, Array{T}(undef, dim))
    if swap_bytes
        for i = 1:length(d)
            @inbounds d[i] = bswap(d[i])
        end
    end
    d
end

function read_bswap(f::IO, swap_bytes::Bool, d::AbstractArray{T}) where T
    readbytes!(f, reinterpret(UInt8, d))
    if swap_bytes
        for i = 1:length(d)
            @inbounds d[i] = bswap(d[i])
        end
    end
    d
end

function checkv4(f::IO)
    M, O, P, T, mrows, ncols, imagf, namlen = MAT_v4.read_header(f, false)
    if 0<=M<=4 && O == 0 && 0<=P<=5 && 0<=T<=2 && mrows>=0 && ncols>=0 && 0<=imagf<=1 && namlen>0
        swap_bytes = false
        return (true, swap_bytes)
    else
        seek(f, 0)   
        M, O, P, T, mrows, ncols, imagf, namlen = MAT_v4.read_header(f, true)
        if 0<=M<=4 && O == 0 && 0<=P<=5 && 0<=T<=2 && mrows>=0 && ncols>=0 && 0<=imagf<=1 && namlen>0
            swap_bytes = true
            return (true, swap_bytes)
        end
    end
    return (false, false)
end

# Read data type and number of bytes at the start of a data element
function read_header(f::IO, swap_bytes::Bool)
    dtype = read_bswap(f, swap_bytes, Int32)

    M = div(rem(dtype, 10000), 1000)
    O = div(rem(dtype, 1000), 100)
    P = div(rem(dtype, 100), 10)
    T = div(rem(dtype, 10), 1)

    mrows = read_bswap(f, swap_bytes, Int32)
    ncols = read_bswap(f, swap_bytes, Int32)
    imagf = read_bswap(f, swap_bytes, Int32)
    namlen = read_bswap(f, swap_bytes, Int32)
    
    M, O, P, T, mrows, ncols, imagf, namlen
end

# Read matrix data
function read_matrix(f::IO, swap_bytes::Bool)
    M, O, P, T, mrows, ncols, imagf, namlen = read_header(f, swap_bytes)
    if ncols == 0 || mrows == 0
        # If one creates a cell array using
        #     y = cell(m, n)
        # then MATLAB will save the empty cells as zero-byte matrices. If one creates a
        # empty cells using
        #     a = {[], [], []}
        # then MATLAB does not save the empty cells as zero-byte matrices. To avoid
        # surprises, we produce an empty array in both cases.
        return ("", Matrix{Union{}}(undef, 0, 0))
    end
    name = String(read_bswap(f, M==mBIG_ENDIAN, Vector{UInt8}(undef, namlen))[1:end-1])
    if T == tNUMERIC || T == tSPARSE
        real_data = read_bswap(f, M==mBIG_ENDIAN, Vector{pTYPE[P]}(undef, ncols*mrows))
        if T == tNUMERIC && imagf == imagfCOMPLEX
            imag_data = read_bswap(f, M==mBIG_ENDIAN, Vector{pTYPE[P]}(undef, ncols*mrows))
            data = complex.(real_data, imag_data)
        else
            data = real_data
        end
        datamat = reshape(data, Int(mrows), Int(ncols))
        if T == tNUMERIC
            return (name, datamat)
        elseif T == tSPARSE
            if size(datamat,2) == 3
                return (name, sparse(datamat[:,1], datamat[:,2], datamat[:,3]))
            else
                return (name, sparse(datamat[:,1], datamat[:,2], complex.(datamat[:,3], datamat[:,4])))
            end
        end
    elseif T == tTEXT
        if mrows > 1
            charvec = UInt8.(read_bswap(f, M==mBIG_ENDIAN, Vector{pTYPE[P]}(undef, ncols*mrows)))
            charmat = reshape(charvec, Int(mrows), Int(ncols))
            data = [String(charmat[i,:]) for i in 1:mrows]
        else
            data = String(UInt8.(read_bswap(f, M==mBIG_ENDIAN, Vector{pTYPE[P]}(undef, ncols))))
        end
        return (name, data)
    end
end

# Open MAT file for reading
matopen(ios::IOStream, endian_indicator::Bool) = Matlabv4File(ios, endian_indicator)

# Read whole MAT file
function read(matfile::Matlabv4File)
    seek(matfile.ios, 0)
    vars = Dict{String, Any}()
    while !eof(matfile.ios)
        (name, data) = read_matrix(matfile.ios, matfile.swap_bytes)
        vars[name] = data
    end
    vars
end

# Read only variable names
function getvarnames(matfile::Matlabv4File)
    if !isdefined(matfile, :varnames)
        seek(matfile.ios, 0)
        matfile.varnames = varnames = Dict{String, Int64}()
        while !eof(matfile.ios)
            offset = position(matfile.ios)
            M, O, P, T, mrows, ncols, imagf, namlen = read_header(matfile.ios, matfile.swap_bytes)
            f = matfile.ios
            
            name = String(read_bswap(f, M==mBIG_ENDIAN, Vector{UInt8}(undef, namlen))[1:end-1])
            varnames[name] = offset
            imag_offset = 0
            skip(f, mrows*ncols*sizeof(pTYPE[P]))
            if imagf == imagfCOMPLEX
                skip(f, mrows*ncols*sizeof(pTYPE[P]))
            end
        end
    end
    matfile.varnames
end

names(matfile::Matlabv4File) = keys(getvarnames(matfile))

# Read a variable from a MAT file
function read(matfile::Matlabv4File, varname::String)
    varnames = getvarnames(matfile)
    if !haskey(varnames, varname)
        error("no variable $varname in file")
    end
    seek(matfile.ios, varnames[varname])
    (name, data) = read_matrix(matfile.ios, matfile.swap_bytes)
    data
end

function colvals(A::AbstractSparseMatrix)
    rows = rowvals(A)
    cols = similar(rows)
    m,n = size(A)
    for i=1:n
        for j in nzrange(A,i)
            cols[j] = i
        end
    end
    cols
end

function write(parent::Matlabv4File, name::String, s)
    M = Int(parent.swap_bytes)
    O = 0
    P = 0
    for p=keys(pTYPE)
        if eltype(s) == pTYPE[p] || eltype(s) == Complex{pTYPE[p]}
            P = p
        end
    end
    if pTYPE[P] != eltype(s) && Complex{pTYPE[P]} != eltype(s) && 
        !(s isa AbstractString && pTYPE[P] == Float64) &&
        !(s isa Vector{String} && pTYPE[P] == Float64)
        error("invalid value type when writing v4 file")
    end
    if s isa AbstractSparseMatrix
        T = tSPARSE
    elseif s isa AbstractString || s isa Vector{String}
        T = tTEXT
    else
        T = tNUMERIC
    end
    write(parent.ios, Int32(1000*M + 100*O + 10*P + T))

    mrows = 1
    ncols = 1
    if s isa AbstractVector && !(s isa Vector{String})
        ncols = length(s)
    elseif s isa AbstractMatrix
        if s isa AbstractSparseMatrix
            mrows = nnz(s)
            ncols = 3
            if eltype(s) <: Complex
                ncols = 4
            end
        else
            mrows, ncols = size(s)
        end
    elseif s isa Vector{String}
        ncols = length(s[1])
        mrows = length(s)
    elseif s isa AbstractString
        ncols = length(s)
    end
    write(parent.ios, Int32(mrows))
    write(parent.ios, Int32(ncols))
    
    imagf = 0
    if eltype(s) <: Complex && T == tNUMERIC
        imagf = 1
    end
    write(parent.ios, Int32(imagf))

    namlen = length(name) + 1
    write(parent.ios, Int32(namlen))

    write(parent.ios, Vector{UInt8}(name))
    write(parent.ios, UInt8(0))

    if s isa AbstractArray && T == tNUMERIC
        write(parent.ios, reshape(real(s), length(s)))
        if imagf == 1
            write(parent.ios, reshape(imag(s), length(s)))
        end
    elseif s isa Number
        write(parent.ios, real(s))
        if imagf == 1
            write(parent.ios, imag(s))
        end
    elseif s isa AbstractString
        floatarray = Float64.(Vector{UInt8}(s))
        write(parent.ios, floatarray)
    elseif s isa Vector{String}
        floatarray = Matrix{Float64}(undef, mrows, ncols)
        for (i,strel) = enumerate(s)
            floatarray[i,:] = Float64.(Vector{UInt8}(strel))
        end
        write(parent.ios, floatarray)
    elseif T == tSPARSE
        rows = rowvals(s)
        cols = colvals(s)
        vals = nonzeros(s)
        write(parent.ios, pTYPE[P].(rows))
        write(parent.ios, pTYPE[P].(cols))
        write(parent.ios, pTYPE[P].(real(vals)))
        if eltype(s) <: Complex
            write(parent.ios, pTYPE[P].(imag(vals)))
        end
    end
end

# Close MAT file
close(matfile::Matlabv4File) = close(matfile.ios)

end
