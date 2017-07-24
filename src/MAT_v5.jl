# MAT_v5.jl
# Tools for reading MATLAB v5 files in Julia
#
# Copyright (C) 2012   Simon Kornblith
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

module MAT_v5
using Libz, BufferedStreams, HDF5
import Base: read, write, close
import HDF5: names, exists

if VERSION < v"0.6.0-dev.1632"
    round_uint8(data) = round(UInt8, data)
    complex_array(a, b) = complex(a, b)
else
    round_uint8(data) = round.(UInt8, data)
    complex_array(a, b) = complex.(a, b)
end

type Matlabv5File <: HDF5.DataFile
    ios::IOStream
    swap_bytes::Bool
    varnames::Dict{String, Int64}

    Matlabv5File(ios, swap_bytes) = new(ios, swap_bytes)
end

const miINT8 = 1
const miUINT8 = 2
const miINT16 = 3
const miUINT16 = 4
const miINT32 = 5
const miUINT32 = 6
const miSINGLE = 7
const miDOUBLE = 9
const miINT64 = 12
const miUINT64 = 13
const miMATRIX = 14
const miCOMPRESSED = 15
const miUTF8 = 16
const miUTF16 = 17
const miUTF32 = 18
const mxCELL_CLASS = 1
const mxSTRUCT_CLASS = 2
const mxOBJECT_CLASS = 3
const mxCHAR_CLASS = 4
const mxSPARSE_CLASS = 5
const mxDOUBLE_CLASS = 6
const mxSINGLE_CLASS = 7
const mxINT8_CLASS = 8
const mxUINT8_CLASS = 9
const mxINT16_CLASS = 10
const mxUINT16_CLASS = 11
const mxINT32_CLASS = 12
const mxUINT32_CLASS = 13
const mxINT64_CLASS = 14
const mxUINT64_CLASS = 15
const READ_TYPES = Type[
    Int8, UInt8, Int16, UInt16, Int32, UInt32, Float32, Union{},
    Float64, Union{}, Union{}, Int64, UInt64]
const CONVERT_TYPES = Type[
    Union{}, Union{}, Union{}, Union{},
    Union{}, Float64, Float32, Int8, UInt8,
    Int16, UInt16, Int32, UInt32, Int64, UInt64]

read_bswap{T}(f::IO, swap_bytes::Bool, ::Type{T}) =
    swap_bytes ? bswap(read(f, T)) : read(f, T)
function read_bswap{T}(f::IO, swap_bytes::Bool, ::Type{T}, dim::Union{Int, Tuple{Vararg{Int}}})
    d = read!(f, Array{T}(dim))
    if swap_bytes
        for i = 1:length(d)
            @inbounds d[i] = bswap(d[i])
        end
    end
    d
end

skip_padding(f::IO, nbytes::Int, hbytes::Int) = if nbytes % hbytes != 0
    skip(f, hbytes-(nbytes % hbytes))
end

# Read data type and number of bytes at the start of a data element
function read_header(f::IO, swap_bytes::Bool)
    dtype = read_bswap(f, swap_bytes, UInt32)

    if (dtype & 0xFFFF0000) != 0
        # Small Data Element Format
        (dtype & 0x0000FFFF, convert(Int, dtype >> 16), 4)
    else
        # Data Element Format
        (dtype, convert(Int, read_bswap(f, swap_bytes, UInt32)), 8)
    end
end

# Read data element as a vector of a given type
function read_element{T}(f::IO, swap_bytes::Bool, ::Type{T})
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    data = read_bswap(f, swap_bytes, T, Int(div(nbytes, sizeof(T))))
    skip_padding(f, nbytes, hbytes)
    data
end

# Read data element as encoded type
function read_data(f::IO, swap_bytes::Bool)
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    read_type = READ_TYPES[dtype]
    data = read_bswap(f, swap_bytes, read_type, Int(div(nbytes, sizeof(read_type))))
    skip_padding(f, nbytes, hbytes)
    data
end

# Read data element as encoded type with given dimensions, converting
# to another type if necessary and collapsing one-element matrices to
# scalars
function read_data{T}(f::IO, swap_bytes::Bool, ::Type{T}, dimensions::Vector{Int32})
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    read_type = READ_TYPES[dtype]

    read_array = any(dimensions .!= 1)
    if sizeof(read_type)*prod(dimensions) != nbytes
        error("Invalid element length")
    end
    if read_array
        data = read_bswap(f, swap_bytes, read_type, tuple(convert(Vector{Int}, dimensions)...))
    else
        data = read_bswap(f, swap_bytes, read_type)
    end
    skip_padding(f, nbytes, hbytes)

    read_array ? convert(Array{T}, data) : convert(T, data)
end

function read_cell(f::IO, swap_bytes::Bool, dimensions::Vector{Int32})
    data = Array{Any}(convert(Vector{Int}, dimensions)...)
    for i = 1:length(data)
        (ignored_name, data[i]) = read_matrix(f, swap_bytes)
    end
    data
end

function read_struct(f::IO, swap_bytes::Bool, dimensions::Vector{Int32}, is_object::Bool)
    if is_object
        class = String(read_element(f, swap_bytes, UInt8))
    end
    field_length = read_element(f, swap_bytes, Int32)[1]
    field_names = read_element(f, swap_bytes, UInt8)
    n_fields = div(length(field_names), field_length)

    # Get field names as strings
    field_name_strings = Vector{String}(n_fields)
    n_el = prod(dimensions)
    for i = 1:n_fields
        sname = field_names[(i-1)*field_length+1:i*field_length]
        index = findfirst(sname, 0)
        field_name_strings[i] = String(index == 0 ? sname : sname[1:index-1])
    end

    data = Dict{String, Any}()
    sizehint!(data, n_fields+1)
    if is_object
        data["class"] = class
    end

    if n_el == 1
        # Read a single struct into a dict
        for field_name in field_name_strings
            data[field_name] = read_matrix(f, swap_bytes)[2]
        end
    else
        # Read multiple structs into a dict of arrays
        for field_name in field_name_strings
            data[field_name] = Array{Any}(dimensions...)
        end
        for i = 1:n_el
            for field_name in field_name_strings
                data[field_name][i] = read_matrix(f, swap_bytes)[2]
            end
        end
    end

    data
end

function plusone!(A)
    for i = 1:length(A)
        @inbounds A[i] += 1
    end
    A
end

function read_sparse(f::IO, swap_bytes::Bool, dimensions::Vector{Int32}, flags::Vector{UInt32})
    local m::Int, n::Int
    if length(dimensions) == 2
        (m, n) = dimensions
    elseif length(dimensions) == 1
        m = dimensions[1]
        n = 1
    elseif length(dimensions) == 0
        m = 0
        n = 0
    else
        error("invalid dimensions encountered for sparse array")
    end

    m = isempty(dimensions) ? 0 : dimensions[1]
    n = length(dimensions) <= 1 ? 0 : dimensions[2]
    ir = plusone!(convert(Vector{Int}, read_element(f, swap_bytes, Int32)))
    jc = plusone!(convert(Vector{Int}, read_element(f, swap_bytes, Int32)))
    if (flags[1] & (1 << 9)) != 0 # logical
        # WTF. For some reason logical sparse matrices are tagged as doubles.
        pr = read_element(f, swap_bytes, Bool)
    else
        pr = read_data(f, swap_bytes)
        if (flags[1] & (1 << 11)) != 0 # complex
            pr = complex_array(pr, read_data(f, swap_bytes))
        end
    end

    SparseMatrixCSC(m, n, jc, ir, pr)
end

if VERSION >= v"0.4.0-dev+1039"
    truncate_to_uint8(x) = x % UInt8
else
    truncate_to_uint8(x) = convert(UInt8, x)
end

function read_string(f::IO, swap_bytes::Bool, dimensions::Vector{Int32})
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    if dtype <= 2 || dtype == 16
        # If dtype <= 2, this may give an error on non-ASCII characters, since the string
        # would be ISO-8859-1 and not UTF-8. However, MATLAB 2012b always saves strings with
        # a 2-byte encoding in v6 format, and saves UTF-8 in v7 format. Thus, this may never
        # happen in the wild.
        chars = read!(f, Vector{UInt8}(nbytes))
        if dimensions[1] <= 1
            data = String(chars)
        else
            data = Vector{String}(dimensions[1])
            for i = 1:dimensions[1]
                data[i] = rstrip(String(chars[i:dimensions[1]:end]))
            end
        end
    elseif dtype <= 4 || dtype == 17
        # Technically, if dtype == 3 or dtype == 4, this is ISO-8859-1 and not Unicode.
        # However, the first 256 Unicode code points are derived from ISO-8859-1, so UCS-2
        # is a superset of 2-byte ISO-8859-1.
        chars = read_bswap(f, swap_bytes, UInt16, convert(Int, div(nbytes, 2)))
        bufs = [IOBuffer() for i = 1:dimensions[1]]
        i = 1
        while i <= length(chars)
            for j = 1:dimensions[1]
                char = convert(Char, chars[i])
                if 255 < convert(UInt32, char)
                    # Newer versions of MATLAB seem to write some mongrel UTF-8...
                    char = String([truncate_to_uint8(chars[i] >> 8), truncate_to_uint8(chars[i])])[1]
                end
                write(bufs[j], char)
                i += 1
            end
        end

        if dimensions[1] == 0
            data = ""
        elseif dimensions[1] == 1
            data = String(take!(bufs[1]))
        else
            data = String[rstrip(String(take!(buf))) for buf in bufs]
        end
    else
        error("Unsupported string type")
    end
    skip_padding(f, nbytes, hbytes)
    data
end

# Read matrix data
function read_matrix(f::IO, swap_bytes::Bool)
    (dtype, nbytes) = read_header(f, swap_bytes)
    if dtype == miCOMPRESSED
        return read_matrix(ZlibInflateInputStream(read!(f, Vector{UInt8}(nbytes))), swap_bytes)
    elseif dtype != miMATRIX
        error("Unexpected data type")
    elseif nbytes == 0
        # If one creates a cell array using
        #     y = cell(m, n)
        # then MATLAB will save the empty cells as zero-byte matrices. If one creates a
        # empty cells using
        #     a = {[], [], []}
        # then MATLAB does not save the empty cells as zero-byte matrices. To avoid
        # surprises, we produce an empty array in both cases.
        return ("", Matrix{Union{}}(0, 0))
    end

    flags = read_element(f, swap_bytes, UInt32)
    dimensions = read_element(f, swap_bytes, Int32)
    name = String(read_element(f, swap_bytes, UInt8))

    class = flags[1] & 0xFF
    local data
    if class == mxCELL_CLASS
        data = read_cell(f, swap_bytes, dimensions)
    elseif class == mxSTRUCT_CLASS || class == mxOBJECT_CLASS
        data = read_struct(f, swap_bytes, dimensions, class == mxOBJECT_CLASS)
    elseif class == mxSPARSE_CLASS
        data = read_sparse(f, swap_bytes, dimensions, flags)
    elseif class == mxCHAR_CLASS && length(dimensions) <= 2
        data = read_string(f, swap_bytes, dimensions)
    else
        convert_type = CONVERT_TYPES[class]
        data = read_data(f, swap_bytes, convert_type, dimensions)
        if (flags[1] & (1 << 11)) != 0 # complex
            data = complex_array(data, read_data(f, swap_bytes, convert_type, dimensions))
        elseif (flags[1] & (1 << 9)) != 0 # logical
            data = reinterpret(Bool, round_uint8(data))
        end
    end

    return (name, data)
end

# Open MAT file for reading
matopen(ios::IOStream, endian_indicator::UInt16) =
    Matlabv5File(ios, endian_indicator == 0x494D)

# Read whole MAT file
function read(matfile::Matlabv5File)
    seek(matfile.ios, 128)
    vars = Dict{String, Any}()
    while !eof(matfile.ios)
        (name, data) = read_matrix(matfile.ios, matfile.swap_bytes)
        vars[name] = data
    end
    vars
end
# Read only variable names from an HDF5 file
function getvarnames(matfile::Matlabv5File)
    if !isdefined(matfile, :varnames)
        seek(matfile.ios, 128)
        matfile.varnames = varnames = Dict{String, Int64}()
        while !eof(matfile.ios)
            offset = position(matfile.ios)
            (dtype, nbytes, hbytes) = read_header(matfile.ios, matfile.swap_bytes)
            if dtype == miCOMPRESSED
                f = ZlibInflateInputStream(matfile.ios)
                read_header(f, matfile.swap_bytes)
            elseif dtype == miMATRIX
                f = matfile.ios
            else
                error("Unexpected data type")
            end

            read_element(f, matfile.swap_bytes, UInt32)
            read_element(f, matfile.swap_bytes, Int32)
            varnames[String(read_element(f, matfile.swap_bytes, UInt8))] = offset

            seek(matfile.ios, offset+nbytes+hbytes)
        end
    end
    matfile.varnames
end

exists(matfile::Matlabv5File, varname::String) =
    haskey(getvarnames(matfile), varname)
names(matfile::Matlabv5File) =
    keys(getvarnames(matfile))

# Read a variable from a MAT file
function read(matfile::Matlabv5File, varname::String)
    varnames = getvarnames(matfile)
    if !haskey(varnames, varname)
        error("no variable $varname in file")
    end
    seek(matfile.ios, varnames[varname])
    (name, data) = read_matrix(matfile.ios, matfile.swap_bytes)
    data
end

# Complain about writing to a MAT file
function write(parent::Matlabv5File, name::String, s)
    error("Writing to a MATLAB v5 file is not currently supported. Create a new file instead.")
end

# Close MAT file
close(matfile::Matlabv5File) = close(matfile.ios)
end
