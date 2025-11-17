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
using CodecZlib, BufferedStreams, HDF5, SparseArrays
import Base: read, write, close
import ..MAT_types: MatlabStructArray, MatlabClassObject

using ..MAT_subsys

round_uint8(data) = round.(UInt8, data)
complex_array(a, b) = complex.(a, b)

mutable struct Matlabv5File <: HDF5.H5DataStore
    ios::IOStream
    swap_bytes::Bool
    subsystem::Subsystem
    subsystem_position::UInt64 # nr of bytes taken by subsystem
    varnames::Dict{String, Int64}

    Matlabv5File(ios, swap_bytes) = new(ios, swap_bytes, Subsystem(), UInt64(0))
end

function Base.show(io::IO, f::Matlabv5File)
    print(io, "Matlabv5File(")
    print(io, f.ios, ", ")
    print(io, f.swap_bytes, ")")
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
const mxFUNCTION_CLASS = 16
const mxOPAQUE_CLASS = 17
const READ_TYPES = Type[
    Int8, UInt8, Int16, UInt16, Int32, UInt32, Float32, Union{},
    Float64, Union{}, Union{}, Int64, UInt64]
const CONVERT_TYPES = Type[
    Union{}, Union{}, Union{}, Union{},
    Union{}, Float64, Float32, Int8, UInt8,
    Int16, UInt16, Int32, UInt32, Int64, UInt64]

read_bswap(f::IO, swap_bytes::Bool, ::Type{T}) where T =
    swap_bytes ? bswap(read(f, T)) : read(f, T)
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
function read_element(f::IO, swap_bytes::Bool, ::Type{T}) where T
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
function read_data(f::IO, swap_bytes::Bool, ::Type{T}, dimensions::Vector{Int32}) where T
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    read_type = READ_TYPES[dtype]
    if (read_type === UInt8) && (T === Bool)
        read_type = Bool
    end

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

function read_cell(f::IO, swap_bytes::Bool, dimensions::Vector{Int32}, subsys::Subsystem)
    data = Array{Any}(undef, convert(Vector{Int}, dimensions)...)
    for i = 1:length(data)
        (ignored_name, data[i]) = read_matrix(f, swap_bytes, subsys)
    end
    data
end

function read_struct(f::IO, swap_bytes::Bool, dimensions::Vector{Int32}, is_object::Bool, subsys::Subsystem)
    if is_object
        class = String(read_element(f, swap_bytes, UInt8))
    else
        class = ""
    end
    field_length = read_element(f, swap_bytes, Int32)[1]
    field_names = read_element(f, swap_bytes, UInt8)
    n_fields = div(length(field_names), field_length)

    # Get field names as strings
    field_name_strings = Vector{String}(undef, n_fields)
    n_el = prod(dimensions)
    for i = 1:n_fields
        sname = field_names[(i-1)*field_length+1:i*field_length]
        index = something(findfirst(iszero, sname), 0)
        field_name_strings[i] = String(index == 0 ? sname : sname[1:index-1])
    end

    local data
    if n_el == 1
        # Read a single struct into a dict
        data = Dict{String, Any}()
        sizehint!(data, n_fields+1)
        for field_name in field_name_strings
            data[field_name] = read_matrix(f, swap_bytes, subsys)[2]
        end
        if is_object
            data = MatlabClassObject(data, class)
        end
    else
        # Read empty or multiple structs
        nfields = length(field_name_strings)
        N = length(dimensions)
        field_values = Array{Any, N}[Array{Any}(undef, dimensions...) for _ in 1:nfields]
        for i = 1:n_el
            for field in 1:nfields
                field_values[field][i] = read_matrix(f, swap_bytes, subsys)[2]
            end
        end
        data = MatlabStructArray{N}(field_name_strings, field_values, class)
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
    if length(ir) > length(pr)
        # Fix for Issue #169, xref https://github.com/JuliaLang/julia/pull/40523
        #=
        # The following expression must be obeyed according to
        # https://github.com/JuliaLang/julia/blob/b3e4341d43da32f4ab6087230d98d00b89c8c004/stdlib/SparseArrays/src/sparsematrix.jl#L86-L90
        @debug "SparseMatrixCSC" m n jc ir pr
        @debug "SparseMatrixCSC check " length(jc) n+1 jc[end]-1 length(ir) length(pr) begin
            length(jc) == n + 1 && jc[end] - 1 == length(ir) == length(pr)
        end
        =#
        # Truncate rowvals (ir) to the be the same length as the non-zero elements (pr)
        resize!(ir, length(pr))
    end
    SparseMatrixCSC(m, n, jc, ir, pr)
end

truncate_to_uint8(x) = x % UInt8

function read_string(f::IO, swap_bytes::Bool, dimensions::Vector{Int32})
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    if dtype <= 2 || dtype == 16
        # If dtype <= 2, this may give an error on non-ASCII characters, since the string
        # would be ISO-8859-1 and not UTF-8. However, MATLAB 2012b always saves strings with
        # a 2-byte encoding in v6 format, and saves UTF-8 in v7 format. Thus, this may never
        # happen in the wild.
        chars = read!(f, Vector{UInt8}(undef, nbytes))
        if dimensions[1] <= 1
            data = String(chars)
        else
            data = Vector{String}(undef, dimensions[1])
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

function read_opaque(f::IO, swap_bytes::Bool, subsys::Subsystem)
    type_name = String(read_element(f, swap_bytes, UInt8))
    classname = String(read_element(f, swap_bytes, UInt8))

    if classname == "FileWrapper__"
        return read_matrix(f, swap_bytes, subsys)
    end

    _, metadata = read_matrix(f, swap_bytes, subsys)
    return MAT_subsys.load_mcos_object(metadata, type_name, subsys)
end

# Read matrix data
function read_matrix(f::IO, swap_bytes::Bool, subsys::Subsystem)
    (dtype, nbytes) = read_header(f, swap_bytes)
    if dtype == miCOMPRESSED
        decompressed_ios = ZlibDecompressorStream(IOBuffer(read!(f, Vector{UInt8}(undef, nbytes))))
        return read_matrix(decompressed_ios, swap_bytes, subsys)
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
        return ("", Matrix{Union{}}(undef, 0, 0))
    end

    flags = read_element(f, swap_bytes, UInt32)
    class = flags[1] & 0xFF

    if class != mxOPAQUE_CLASS
        dimensions = read_element(f, swap_bytes, Int32)
    end

    name = String(read_element(f, swap_bytes, UInt8))

    local data
    if class == mxCELL_CLASS
        data = read_cell(f, swap_bytes, dimensions, subsys)
    elseif class == mxSTRUCT_CLASS || class == mxOBJECT_CLASS
        data = read_struct(f, swap_bytes, dimensions, class == mxOBJECT_CLASS, subsys)
    elseif class == mxSPARSE_CLASS
        data = read_sparse(f, swap_bytes, dimensions, flags)
    elseif class == mxCHAR_CLASS && length(dimensions) <= 2
        data = read_string(f, swap_bytes, dimensions)
    elseif class == mxFUNCTION_CLASS
        data = read_matrix(f, swap_bytes, subsys)
    elseif class == mxOPAQUE_CLASS
        data = read_opaque(f, swap_bytes, subsys)
    else
        if (flags[1] & (1 << 9)) != 0 # logical
            data = read_data(f, swap_bytes, Bool, dimensions)
        else
            convert_type = CONVERT_TYPES[class]
            data = read_data(f, swap_bytes, convert_type, dimensions)
            if (flags[1] & (1 << 11)) != 0 # complex
                data = complex_array(data, read_data(f, swap_bytes, convert_type, dimensions))
            end
        end
    end

    return (name, data)
end

# Open MAT file for reading
function matopen(ios::IOStream, endian_indicator::UInt16)
    matfile = Matlabv5File(ios, endian_indicator == 0x494D)

    seek(matfile.ios, 116)
    subsys_offset = read_bswap(matfile.ios, matfile.swap_bytes, UInt64)
    if subsys_offset == 0x2020202020202020
        subsys_offset = UInt64(0)
    end
    if subsys_offset != 0
        matfile.subsystem_position = subsys_offset
        read_subsystem!(matfile)
    end

    return matfile
end

# Read whole MAT file
function read(matfile::Matlabv5File)
    vars = Dict{String, Any}()

    seek(matfile.ios, 128)
    while !eof(matfile.ios)
        pos = position(matfile.ios)
        if pos == matfile.subsystem_position
            # Skip reading subsystem again
            (_, nbytes) = read_header(matfile.ios, matfile.swap_bytes)
            skip(matfile.ios, nbytes)
            continue
        end
        (name, data) = read_matrix(matfile.ios, matfile.swap_bytes, matfile.subsystem)
        vars[name] = data
    end
    vars
end

function read_subsystem!(matfile::Matlabv5File)
    ios = matfile.ios
    swap_bytes = matfile.swap_bytes
    seek(ios, matfile.subsystem_position)
    (_, subsystem_data) = read_matrix(ios, swap_bytes, matfile.subsystem)
    buf = IOBuffer(vec(subsystem_data))
    seek(buf, 8) # Skip subsystem header
    _, subsys_data = read_matrix(buf, swap_bytes, matfile.subsystem)
    MAT_subsys.load_subsys!(matfile.subsystem, subsys_data, swap_bytes)
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
                f = ZlibDecompressorStream(matfile.ios)
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

Base.haskey(matfile::Matlabv5File, varname::String) =
    haskey(getvarnames(matfile), varname)
Base.keys(matfile::Matlabv5File) =
    keys(getvarnames(matfile))

# Read a variable from a MAT file
function read(matfile::Matlabv5File, varname::String)
    varnames = getvarnames(matfile)
    if !haskey(varnames, varname)
        error("no variable $varname in file")
    end
    seek(matfile.ios, varnames[varname])
    (name, data) = read_matrix(matfile.ios, matfile.swap_bytes, matfile.subsystem)
    data
end

# Complain about writing to a MAT file
function write(parent::Matlabv5File, name::String, s)
    error("Writing to a MATLAB v5 file is not currently supported. Create a new file instead.")
end

# Close MAT file
close(matfile::Matlabv5File) = close(matfile.ios)
end
