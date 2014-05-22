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
using Zlib, HDF5
import Base: read, write, close
import HDF5: names, exists

type Matlabv5File{I<:IO} <: HDF5.DataFile
    ios::I
    swap_bytes::Bool
    subsys_offset::FileOffset
    varnames::Dict{ASCIIString, FileOffset}
    subsystem::Dict{ASCIIString, Any}

    Matlabv5File(ios::I, swap_bytes::Bool, subsys_offset) = new(ios, swap_bytes, subsys_offset)
end
Matlabv5File{I<:IO}(ios::I, swap_bytes::Bool, subsys_offset) = Matlabv5File{I}(ios, swap_bytes, subsys_offset)
# TODO: Remove if/when github.com/JuliaLang/julia/issues/5333 lands for mutables
function Matlabv5File(matfile::Matlabv5File, f::IO)
    v = Matlabv5File(f, matfile.swap_bytes, 0)
    if isdefined(matfile, :varnames)
        v.varnames = matfile.varnames
    end
    if isdefined(matfile, :subsystem)
        v.subsystem = matfile.subsystem
    end
    v
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
const mxFUNCTION_CLASS = 16 # undocumented (function handles)
const mxOPAQUE_CLASS = 17   # undocumented (classdef objects)
const READ_TYPES = Type[Int8, Uint8, Int16, Uint16, Int32, Uint32, Float32, None, Float64,
    None, None, Int64, Uint64, Uint8, Uint8]
const CONVERT_TYPES = Type[None, None, None, None, None, Float64, Float32, Int8, Uint8,
    Int16, Uint16, Int32, Uint32, Int64, Uint64, None, None]

const SUBSYS_HEADER_PADDING = 120

read_bswap{T}(f::IO, swap_bytes::Bool, ::Type{T}) =
    swap_bytes ? bswap(read(f, T)) : read(f, T)
read_bswap{T}(f::IO, swap_bytes::Bool, ::Type{T}, dim::Union(Int, (Int...))) =
    swap_bytes ? [bswap(x) for x in read(f, T, dim)] : read(f, T, dim)

skip_padding(f::IO, nbytes::Int64, hbytes::Int) = if nbytes % hbytes != 0
    skip(f, hbytes-(nbytes % hbytes))
end

# Read data type and number of bytes at the start of a data element
function read_header(f::IO, swap_bytes::Bool)
    dtype = read_bswap(f, swap_bytes, Uint32)

    if (dtype & 0xFFFF0000) != 0
        # Small Data Element Format
        (dtype & 0x0000FFFF, int64(dtype >> 16), 4)
    else
        # Data Element Format
        (dtype, int64(read_bswap(f, swap_bytes, Uint32)), 8)
    end
end

# Read data element as a vector of a given type
function read_element{T}(f::IO, swap_bytes::Bool, ::Type{T})
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    data = read_bswap(f, swap_bytes, T, int(div(nbytes, sizeof(T))))
    skip_padding(f, nbytes, hbytes)
    data
end

# Read data element as encoded type
function read_data(f::IO, swap_bytes::Bool)
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    read_type = READ_TYPES[dtype]
    data = read_bswap(f, swap_bytes, read_type, int(div(nbytes, sizeof(read_type))))
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
        data = read_bswap(f, swap_bytes, read_type, tuple(int(dimensions)...))
    else
        data = read_bswap(f, swap_bytes, read_type)
    end
    skip_padding(f, nbytes, hbytes)

    read_array ? convert(Array{T}, data) : convert(T, data)
end

function read_cell(matfile::Matlabv5File, dimensions::Vector{Int32})
    data = cell(int(dimensions)...)
    for i = 1:length(data)
        (ignored_name, data[i]) = read_matrix(matfile)
    end
    data
end

function read_struct(matfile::Matlabv5File, dimensions::Vector{Int32}, is_object::Bool)
    f, swap_bytes = matfile.ios, matfile.swap_bytes
    field_length = read_element(f, swap_bytes, Int32)[1]
    field_names = read_element(f, swap_bytes, Uint8)
    n_fields = div(length(field_names), field_length)
    if is_object
        class = ascii(read_element(f, swap_bytes, Uint8))
    end


    # Get field names as strings
    field_name_strings = Array(String, n_fields)
    n_el = prod(dimensions)
    for i = 1:n_fields
        sname = field_names[(i-1)*field_length+1:i*field_length]
        index = findfirst(sname, 0)
        field_name_strings[i] = ascii(index == 0 ? sname : sname[1:index-1])
    end

    data = Dict{ASCIIString, Any}()
    sizehint(data, n_fields+1)
    if is_object
        data["class"] = class
    end

    if n_el == 1
        # Read a single struct into a dict
        for field_name in field_name_strings
            data[field_name] = read_matrix(matfile)[2]
        end
    else
        # Read multiple structs into a dict of arrays
        for field_name in field_name_strings
            data[field_name] = cell(dimensions...)
        end
        for i = 1:n_el
            for field_name in field_name_strings
                data[field_name][i] = read_matrix(matfile)[2]
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

function read_sparse(matfile::Matlabv5File, dimensions::Vector{Int32}, flags::Vector{Uint32})
    f, swap_bytes = matfile.ios, matfile.swap_bytes
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
    ir = plusone!(int(read_element(f, swap_bytes, Int32)))
    jc = plusone!(int(read_element(f, swap_bytes, Int32)))
    if (flags[1] & (1 << 9)) != 0 # logical
        # WTF. For some reason logical sparse matrices are tagged as doubles.
        pr = read_element(f, swap_bytes, Bool)
    else
        pr = read_data(f, swap_bytes)
        if (flags[1] & (1 << 11)) != 0 # complex
            pr = complex(pr, read_data(f, swap_bytes))
        end
    end

    SparseMatrixCSC(m, n, jc, ir, pr)
end

function read_string(matfile::Matlabv5File, dimensions::Vector{Int32})
    f, swap_bytes = matfile.ios, matfile.swap_bytes
    (dtype, nbytes, hbytes) = read_header(f, swap_bytes)
    if dtype <= 2 || dtype == 16
        # If dtype <= 2, this may give an error on non-ASCII characters, since the string
        # would be ISO-8859-1 and not UTF-8. However, MATLAB 2012b always saves strings with
        # a 2-byte encoding in v6 format, and saves UTF-8 in v7 format. Thus, this may never
        # happen in the wild.
        chars = read(f, Uint8, nbytes)
        if dimensions[1] <= 1
            data = bytestring(chars)
        else
            data = Array(ByteString, dimensions[1])
            for i = 1:dimensions[1]
                data[i] = rstrip(bytestring(chars[i:dimensions[1]:end]))
            end
        end
    elseif dtype <= 4 || dtype == 17
        # Technically, if dtype == 3 or dtype == 4, this is ISO-8859-1 and not Unicode.
        # However, the first 256 Unicode code points are derived from ISO-8859-1, so UCS-2
        # is a superset of 2-byte ISO-8859-1.
        chars = read_bswap(f, swap_bytes, Uint16, int(div(nbytes, 2)))
        for i = 1:length(chars)
            if chars[i] > 255
                # Newer versions of MATLAB seem to write some mongrel UTF-8...
                chars[i] = bytestring(reverse(reinterpret(Uint8, [chars[i]])))[1]
            end
        end
        if dimensions[1] <= 1
            data = bytestring(CharString(chars))
        else
            data = Array(ByteString, dimensions[1])
            for i = 1:dimensions[1]
                data[i] = rstrip(bytestring(CharString(chars[i:dimensions[1]:end])))
            end
        end
    else
        error("Unsupported string type")
    end
    skip_padding(f, nbytes, hbytes)
    data
end

function read_opaque(matfile::Matlabv5File)
    f, swap_bytes = matfile.ios, matfile.swap_bytes
    # This is undocumented. There are two strings followed by an array
    # The first string is always (?) MCOS: "Matlab Common/Class Object System"
    class_system = ascii(read_element(f, swap_bytes, Uint8))
    class_system == "MCOS" || error("unknown opaque class system: ", class_system)
    # The second string is the object's class name
    class_name = ascii(read_element(f, swap_bytes, Uint8))

    # The subsystem contains a FileWrapper__; an opaque object unlike all others
    class_name == "FileWrapper__" && return read_matrix(matfile)[2]

    # This is an unnamed array of uint32s: indexes into the subsystem
    ids = read_matrix(matfile)[2]
    ids[1] == 0xdd00_0000 || error("unknown opaque sentinal value: ", idxs[1])
    ndims = ids[2]
    dims = ids[3:2+ndims]
    n_el = isempty(dims) ? 0 : prod(dims)
    object_ids = ids[3+ndims:2+ndims+n_el]
    class_id = ids[3+ndims+n_el]
    @assert 3+ndims+n_el == length(ids) "unexpected data at end of opaque id"
    
    # Return like a struct or struct array
    obj = Dict{ASCIIString,Any}()
    obj["_class"] = matfile.subsystem["_classes"][class_id]
    if n_el == 1
        obj = merge!(obj,matfile.subsystem["_objects"][object_ids[1]])
    elseif n_el > 1
        # Multiple objects as a dict of arrays, use the first as a template?
        fields = keys(matfile.subsystem["_objects"][object_ids[1]])
        for fld in fields
            obj[fld] = cell(dims...)
        end
        for i = 1:n_el
            sys_obj = matfile.subsystem["_objects"][object_ids[i]]
            for fld in fields
                obj[fld][i] = sys_obj[fld]
            end
        end
    end
    
    return obj
end

# Read matrix data
function read_matrix(matfile::Matlabv5File)
    f, swap_bytes = matfile.ios, matfile.swap_bytes
    (dtype, nbytes) = read_header(f, swap_bytes)
    if dtype == miCOMPRESSED
        bytes = decompress(read(f, Uint8, nbytes))
        mi = Matlabv5File(matfile,IOBuffer(bytes))
        output = read_matrix(mi)
        close(mi)
        return output
    elseif dtype != miMATRIX
        error("Unexpected data type")
    elseif nbytes == 0
        # If one creates a cell array using
        #     y = cell(matfile, n)
        # then MATLAB will save the empty cells as zero-byte matrices. If one creates a
        # empty cells using
        #     a = {[], [], []}
        # then MATLAB does not save the empty cells as zero-byte matrices. To avoid
        # surprises, we produce an empty array in both cases.
        return ("", Array(None, 0, 0))
    end

    flags = read_element(f, swap_bytes, Uint32)
    class = flags[1] & 0xFF
    # Opaque objects are dimensionless
    dimensions = (class == mxOPAQUE_CLASS) ? Int32[] : read_element(f, swap_bytes, Int32)
    name = ascii(read_element(f, swap_bytes, Uint8))

    local data
    if class == mxCELL_CLASS
        data = read_cell(matfile, dimensions)
    elseif class == mxSTRUCT_CLASS || class == mxOBJECT_CLASS
        data = read_struct(matfile, dimensions, class == mxOBJECT_CLASS)
    elseif class == mxSPARSE_CLASS
        data = read_sparse(matfile, dimensions, flags)
    elseif class == mxCHAR_CLASS && length(dimensions) <= 2
        data = read_string(matfile, dimensions)
    elseif class == mxFUNCTION_CLASS
        # Unnamed submatrix with a struct of function implementation details
        data = read_matrix(matfile)[2]
    elseif class == mxOPAQUE_CLASS
        data = read_opaque(matfile)
    else
        convert_type = CONVERT_TYPES[class]
        data = read_data(f, swap_bytes, convert_type, dimensions)
        if (flags[1] & (1 << 11)) != 0 # complex
            data = complex(data, read_data(f, swap_bytes, convert_type, dimensions))
        elseif (flags[1] & (1 << 9)) != 0 # logical
            data = reinterpret(Bool, uint8(data))
        end
    end

    return (name, data)
end

# Open MAT file for reading
function matopen(filename::String, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool)
    if wr || cr || tr || ff
        error("Creating or appending to MATLAB v5 files is not supported")
    end

    ios = open(filename, "r")
    header = read(ios, Uint8, 116)
    subsys_offset = read(ios, Uint64)
    if subsys_offset == 0x2020_2020_2020_2020
        subsys_offset = zero(subsys_offset)
    end
    version = read(ios, Uint16)
    endian_indicator = read(ios, Uint16)

    local swap_bytes
    if endian_indicator == 0x4D49
        swap_bytes = false
    elseif endian_indicator == 0x494D
        swap_bytes = true
    else
        error("Invalid endian indicator")
    end

    if swap_bytes
        version = bswap(version)
    end
    if version != 0x0100
        error("Unsupported MATLAB file version")
    end

    matfile = Matlabv5File(ios, swap_bytes, subsys_offset)

    if subsys_offset > 0
        seek(ios, subsys_offset)
        matfile.subsystem = read_subsystem(matfile)
    end

    matfile
end

# Undocumented methods for reading the subsystem
is_subsystem_matfile(x::Any) = false;
function is_subsystem_matfile{N}(x::Array{Uint8,N})
    N > 2          && return false
    length(x) < 12 && return false

    # Check first 4 characters in either endianness
    (x[1:4] == [0x00,0x01,0x49,0x4d] || x[1:4] == [0x01,0x00,0x4d,0x49])
end

function read_subsystem(matfile::Matlabv5File)
    # Read the final, unnamed matrix hidden at the end of the file
    name, data = read_matrix(matfile)
    @assert isempty(name) && is_subsystem_matfile(data) "invalid subsystem"
    @assert eof(matfile.ios) "unread data at end of subsystem"

    sys = read_subsystem_matfile(data);

    # Parse the FileWrapper__ object - might there be more than one?
    mcos = sys["_i1"]["MCOS"]
    sys["_classes"], sys["_objects"] = parse_filewrapper(mcos)

    sys
end

function read_subsystem_matfile{N}(data::Array{Uint8,N})
    # A Matfile is stored within this matrix's data
    f = IOBuffer(vec(data))
    # Check endianness and versioning.
    seek(f,2)
    endian_flag = read(f,Uint16)
    swap_bytes = endian_flag == 0x494D
    @assert swap_bytes || endian_flag == 0x4D49 "unknown endian flags in subsystem"
    seek(f,0)
    @assert read_bswap(f,swap_bytes,Uint16) == 0x0100 "unsupported MATLAB file version in subsystem"
    seek(f,8)

    matfile = Matlabv5File(f,swap_bytes,0);
    svars = Dict{ASCIIString,Any}()
    i=0;
    while !eof(matfile.ios)
        name,data = read_matrix(matfile)

        if isempty(name) && is_subsystem_matfile(data)
            # There are sometimes nested subsystems? Are they ever non-empty?
            data = read_subsystem_matfile(data)
        end
        name = isempty(name) ? "_i$(i+=1)" : name;
        svars[name] = data
    end
    close(f)

    svars
end

# Returns a tuple of parsed (classes, objects)
function parse_filewrapper(mcos::Array{Any})
    # Parse the first element as an IOBuffer data stream
    f = IOBuffer(vec(mcos[1]))
    offsets,strs = parse_filewrapper_header(f)
    classes = parse_filewrapper_classes(f,strs,offsets)
    objects = parse_filewrapper_objects(f,strs,mcos,offsets)

    (classes, objects)
end

function parse_filewrapper_header(f)
    id = read(f,Uint32) # First element is a version number? Always 2?
    id == 2 || error("unknown first field (version/id?): ", id)

    n_strs = read(f,Uint32) # Second element is the number of strings
    offsets = read(f,Uint32,6) # Followed by up to 6 section offsets

    # And two unknown/reserved fields
    all(read(f,Uint32,2) .== 0) || error("reserved header fields nonzero")

    # The string data section: a stream of null-delimited strings
    @assert position(f) == 0x28
    strs = Array(ASCIIString,n_strs)
    for i = 1:n_strs
        strs[i] = readuntil(f, '\0')[1:end-1] # drop the trailing null byte
    end

    (offsets,strs)
end

function parse_filewrapper_classes(f,strs,offsets)
    # Class information is a set of four Int32s per class, but the first four
    # values always seem to be 0. Are they reserved? Or simply never used?
    seek(f,offsets[1])
    all(read(f,Int32,4) .== 0) || error("unknown header to class information")

    n = div(offsets[2]-offsets[1],4*4) - 1
    classes = Array((ASCIIString,ASCIIString),n)
    for i=1:n
        package_idx = read(f,Int32)
        package = package_idx > 0 ? strs[package_idx] : ""
        name_idx = read(f,Int32)
        name = name_idx > 0 ? strs[name_idx] : ""
        all(read(f,Uint32,2) .== 0) || error("discovered a nonzero class property for ",name)
        classes[i] = (package, name)
    end
    @assert position(f) == offsets[2]

    classes
end

function parse_filewrapper_objects(f,names,heap,offsets)
    seek(f,offsets[2])
    seg2_props = parse_filewrapper_props(f,names,heap,offsets[3])

    seek(f,offsets[3])
    # Again, first set of elements are all zero
    all(read(f,Int32,6) .== 0) || error("unknown header to object information")

    n = div(offsets[4]-position(f),6*4)
    # 6 values per obj: class_idx, unknown, unknown, seg2_idx, seg4_idx, obj_id
    object_info = Array(Int32,6,n)
    read(f,object_info)
    object_info = object_info.' # Transpose to column-major
    @assert all(object_info[:,2:3] .== 0) "discovered a nonzero unknown object property"

    @assert position(f) == offsets[4]
    seg4_props = parse_filewrapper_props(f,names,heap,offsets[5])

    objects = Array(Dict{ASCIIString,Any},n)
    for i=1:n
        class_idx = object_info[i,1]
        seg2_idx, seg4_idx = object_info[i,4:5]
        obj_id = object_info[i,6]

        if seg2_idx > 0 && seg4_idx == 0
            props = seg2_props[seg2_idx]
        elseif seg4_idx > 0 && seg2_idx == 0
            props = seg4_props[seg4_idx]
        else
            error("unable to find property for object with id ", obj_id)
        end

        # Merge it with the matfile defaults for this class
        @assert obj_id <= n "obj_ids are assumed to be continuous"
        objects[obj_id] = merge(heap[end][class_idx+1],props)
    end
    objects
end

function parse_filewrapper_props(f,names,heap,section_end)
    props = Array(Dict{ASCIIString,Any},0)
    position(f) >= section_end && return props

    # Another first element that's always zero; reserved? or ignored?
    all(read(f,Int32,2) .== 0) || error("unknown header to properties segment")

    # We have to guess on the array size; 8 int32s would be 2 props per object
    sizehint(props,iceil((section_end-position(f))/(8*4)))

    while position(f) < section_end
        # For each object, there is first an Int32 for the number of properties
        nprops = read(f,Int32)
        d = Dict{ASCIIString,Any}()
        sizehint(d,nprops)
        for i=1:nprops
            # And then three Int32s for each property
            name  = names[read(f,Int32)]
            flag  = read(f,Int32) # A flag describing how the value is used
            value = read(f,Int32)

            if flag == 0  # The prop is stored in the names array
                d[name] = names[value]
            elseif flag == 1 # The prop is stored in the MCOS FileWrapper__ heap
                d[name] = heap[value+3] # But… it's off by 3!? Crazy.
            elseif flag == 2 # The prop is a boolean, and it's the value itself
                @assert 0 <= value <= 1 "boolean flag has a value other than 0 or 1"
                d[name] = bool(value)
            else
                error("unknown flag ",flag," for property ",name," with value ",value)
            end
        end
        push!(props,d)

        # Jump to the next 8-byte aligned offset
        if position(f) % 8 != 0
            seek(f,iceil(position(f)/8)*8)
        end
    end
    props
end

# Read whole MAT file
function read(matfile::Matlabv5File)
    seek(matfile.ios, 128)
    vars = Dict{ASCIIString, Any}()
    i = 0;
    while !eof(matfile.ios) && (matfile.subsys_offset==0 || position(matfile.ios) < matfile.subsys_offset)
        (name, data) = read_matrix(matfile)
        name = isempty(name) ? "_i$(i+=1)" : name;
        vars[name] = data
    end

    vars
end

# Read only variable names from an HDF5 file
function getvarnames(matfile::Matlabv5File)
    if !isdefined(matfile, :varnames)
        seek(matfile.ios, 128)
        matfile.varnames = varnames = Dict{ASCIIString, FileOffset}()
        while !eof(matfile.ios) && (matfile.subsys_offset==0 || position(matfile.ios) < matfile.subsys_offset)
            offset = position(matfile.ios)
            (dtype, nbytes, hbytes) = read_header(matfile.ios, matfile.swap_bytes)
            if dtype == miCOMPRESSED
                # Read 1K or less of data
                read_bytes = min(nbytes, 1024)
                source = read(matfile.ios, Uint8, read_bytes)
                skip(matfile.ios, nbytes-read_bytes)

                # Uncompress to 1K byte buffer
                dest = zeros(Uint8, 1024)
                dest_buf_size = Int[1024]

                ret = ccall((:uncompress, @windows ? "zlib1" : "libz"), Int32, (Ptr{Uint8}, Ptr{Int}, Ptr{Uint8}, Uint),
                    dest, dest_buf_size, source, length(source))

                # Zlib may complain because the buffer is small or the data are incomplete
                if ret != Zlib.Z_OK && ret != Zlib.Z_BUF_ERROR && ret != Zlib.Z_DATA_ERROR
                    throw(ZError(ret))
                end

                # Create IOBuffer from uncompressed buffer
                f = IOBuffer(dest)

                # Read header
                read_header(f, matfile.swap_bytes)
            elseif dtype == miMATRIX
                f = matfile.ios
            else
                error("Unexpected data type")
            end

            class = read_element(f, matfile.swap_bytes, Uint32)[1] & 0xFF
            (class != mxOPAQUE_CLASS) && read_element(f, matfile.swap_bytes, Int32)
            varnames[ascii(read_element(f, matfile.swap_bytes, Uint8))] = offset

            if dtype == miCOMPRESSED
                # Close IOBuffer
                close(f)
            else
                # Seek past
                seek(f, offset+nbytes+hbytes)
            end
        end
    end
    matfile.varnames
end

exists(matfile::Matlabv5File, varname::ASCIIString) =
    haskey(getvarnames(matfile), varname)
names(matfile::Matlabv5File) =
    keys(getvarnames(matfile))

# Read a variable from a MAT file
function read(matfile::Matlabv5File, varname::ASCIIString)
    varnames = getvarnames(matfile)
    if !haskey(varnames, varname)
        error("no variable $varname in file")
    end
    seek(matfile.ios, varnames[varname])
    (name, data) = read_matrix(matfile)
    data
end

# Complain about writing to a MAT file
function write(parent::Matlabv5File, name::ByteString, s)
    error("Writing to a MATLAB v5 file is not currently supported. Create a new file instead.")
end

# Close MAT file
close(matfile::Matlabv5File) = close(matfile.ios)
end
