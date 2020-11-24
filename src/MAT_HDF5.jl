# MAT_HDF5.jl
# Tools for reading MATLAB HDF5 (v7.3) files in Julia
#
# Copyright (C) 2012   Timothy E. Holy
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

###########################################
## Reading and writing MATLAB .mat files ##
###########################################

module MAT_HDF5

using HDF5, SparseArrays
# deprecated for HDF5 v0.14+, but use deprecated binding to have common function with
# e.g. JLD.jl
import HDF5: exists

import Base: names, read, write, close
import HDF5: Reference

const HDF5Parent = Union{HDF5.File, HDF5.Group}
const HDF5BitsOrBool = Union{HDF5.BitsType,Bool}

mutable struct MatlabHDF5File <: HDF5.H5DataStore
    plain::HDF5.File
    toclose::Bool
    writeheader::Bool
    refcounter::Int
    compress::Bool

    function MatlabHDF5File(plain, toclose::Bool=true, writeheader::Bool=false, refcounter::Int=0, compress::Bool=false)
        f = new(plain, toclose, writeheader, refcounter, compress)
        if toclose
            finalizer(close, f)
        end
        f
    end
end

"""
    close(matfile_handle)

Close a Matlab file.
"""
function close(f::MatlabHDF5File)
    if f.toclose
        close(f.plain)
        if f.writeheader
            magic = zeros(UInt8, 512)
            identifier = "MATLAB 7.3 MAT-file" # minimal but sufficient
            GC.@preserve magic identifier begin
                magicptr = pointer(magic)
                idptr = pointer(identifier)
                unsafe_copyto!(magicptr, idptr, length(identifier))
            end
            magic[126] = 0x02
            magic[127] = 0x49
            magic[128] = 0x4d
            rawfid = open(f.plain.filename, "r+")
            write(rawfid, magic)
            close(rawfid)
        end
        f.toclose = false
    end
    nothing
end

function matopen(filename::AbstractString, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool, compress::Bool)
    local f
    if ff && !wr
        error("Cannot append to a write-only file")
    end
    if !cr && !isfile(filename)
        error("File ", filename, " cannot be found")
    end
    pa = create_property(HDF5.H5P_FILE_ACCESS; fclose_degree = HDF5.H5F_CLOSE_STRONG)
    if cr && (tr || !isfile(filename))
        # We're truncating, so we don't have to check the format of an existing file
        # Set the user block to 512 bytes, to save room for the header
        p = create_property(HDF5.H5P_FILE_CREATE; userblock = 512)
        f = HDF5.h5f_create(filename, HDF5.H5F_ACC_TRUNC, p.id, pa.id)
        writeheader = true
    else
        f = HDF5.h5f_open(filename, wr ? HDF5.H5F_ACC_RDWR : HDF5.H5F_ACC_RDONLY, pa.id)
        writeheader = false
    end
    close(pa)
    fid = MatlabHDF5File(HDF5.File(f, filename), true, writeheader, 0, compress)
    pathrefs = "/#refs#"
    if haskey(fid.plain, pathrefs)
        g = fid.plain[pathrefs]
        fid.refcounter = length(g)-1
        close(g)
    end
    fid
end

### Matlab file format specification ###

const name_type_attr_matlab = "MATLAB_class"
const empty_attr_matlab = "MATLAB_empty"
const sparse_attr_matlab = "MATLAB_sparse"
const int_decode_attr_matlab = "MATLAB_int_decode"

### Reading
function read_complex(dtype::HDF5.Datatype, dset::HDF5.Dataset, ::Type{T}) where T
    if !check_datatype_complex(dtype)
        close(dtype)
        error("Unrecognized compound data type when reading ", HDF5.name(dset))
    end
    return read(dset, Complex{T})
end

function m_read(dset::HDF5.Dataset)
    if haskey(dset, empty_attr_matlab)
        # Empty arrays encode the dimensions as the dataset
        dims = convert(Vector{Int}, read(dset))
        mattype = read_attribute(dset, name_type_attr_matlab)
        if mattype == "char"
            return ""
        elseif mattype == "struct"
            # Not sure if this check is necessary but it is checked in
            # `m_read(g::HDF5.Group)`
            if haskey(dset, "MATLAB_fields")
                return Dict{String,Any}(join(n)=>[] for n in read_attribute(dset, "MATLAB_fields"))
            else
                return Dict{String,Any}()
            end
        else
            T = mattype == "canonical empty" ? Union{} : str2type_matlab[mattype]
            return Array{T}(undef, dims...)
        end
    end

    mattype = haskey(dset, name_type_attr_matlab) ? read_attribute(dset, name_type_attr_matlab) : "cell"

    if mattype == "cell"
        # Cell arrays, represented as an array of refs
        refs = read(dset, Reference)
        out = Array{Any}(undef, size(refs))
        f = HDF5.file(dset)
        for i = 1:length(refs)
            dset = f[refs[i]]
            try
                out[i] = m_read(dset)
            finally
                close(dset)
            end
        end
        return out
    elseif !haskey(str2type_matlab,mattype)
        @warn "MATLAB $mattype values are currently not supported"
        return missing
    end

    # Regular arrays of values
    # Convert to Julia type
    T = str2type_matlab[mattype]

    # Check for a COMPOUND data set, and if so handle complex numbers specially
    dtype = datatype(dset)
    try
        class_id = HDF5.h5t_get_class(dtype.id)
        d = class_id == HDF5.H5T_COMPOUND ? read_complex(dtype, dset, T) : read(dset, T)
        length(d) == 1 ? d[1] : d
    finally
        close(dtype)
    end
end

function add!(A, x)
    for i = 1:length(A)
        @inbounds A[i] += x
    end
    A
end

# reading a struct, struct array, or sparse matrix
function m_read(g::HDF5.Group)
    mattype = read_attribute(g, name_type_attr_matlab)
    if mattype != "struct"
        # Check if this is a sparse matrix.
        fn = keys(g)
        if haskey(attributes(g), sparse_attr_matlab)
            # This is a sparse matrix.
            # ir is the row indices, jc is the column boundaries.
            # We add one to account for the zero-based (MATLAB) to one-based (Julia) transition
            jc = add!(convert(Vector{Int}, read(g, "jc")), 1)
            if "data" in fn && "ir" in fn && "jc" in fn
                # This matrix is not empty.
                ir = add!(convert(Vector{Int}, read(g, "ir")), 1)
                dset = g["data"]
                T = str2type_matlab[mattype]
                try
                    dtype = datatype(dset)
                    class_id = HDF5.h5t_get_class(dtype.id)
                    try
                        data = class_id == HDF5.H5T_COMPOUND ? read_complex(dtype, dset, T) : read(dset, T)
                    finally
                        close(dtype)
                    end
                finally
                    close(dset)
                end
            else
                # This matrix is empty.
                ir = Int[]
                data = str2type_matlab[mattype][]
            end
            return SparseMatrixCSC(convert(Int, read_attribute(g, sparse_attr_matlab)), length(jc)-1, jc, ir, data)
        elseif mattype == "function_handle"
            @warn "MATLAB $mattype values are currently not supported"
            return missing
        else
            error("Cannot read from a non-struct group, type was $mattype")
        end
    end
    if haskey(g, "MATLAB_fields")
        fn = [join(f) for f in read_attribute(g, "MATLAB_fields")]
    else
        fn = keys(g)
    end
    s = Dict{String, Any}()
    for i = 1:length(fn)
        dset = g[fn[i]]
        try
            s[fn[i]] = m_read(dset)
        finally
            close(dset)
        end
    end
    s
end

"""
    read(matfile_handle, varname) -> value

Read a variable from an opened Matlab file and return its value.

See `matopen` and `matread`.
"""
function read(f::MatlabHDF5File, name::String)
    local val
    obj = f.plain[name]
    try
        val = m_read(obj)
    finally
        close(obj)
    end
    val
end

"""
    names(matfile_handle) -> Vector{String}

Return a list of variables in an opened Matlab file.

See `matopen`.
"""
names(f::MatlabHDF5File) = keys(f)

"""
    exists(matfile_handle, varname) -> Bool

Return true if a variable is present in an opened Matlab file.

See `matopen`.
"""
exists(p::MatlabHDF5File, path::String) = haskey(p, path)

# HDF5v0.14+ H5DataStore uses keys/haskey
Base.keys(f::MatlabHDF5File) = filter!(x -> x!="#refs#" && x!="#subsystem#", keys(f.plain))
Base.haskey(p::MatlabHDF5File, path::String) = haskey(p.plain, path)

### Writing

# Check whether a varname is valid for MATLAB
check_valid_varname(s::AbstractString) = if match(r"^[a-zA-Z][a-zA-Z0-9_]*$", s) == nothing
    error("Invalid variable name or key \"$s\": variable names must start with a letter and contain only alphanumeric characters and underscore")
elseif length(s) > 63
    error("Invalid variable name or key \"$s\": variable names must be less than 64 characters")
end

toarray(x::Array) = x
toarray(x::Array{Bool}) = reinterpret(UInt8, x)
toarray(x::Bool) = UInt8[x]
toarray(x) = [x]

# Write the MATLAB type string for dset
m_writetypeattr(dset, ::Type{Complex{T}}) where T = m_writetypeattr(dset, T)
function m_writetypeattr(dset, T)
    if !haskey(type2str_matlab, T)
        error("Type ", T, " is not (yet) supported")
    end
    typename = type2str_matlab[T]

    # Write the attribute
    write_attribute(dset, name_type_attr_matlab, typename)
    if T == Bool
        write_attribute(dset, int_decode_attr_matlab, Int32(1))
    end
end

# Writes an empty scalar or array
function m_writeempty(parent::HDF5Parent, name::String, data::AbstractArray)
    adata = [size(data)...]
    dset, dtype = create_dataset(parent, name, adata)
    try
        write_attribute(dset, empty_attr_matlab, 0x01)
        m_writetypeattr(dset, eltype(data))
        write_dataset(dset, dtype, adata)
    finally
        close(dset)
        close(dtype)
    end
end

# Write an array to a dataset in a MATLAB file, returning the dataset
function m_writearray(parent::HDF5Parent, name::String, adata::AbstractArray{T}, compress::Bool) where {T<:HDF5BitsOrBool}
    if compress
        dset, dtype = create_dataset(parent, name, adata;
                               compress = 3, chunk = HDF5.heuristic_chunk(adata))
    else
        dset, dtype = create_dataset(parent, name, adata)
    end
    try
        write_dataset(dset, dtype, adata)
        dset
    catch e
        close(dset)
        rethrow(e)
    finally
        close(dtype)
    end
end
function m_writearray(parent::HDF5Parent, name::String, adata::AbstractArray{Complex{T}}, compress::Bool) where {T<:HDF5BitsOrBool}
    dtype = build_datatype_complex(T)
    try
        stype = dataspace(adata)
        if compress
            dset = create_dataset(parent, name, dtype, stype;
                              compress = 3, chunk = HDF5.heuristic_chunk(adata))
        else
            dset = create_dataset(parent, name, dtype, stype)
        end
        try
            arr = reshape(reinterpret(T, adata), tuple(2, size(adata)...))
            write_dataset(dset, dtype, arr)
        catch e
            close(dset)
            rethrow(e)
        finally
            close(stype)
        end
        return dset
    finally
        close(dtype)
    end
end

# Write a scalar or array
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, data::Union{T, Complex{T}, Array{T}, Array{Complex{T}}}) where {T<:HDF5BitsOrBool}
    if isempty(data)
        m_writeempty(parent, name, data)
        return
    end
    dset = m_writearray(parent, name, toarray(data), mfile.compress)
    try
        m_writetypeattr(dset, T)
    finally
        close(dset)
    end
end

# Write sparse arrays
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, data::SparseMatrixCSC{T}) where T
    g = create_group(parent, name)
    try
        m_writetypeattr(g, T)
        write_attribute(g, sparse_attr_matlab, UInt64(size(data, 1)))
        if !isempty(data.nzval)
            close(m_writearray(g, "data", toarray(data.nzval), mfile.compress))
            close(m_writearray(g, "ir", add!(isa(data.rowval, Vector{UInt64}) ? copy(data.rowval) : convert(Vector{UInt64}, data.rowval), typemax(UInt64)), mfile.compress))
        end
        close(m_writearray(g, "jc", add!(isa(data.colptr, Vector{UInt64}) ? copy(data.colptr) : convert(Vector{UInt64}, data.colptr), typemax(UInt64)), mfile.compress))
    finally
        close(g)
    end
end

# Write BitArray as Array{Bool}. Would be better not to require the conversion, but this is easy
m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, s::BitArray) =
    m_write(mfile, parent, name, convert(Array{Bool}, s))

# Write a string
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, str::AbstractString)
    if isempty(str)
        data = UInt64[0, 0]

        # Create the dataset
        dset, dtype = create_dataset(parent, name, data)
        try
            write_attribute(dset, name_type_attr_matlab, "char")
            write_attribute(dset, empty_attr_matlab, 0x01)
            write_dataset(dset, dtype, data)
        finally
            close(dset)
            close(dtype)
        end
    else
        # Here we assume no UTF-16
        data = zeros(UInt16, 1, length(str))
        i = 1
        for c in str
            data[i] = c
            i += 1
        end

        # Create the dataset
        dset, dtype = create_dataset(parent, name, data)
        try
            write_attribute(dset, name_type_attr_matlab, "char")
            write_attribute(dset, int_decode_attr_matlab, Int32(2))
            write_dataset(dset, dtype, data)
        finally
            close(dset)
            close(dtype)
        end
    end
end

# Write cell arrays
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, data::Array{T}) where T
    pathrefs = "/#refs#"
    fid = HDF5.file(parent)
    local g
    local refs
    if !haskey(fid, pathrefs)
        g = create_group(fid, pathrefs)
    else
        g = fid[pathrefs]
    end
    try
        # If needed, create the "empty" item
        if !haskey(g, "a")
            edata = zeros(UInt64, 2)
            eset, etype = create_dataset(g, "a", edata)
            try
                write_dataset(eset, etype, edata)
                write_attribute(eset, name_type_attr_matlab, "canonical empty")
                write_attribute(eset, "MATLAB_empty", 0x00)
            finally
                close(etype)
                close(eset)
            end
        else
            a = g["a"]
            if !haskey(attributes(a), "MATLAB_empty")
                error("Must create the empty item, with name a, first")
            end
            close(a)
        end
        # Write the items to the reference group
        refs = Array{Reference}(undef, size(data))
        for i = 1:length(data)
            mfile.refcounter += 1
            itemname = string(mfile.refcounter)
            m_write(mfile, g, itemname, data[i])
            # Extract references
            tmp = g[itemname]
            refs[i] = Reference(tmp, pathrefs*"/"*itemname)
            close(tmp)
        end
    finally
        close(g)
    end
    # Write the references as the chosen variable
    cset, ctype = create_dataset(parent, name, refs)
    try
        write_dataset(cset, ctype, refs)
        write_attribute(cset, name_type_attr_matlab, "cell")
    finally
        close(ctype)
        close(cset)
    end
end

# Check that keys are valid for a struct, and convert them to an array of ASCIIStrings
function check_struct_keys(k::Vector)
    asckeys = Vector{String}(undef, length(k))
    for i = 1:length(k)
        key = k[i]
        if !isa(key, AbstractString)
            error("Only Dicts with string keys may be saved as MATLAB structs")
        end
        check_valid_varname(key)
        asckeys[i] = convert(String, key)
    end
    asckeys
end

# Write a struct from arrays of keys and values
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, k::Vector{String}, v::Vector)
    g = create_group(parent, name)
    write_attribute(g, name_type_attr_matlab, "struct")
    for i = 1:length(k)
        m_write(mfile, g, k[i], v[i])
    end
    write_attribute(g, "MATLAB_fields", HDF5.VLen(k))
end

# Write Associative as a struct
m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, s::AbstractDict) =
    m_write(mfile, parent, name, check_struct_keys(collect(keys(s))), collect(values(s)))

# Write generic CompositeKind as a struct
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, s)
    if isbits(s)
        error("This is the write function for CompositeKind, but the input doesn't fit")
    end
    T = typeof(s)
    m_write(mfile, parent, name, check_struct_keys([string(x) for x in fieldnames(T)]), [getfield(s, x) for x in fieldnames(T)])
end

# Check whether a variable name is valid, then write it
"""
    write(matfile_handle, varname, value)

Write the value into an opened Matlab file as the specified variable.

See `matopen` and `matwrite`.
"""
function write(parent::MatlabHDF5File, name::String, thing)
    check_valid_varname(name)
    m_write(parent, parent.plain, name, thing)
end

## Type conversion operations ##

struct MatlabString end

const str2type_matlab = Dict(
    "canonical empty" => nothing,
    "int8"    => Int8,
    "uint8"   => UInt8,
    "int16"   => Int16,
    "uint16"  => UInt16,
    "int32"   => Int32,
    "uint32"  => UInt32,
    "int64"   => Int64,
    "uint64"  => UInt64,
    "single"  => Float32,
    "double"  => Float64,
    "cell"    => Any,
    "char"    => MatlabString,
    "logical" => Bool
)
const type2str_matlab = Dict(
    Int8    => "int8",
    UInt8   => "uint8",
    Int16   => "int16",
    UInt16  => "uint16",
    Int32   => "int32",
    UInt32  => "uint32",
    Int64   => "int64",
    UInt64  => "uint64",
    Float32 => "single",
    Float64 => "double",
    Bool    => "logical"
)


function read(obj::Union{HDF5.Dataset,HDF5.Attribute}, ::Type{MatlabString})
    T = HDF5.get_jl_type(obj)
    data = read(obj, T)
    if size(data, 1) == 1
        sz = size(data)
        data = reshape(data, sz[2:end])
    end
    if ndims(data) == 1
        return String(convert(Vector{Char}, data))
    elseif ndims(data) == 2
        return datap = String[rstrip(String(convert(Vector{Char}, vec(data[i, :])))) for i = 1:size(data, 1)]
    else
        return data
    end
end

## Utilities for handling complex numbers
function build_datatype_complex(T::Type)
    memtype = create_datatype(HDF5.H5T_COMPOUND, 2*sizeof(T))
    HDF5.h5t_insert(memtype, "real", 0, HDF5.hdf5_type_id(T))
    HDF5.h5t_insert(memtype, "imag", sizeof(T), HDF5.hdf5_type_id(T))
    return memtype
end

function check_datatype_complex(dtype::HDF5.Datatype)
    n = HDF5.h5t_get_nmembers(dtype.id)
    if n != 2
        return false
    end
    if HDF5.h5t_get_member_name(dtype.id, 0) != "real" ||
       HDF5.h5t_get_member_name(dtype.id, 1) != "imag"
        return false
    end
    true
end

end
