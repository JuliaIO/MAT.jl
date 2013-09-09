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
using HDF5
# Add methods to...
import HDF5: close, names, read, write

# Debugging: comment this block out if you un-modulize hdf5.jl
# Types
Hid = HDF5.Hid
HDF5ReferenceObj = HDF5.HDF5ReferenceObj
HDF5BitsKind = HDF5.HDF5BitsKind
# Constants
H5F_ACC_RDONLY = HDF5.H5F_ACC_RDONLY
H5F_ACC_RDWR = HDF5.H5F_ACC_RDWR
H5F_ACC_TRUNC = HDF5.H5F_ACC_TRUNC
H5F_CLOSE_STRONG = HDF5.H5F_CLOSE_STRONG
H5P_DEFAULT = HDF5.H5P_DEFAULT
H5P_FILE_ACCESS = HDF5.H5P_FILE_ACCESS
H5P_FILE_CREATE = HDF5.H5P_FILE_CREATE
# Functions
h5f_close  = HDF5.h5f_close
h5f_create = HDF5.h5f_create
h5f_open   = HDF5.h5f_open
writearray = HDF5.writearray
hdf5_to_julia = HDF5.hdf5_to_julia

type MatlabHDF5File
    plain::HDF5File
    toclose::Bool
    writeheader::Bool
    refcounter::Int

    function MatlabHDF5File(plain, toclose::Bool=true, writeheader::Bool=false, refcounter::Int=0)
        f = new(plain, toclose, writeheader, refcounter)
        if toclose
            finalizer(f, close)
        end
        f
    end
end

function close(f::MatlabHDF5File)
    if f.toclose
        close(f.plain)
        if f.writeheader
            magic = zeros(Uint8, 512)
            const identifier = "MATLAB 7.3 MAT-file" # minimal but sufficient
            magic[1:length(identifier)] = identifier.data
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

function matopen(filename::String, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool)
    local f
    if ff && !wr
        error("Cannot append to a write-only file")
    end
    if !cr && !isfile(filename)
        error("File ", filename, " cannot be found")
    end
    pa = p_create(H5P_FILE_ACCESS)
    pa["fclose_degree"] = H5F_CLOSE_STRONG
    if cr && (tr || !isfile(filename))
        # We're truncating, so we don't have to check the format of an existing file
        # Set the user block to 512 bytes, to save room for the header
        p = p_create(H5P_FILE_CREATE)
        p["userblock"] = 512
        f = h5f_create(filename, H5F_ACC_TRUNC, p.id, pa.id)
        writeheader = true
    else
        f = h5f_open(filename, wr ? H5F_ACC_RDWR : H5F_ACC_RDONLY, pa.id)
        writeheader = false
    end
    close(pa)
    fid = MatlabHDF5File(HDF5File(f, filename), true, writeheader, 0)
    pathrefs = "/#refs#"
    if exists(fid.plain, pathrefs)
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

### Reading
function read_complex{T}(dtype::HDF5Datatype, dset::HDF5Dataset, ::Type{Array{T}})
    if !check_datatype_complex(dtype)
        close(dtype)
        error("Unrecognized compound data type when reading ", name(dset))
    end
    memtype = build_datatype_complex(T)
    sz = size(dset)
    st = sizeof(T)
    buf = Array(Uint8, 2*st, sz...)
    HDF5.h5d_read(dset.id, memtype.id, HDF5.H5S_ALL, HDF5.H5S_ALL, HDF5.H5P_DEFAULT, buf)

    if T == Float32
        d = reinterpret(Complex64, buf, sz)
    elseif T == Float64
        d = reinterpret(Complex128, buf, sz)
    else
        d = reinterpret(T, slicedim(buf, 1, 1:st), sz) + im*reinterpret(T, slicedim(buf, 1, st+1:2*st), sz)
    end
    length(d) == 1 ? d[1] : d
end

function m_read(dset::HDF5Dataset)
    if exists(dset, empty_attr_matlab)
        # Empty arrays encode the dimensions as the dataset
        dims = int(read(dset))
        mattype = a_read(dset, name_type_attr_matlab)
        T = mattype == "canonical empty" ? None : str2type_matlab[mattype]
        return Array(T, dims...)
    end

    mattype = exists(dset, name_type_attr_matlab) ? a_read(dset, name_type_attr_matlab) : "cell"

    if mattype == "cell"
        # Cell arrays, represented as an array of refs
        refs = read(dset, Array{HDF5ReferenceObj})
        out = Array(Any, size(refs))
        f = file(dset)
        for i = 1:length(refs)
            dset = f[refs[i]]
            try
                out[i] = m_read(dset)
            finally
                close(dset)
            end
        end
        return out
    end

    # Regular arrays of values
    # Convert to Julia type
    T = str2type_matlab[mattype]

    # Check for a COMPOUND data set, and if so handle complex numbers specially
    dtype = datatype(dset)
    class_id = HDF5.h5t_get_class(dtype.id)
    d = class_id == HDF5.H5T_COMPOUND ? read_complex(dtype, dset, T) : read(dset, T)
    
    close(dtype)
    length(d) == 1 ? d[1] : d
end

# reading a struct or struct array
function m_read(g::HDF5Group)
    mattype = a_read(g, name_type_attr_matlab)
    if mattype != "struct"
        # Check if this is a sparse matrix.
        fn = names(g)
        if exists(attrs(g), sparse_attr_matlab)
            
            # This is a sparse matrix.
            # ir is the row indices, jc is the column boundaries.
            # We add one to account for the zero-based (MATLAB) to one-based (Julia) transition
            jc = read(g, "jc") + 1
            if fn == ["data", "ir", "jc"]
                # This matrix is not empty.
                ir = read(g, "ir") + 1
                data = read(g, "data")
            else
                # This matrix is empty.
                ir = HDF5.hdf5_to_julia_eltype(HDF5.datatype(g["jc"]))[]
                data = str2eltype_matlab[mattype][]
            end
            return SparseMatrixCSC(HDF5.a_read(g, "MATLAB_sparse"), length(jc) -1, jc, ir, data)
        else
            error("Cannot read from a non-struct group, type was $mattype")
        end
    end
    if exists(g, "MATLAB_fields")
        fn = a_read(g, "MATLAB_fields")
    else
        fn = names(g)
    end
    s = Dict{ASCIIString, Any}()
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

function read(f::MatlabHDF5File, name::ASCIIString)
    local val
    obj = f.plain[name]
    try
        val = m_read(obj)
    finally
        close(obj)
    end
    val
end

names(f::MatlabHDF5File) = names(f.plain)

# read every variable in the file
function read(f::MatlabHDF5File)
    vars = names(f)
    vars = vars[!(vars .== "#refs#")]  # delete "#refs#"
    vals = Array(Any, length(vars))
    for i = 1:length(vars)
        vals[i] = read(f, vars[i])
    end
    Dict(vars, vals)
end

### Writing

# Check whether a varname is valid for MATLAB
check_valid_varname(s::String) = if match(r"^[a-zA-Z][a-zA-Z0-9_]*$", s) == nothing 
    error("Invalid variable name or key \"$s\": variable names must start with a letter and contain only alphanumeric characters and underscore")
elseif length(s) > 63
    error("Invalid variable name or key \"$s\": variable names must be less than 64 characters")
end

toarray(x::Array) = x
toarray(x) = [x]

# Writes an array given a dataset and datatype, and then closes them
function m_writearray{T <: HDF5BitsKind}(dset::HDF5Dataset, dtype::HDF5Datatype, name::ByteString, data::Array{T})
    try
        # Determine the Matlab type
        if !haskey(type2str_matlab, T)
            error("Type ", T, " is not (yet) supported")
        end
        typename = type2str_matlab[T]

        # Write the attribute
        a_write(dset, name_type_attr_matlab, typename)
        # Write the data
        writearray(dset, dtype.id, data)
    finally
        close(dset)
        close(dtype)
    end
end

# Writes an empty scalar or array
function m_writeempty(parent::Union(HDF5File, HDF5Group), name::ByteString, data::Array)
    adata = [size(data)...]
    dset, dtype = d_create(parent, name, adata)
    a_write(dset, empty_attr_matlab, uint8(1))
    m_writearray(dset, dtype, name, adata)
end

# Write a scalar or array
function m_write{T<:HDF5BitsKind}(mfile::MatlabHDF5File, parent::Union(HDF5File, HDF5Group), name::ByteString, data::Union(T, Array{T}))
    if isempty(data)
        m_writeempty(parent, name, data)
        return
    end
    adata = toarray(data)

    dset, dtype = d_create(parent, name, adata)
    m_writearray(dset, dtype, name, adata)
end

# Write a complex scalar or array
function m_write{T}(mfile::MatlabHDF5File, parent::Union(HDF5File, HDF5Group), name::ByteString, data::Union(Complex{T}, Array{Complex{T}}))
    if isempty(data)
        m_writeempty(parent, name, data)
        return
    end
    adata = toarray(data)

    dtype = build_datatype_complex(T)
    stype = dataspace(adata)
    obj_id = HDF5.h5d_create(parent.id, name, dtype.id, stype.id)
    dset = HDF5Dataset(obj_id, file(parent))
    arr = isbits(T) ? reinterpret(T, adata, tuple(2, size(adata)...)) : [real(adata), imag(adata)]
    m_writearray(dset, dtype, name, arr)
end

# Write a string
function m_write(mfile::MatlabHDF5File, parent::Union(HDF5File, HDF5Group), name::ByteString, str::String)
    # Here we assume no UTF-16
    data = zeros(Uint16, length(str))
    i = 1
    for c in str
        data[i] = c
        i += 1
    end

    # Create the dataset
    dset, dtype = d_create(parent, name, data)
    try
        # Write the attribute
        a_write(dset, name_type_attr_matlab, "char")
        # Write the data
        writearray(dset, dtype.id, data)
    finally
        close(dset)
        close(dtype)
    end
end

# Write cell arrays
function m_write{T}(mfile::MatlabHDF5File, parent::Union(HDF5File, HDF5Group), name::ByteString, data::Array{T})
    pathrefs = "/#refs#"
    fid = file(parent)
    local g
    local refs
    if !exists(fid, pathrefs)
        g = g_create(fid, pathrefs)
    else
        g = fid[pathrefs]
    end
    try
        # If needed, create the "empty" item
        if !exists(g, "a")
            edata = zeros(Uint64, 2)
            eset, etype = d_create(g, "a", edata)
            try
                writearray(eset, etype.id, edata)
                a_write(eset, name_type_attr_matlab, "canonical empty")
                a_write(eset, "MATLAB_empty", uint8(0))
            finally
                close(etype)
                close(eset)
            end
        else
            a = g["a"]
            if !exists(attrs(a), "MATLAB_empty")
                error("Must create the empty item, with name a, first")
            end
            close(a)
        end
        # Write the items to the reference group
        refs = Array(HDF5ReferenceObj, size(data)...)
        for i = 1:length(data)
            mfile.refcounter += 1
            itemname = string(mfile.refcounter)
            m_write(mfile, g, itemname, data[i])
            # Extract references
            tmp = g[itemname]
            refs[i] = HDF5ReferenceObj(tmp, pathrefs*"/"*itemname)
            close(tmp)
        end
    finally
        close(g)
    end
    # Write the references as the chosen variable
    cset, ctype = d_create(parent, name, refs)
    try
        writearray(cset, ctype.id, refs)
        a_write(cset, name_type_attr_matlab, "cell")
    finally
        close(ctype)
        close(cset)
    end
end

# Check that keys are valid for a struct, and convert them to an array of ASCIIStrings
function check_struct_keys(k::Vector)
    asckeys = Array(ASCIIString, length(k))
    for i = 1:length(k)
        key = k[i]
        if !isa(key, String)
            error("Only Dicts with string keys may be saved as MATLAB structs")
        end
        check_valid_varname(key)
        asckeys[i] = convert(ASCIIString, key)
    end
    asckeys
end

# Write a struct from arrays of keys and values
function m_write(mfile::MatlabHDF5File, parent::Union(HDF5File, HDF5Group), name::ByteString, k::Vector{ASCIIString}, v::Vector)
    g = g_create(parent, name)
    a_write(g, name_type_attr_matlab, "struct")
    for i = 1:length(k)
        m_write(mfile, g, k[i], v[i])
    end
    a_write(g, "MATLAB_fields", HDF5Vlen(k))
end

# Write Associative as a struct
m_write(mfile::MatlabHDF5File, parent::Union(HDF5File, HDF5Group), name::ByteString, s::Associative) =
    m_write(mfile, parent, name, check_struct_keys(collect(keys(s))), collect(values(s)))

# Write generic CompositeKind as a struct
function m_write(mfile::MatlabHDF5File, parent::Union(HDF5File, HDF5Group), name::ByteString, s)
    if isbits(s)
        error("This is the write function for CompositeKind, but the input doesn't fit")
    end
    T = typeof(s)
    m_write(mfile, parent, name, check_struct_keys([string(x) for x in T.names]), [getfield(s, x) for x in T.names])
end

# Check whether a variable name is valid, then write it
function write(parent::MatlabHDF5File, name::ByteString, thing)
    check_valid_varname(name)
    m_write(parent, parent.plain, name, thing)
end

## Type conversion operations ##

type MatlabString; end

const str2type_matlab = {
    "canonical empty" => nothing,
    "int8"    => Array{Int8},
    "uint8"   => Array{Uint8},
    "int16"   => Array{Int16},
    "uint16"  => Array{Uint16},
    "int32"   => Array{Int32},
    "uint32"  => Array{Uint32},
    "int64"   => Array{Int64},
    "uint64"  => Array{Uint64},
    "single"  => Array{Float32},
    "double"  => Array{Float64},
    "cell"    => Array{Any},
    "char"    => MatlabString,
    "logical" => Bool,
}
# These operate on the element type rather than the whole type
const str2eltype_matlab = {
    "canonical empty" => nothing,
    "int8"    => Int8,
    "uint8"   => Uint8,
    "int16"   => Int16,
    "uint16"  => Uint16,
    "int32"   => Int32,
    "uint32"  => Uint32,
    "int64"   => Int64,
    "uint64"  => Uint64,
    "single"  => Float32,
    "double"  => Float64,
    "cell"    => Any,
    "char"    => MatlabString,
    "logical" => Bool,
}
const type2str_matlab = {
    Int8    => "int8",
    Uint8   => "uint8",
    Int16   => "int16",
    Uint16  => "uint16",
    Int32   => "int32",
    Uint32  => "uint32",
    Int64   => "int64",
    Uint64  => "uint64",
    Float32 => "single",
    Float64 => "double"
}


function read(obj::HDF5Object, ::Type{MatlabString})
    T = hdf5_to_julia(obj)
    data = read(obj, T)
    if size(data, 1) == 1
        sz = size(data)
        data = reshape(data, sz[2:end])
    end
    if ndims(data) == 1
        return bytestring(CharString(data))
    elseif ndims(data) == 2
        datap = Array(String, size(data, 1))
        for i = 1:length(datap)
            datap[i] = rstrip(bytestring(CharString(reshape(data[i, :], size(data,2)))))
        end
        return datap
    else
        return data
    end
end
function read(obj::HDF5Object, ::Type{Bool})
    tf = read(obj, Uint8)
    tf > 0
end

## Utilities for handling complex numbers
function build_datatype_complex(T::Type)
    memtype_id = HDF5.h5t_create(HDF5.H5T_COMPOUND, 2*sizeof(T))
    HDF5.h5t_insert(memtype_id, "real", 0, HDF5.hdf5_type_id(T))
    HDF5.h5t_insert(memtype_id, "imag", sizeof(T), HDF5.hdf5_type_id(T))
    HDF5Datatype(memtype_id)
end

function check_datatype_complex(dtype::HDF5Datatype)
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
