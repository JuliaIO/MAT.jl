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
import HDF5.close, HDF5.read, HDF5.write

# Debugging: comment this block out if you un-modulize hdf5.jl
# Types
Hid = HDF5.Hid
HDF5ReferenceObj = HDF5.HDF5ReferenceObj
HDF5ReferenceObjArray = HDF5.HDF5ReferenceObjArray
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

type MatlabHDF5File <: HDF5File
    id::Hid
    filename::String
    toclose::Bool
    writeheader::Bool
    refcounter::Int

    function MatlabHDF5File(id, filename, toclose::Bool, writeheader::Bool, refcounter::Int)
        f = new(id, filename, toclose, writeheader, refcounter)
        if toclose
            finalizer(f, close)
        end
        f
    end
end
MatlabHDF5File(id, filename, toclose) = MatlabHDF5File(id, filename, toclose, false, 0)
MatlabHDF5File(id, filename) = MatlabHDF5File(id, filename, true, false, 0)
function close(f::MatlabHDF5File)
    if f.toclose
        h5f_close(f.id)
        if f.writeheader
            magic = zeros(Uint8, 512)
            const identifier = "MATLAB 7.3 MAT-file" # minimal but sufficient
            magic[1:length(identifier)] = identifier.data
            magic[126] = 0x02
            magic[127] = 0x49
            magic[128] = 0x4d
            rawfid = open(f.filename, "r+")
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
    fid = MatlabHDF5File(f, filename, true, writeheader, 0)
    pathrefs = "/#refs#"
    if has(fid, pathrefs)
        g = fid[pathrefs]
        fid.refcounter = length(g)-1
        close(g)
    end
    fid
end

### Matlab file format specification ###

const name_type_attr_matlab = "MATLAB_class"
const empty_attr_matlab = "MATLAB_empty"

### Reading
function read_complex(dtype::HDF5Datatype, dset::HDF5Dataset{MatlabHDF5File}, T::Type)
    if !check_datatype_complex(dtype)
        close(dtype)
        error("Unrecognized compound data type when reading ", name(dset))
    end
    T = abstr_eltype(T)
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

function read(dset::HDF5Dataset{MatlabHDF5File})
    if exists(dset, empty_attr_matlab)
        # Empty arrays encode the dimensions as the dataset
        dims = int(read(plain(dset)))
        mattype = a_read(dset, name_type_attr_matlab)
        T = mattype == "canonical empty" ? None : str2type_matlab[mattype]
        return Array(T, dims...)
    end

    mattype = exists(dset, name_type_attr_matlab) ? a_read(dset, name_type_attr_matlab) : "cell"

    if mattype == "cell"
        # Cell arrays, represented as an array of refs
        refs = read(plain(dset), Array{HDF5ReferenceObj})
        out = Array(Any, size(refs))
        f = file(dset)
        for i = 1:length(refs)
            out[i] = read(f[refs[i]])
        end
        return out
    end

    # Regular arrays of values
    # Convert to Julia type
    T = str2type_matlab[mattype]

    # Check for a COMPOUND data set, and if so handle complex numbers specially
    dtype = datatype(dset)
    class_id = HDF5.h5t_get_class(dtype.id)
    d = class_id == HDF5.H5T_COMPOUND ? read_complex(dtype, dset, T) : read(plain(dset), T)
    
    close(dtype)
    length(d) == 1 ? d[1] : d
end

# reading a struct or struct array
function read(g::HDF5Group{MatlabHDF5File})
    mattype = a_read(g, name_type_attr_matlab)
    if mattype != "struct"
        error("Cannot read from a non-struct group")
    end
    if exists(g, "MATLAB_fields")
        fn = a_read(g, "MATLAB_fields")
    else
        fn = names(g)
    end
    s = Dict{ASCIIString, Any}()
    for i = 1:length(fn)
        s[fn[i]] = read(g, fn[i])
    end
    s
end

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
function m_writeempty(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, data::Array)
    adata = [size(data)...]
    dset, dtype = d_create(plain(parent), name, adata)
    a_write(dset, empty_attr_matlab, uint8(1))
    m_writearray(dset, dtype, name, adata)
end

# Write a scalar or array
function m_write{T<:HDF5BitsKind}(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, data::Union(T, Array{T}))
    if isempty(data)
        m_writeempty(parent, name, data)
        return
    end
    adata = toarray(data)

    dset, dtype = d_create(plain(parent), name, adata)
    m_writearray(dset, dtype, name, adata)
end

# Write a complex scalar or array
function m_write{C<:Complex}(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, data::Union(C, Array{C}))
    if isempty(data)
        m_writeempty(parent, name, data)
        return
    end
    adata = toarray(data)

    T = realtype(C)
    dtype = build_datatype_complex(T)
    stype = dataspace(adata)
    obj_id = HDF5.h5d_create(parent.id, name, dtype.id, stype.id)
    dset = HDF5Dataset(obj_id, plain(file(parent)))
    arr = isbits(C) ? reinterpret(T, adata, tuple(2, size(adata)...)) : [real(adata), imag(adata)]
    m_writearray(dset, dtype, name, arr)
end

# Write a string
function m_write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, str::String)
    # Here we assume no UTF-16
    data = zeros(Uint16, length(str))
    i = 1
    for c in str
        data[i] = c
        i += 1
    end

    # Create the dataset
    dset, dtype = d_create(plain(parent), name, data)
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
function m_write{T}(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, data::Array{T})
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
            pg = plain(g)
            edata = zeros(Uint64, 2)
            eset, etype = d_create(pg, "a", edata)
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
        refs = HDF5ReferenceObjArray(size(data)...)
        for i = 1:length(data)
            fid.refcounter += 1
            itemname = string(fid.refcounter)
            m_write(g, itemname, data[i])
            # Extract references
            tmp = g[itemname]
            refs[i] = (tmp, pathrefs*"/"*itemname)
            close(tmp)
        end
    finally
        close(g)
    end
    # Write the references as the chosen variable
    cset, ctype = d_create(plain(parent), name, refs)
    try
        writearray(cset, ctype.id, refs.r)
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
function m_write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, k::Vector{ASCIIString}, v::Vector)
    g = g_create(parent, name)
    gplain = plain(g)
    a_write(gplain, name_type_attr_matlab, "struct")
    for i = 1:length(k)
        m_write(g, k[i], v[i])
    end
    a_write(gplain, "MATLAB_fields", HDF5Vlen(k))
end

# Write Associative as a struct
m_write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, s::Associative) =
    m_write(parent, name, check_struct_keys(collect(keys(s))), collect(values(s)))

# Write generic CompositeKind as a struct
function m_write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, s)
    if isbits(s)
        error("This is the write function for CompositeKind, but the input doesn't fit")
    end
    T = typeof(s)
    m_write(parent, name, check_struct_keys([string(x) for x in T.names]), [getfield(s, x) for x in T.names])
end

# Check whether a variable name is valid, then write it
function write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, thing)
    check_valid_varname(name)
    m_write(parent, name, thing)
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


function read(obj::HDF5Object{PlainHDF5File}, ::Type{MatlabString})
    T = hdf5_to_julia(obj)
    data = read(obj, T)
    if size(data, 1) == 1
        sz = size(data)
        data = reshape(data, sz[2:end])
    end
    if ndims(data) == 1
        return CharString(data)
    elseif ndims(data) == 2
        datap = Array(String, size(data, 1))
        for i = 1:length(datap)
            datap[i] = rstrip(CharString(reshape(data[i, :], size(data,2))))
        end
        return datap
    else
        return data
    end
end
function read(obj::HDF5Object{PlainHDF5File}, ::Type{Bool})
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

abstr_eltype{T}(::Type{Array{T}}) = T

realtype(::Type{Complex64}) = Float32
realtype(::Type{Complex128}) = Float64
realtype{T}(::Type{Complex{T}}) = T
realtype{T}(::Type{ComplexPair{T}}) = T

end
