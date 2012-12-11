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

load("hdf5")
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

    function MatlabHDF5File(id, filename, toclose::Bool, writeheader::Bool)
        f = new(id, filename, toclose, writeheader)
        if toclose
            finalizer(f, close)
        end
        f
    end
end
MatlabHDF5File(id, filename, toclose) = MatlabHDF5File(id, filename, toclose, false)
MatlabHDF5File(id, filename) = MatlabHDF5File(id, filename, true, false)
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
    MatlabHDF5File(f, filename, true, writeheader)
end

### Matlab file format specification ###

const name_type_attr_matlab = "MATLAB_class"
const empty_attr_matlab = "MATLAB_empty"

function read(dset::HDF5Dataset{MatlabHDF5File})
    if exists(dset, empty_attr_matlab)
        # Empty arrays encode the dimensions as the dataset
        dims = int(read(plain(dset)))
        mattype = a_read(dset, name_type_attr_matlab)
        return Array(str2type_matlab[mattype], dims...)
    end
    # Read the MATLAB class
    mattype = "cell"
    local T
    if exists(dset, name_type_attr_matlab)
        mattype = a_read(dset, name_type_attr_matlab)
        # Convert to Julia type
        T = str2type_matlab[mattype]
    end
    # Check for a COMPOUND data set, and if so handle complex numbers specially
    dtype = datatype(dset)
    class_id = HDF5.h5t_get_class(dtype.id)
    if class_id == HDF5.H5T_COMPOUND
        if !check_datatype_complex(dtype)
            close(dtype)
            error("Unrecognized compound data type when reading ", name(dset))
        end
        close(dtype)
        T = abstr_eltype(T)
        if !(T == Float32 || T == Float64)
            error("For complex numbers, only Float32 and Float64 supported")
            # only Complex64 and Complex128 are defined as BitsKinds
        end
        memtype = build_datatype_complex(T)
        sz = size(dset)
        buf = Array(Uint8, 2*sizeof(T), sz...)
        HDF5.h5d_read(dset.id, memtype.id, HDF5.H5S_ALL, HDF5.H5S_ALL, HDF5.H5P_DEFAULT, buf)
        C = (T == Float32) ? Complex64 : Complex128
        d = reinterpret(C, buf, sz)
        return numel(d) == 1 ? d[1] : d
    end
    close(dtype)
    # Read the dataset
    if mattype == "cell"
        # Represented as an array of refs
        refs = read(plain(dset), Array{HDF5ReferenceObj})
        out = Array(Any, size(refs))
        f = file(dset)
        for i = 1:numel(refs)
            out[i] = read(f[refs[i]])
        end
        return out
    end
    d = read(plain(dset), T)
    numel(d) == 1 ? d[1] : d
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

for (fsym, dsym) in
    ((:(write{T<:HDF5BitsKind}), :T),
     (:(write{T<:HDF5BitsKind}), :(Array{T})))
    @eval begin
        function ($fsym)(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, data::$dsym)
            local typename
            # Determine the Matlab type
            if has(type2str_matlab, T)
                typename = type2str_matlab[T]
            else
                error("Type ", T, " is not (yet) supported")
            end
            # Everything in Matlab is an array
            empty = false
            if !isa(data, Array)
                data = [data]
            elseif isempty(data)
                empty = true
                data = [size(data)...]
            end
            # Create the dataset
            dset, dtype = d_create(plain(parent), name, data)
            try
                if empty
                    a_write(dset, empty_attr_matlab, uint8(1))
                end
                # Write the attribute
                a_write(dset, name_type_attr_matlab, typename)
                # Write the data
                writearray(dset, dtype.id, data)
            finally
                close(dset)
                close(dtype)
            end
        end
    end
end

# Write a string
function write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, str::String)
    # Here we assume no UTF-16
    data = zeros(Uint16, strlen(str))
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
function write{T}(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, data::Array{T}, l::Int)
    pathrefs = "/#refs#"
    local g
    local refs
    if !exists(parent, pathrefs)
        g = g_create(file(parent), pathrefs)
    else
        g = parent[pathrefs]
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
        if l == -1
            l = length(g)-1
        end
        for i = 1:length(data)
            itemname = string(l+i)
            if isa(data[i], Array)
                # Need to pass level so that we don't create items with the same names
                write(g, itemname, data[i], l+length(data))
            else
                write(g, itemname, data[i])
            end
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
write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, data::Array) = write(parent, name, data, -1)

# Check that keys are valid for a struct, and convert them to an array of ASCIIStrings
function check_struct_keys(k::Vector)
    asckeys = Array(ASCIIString, length(k))
    for i = 1:length(k)
        key = k[i]
        if !isa(key, String)
            error("Only Dicts with string keys may be saved as MATLAB structs")
        elseif match(r"^[a-zA-Z][a-zA-Z0-9_]*$", key) == nothing 
            error("Invalid key \"$key\" for MATLAB struct: keys must start with a letter and contain only alphanumeric characters and underscore")
        elseif length(key) > 63
            error("Invalid key \"$key\" for MATLAB struct: keys must be less than 64 characters")
        end
        asckeys[i] = convert(ASCIIString, key)
    end
    asckeys
end

# Write a struct from arrays of keys and values
function write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, k::Vector{ASCIIString}, v::Vector)
    g = g_create(parent, name)
    gplain = plain(g)
    a_write(gplain, name_type_attr_matlab, "struct")
    for i = 1:length(k)
        write(g, k[i], v[i])
    end
    a_write(gplain, "MATLAB_fields", HDF5Vlen(k))
end

# Write Associative as a struct
write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, s::Associative) =
    write(parent, name, check_struct_keys(keys(s)), values(s))

# Write generic CompositeKind as a struct
function write(parent::Union(MatlabHDF5File, HDF5Group{MatlabHDF5File}), name::ByteString, s)
    T = typeof(s)
    if !isa(T, CompositeKind)
        error("This is the write function for CompositeKind, but the input doesn't fit")
    end
    write(parent, name, check_struct_keys([string(x) for x in T.names]), [getfield(s, x) for x in T.names])
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
            datap[i] = rstrip(CharString(squeeze(data[i, :])))
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

end
