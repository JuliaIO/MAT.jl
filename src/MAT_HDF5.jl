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
using ..MAT_subsys

import Base: names, read, write, close
import HDF5: Reference
import Dates
import Tables
import PooledArrays: PooledArray

import ..MAT_types:
    convert_struct_array,
    EmptyStruct,
    MatlabClassObject,
    MatlabOpaque,
    MatlabStructArray,
    MatlabTable,
    ScalarOrArray,
    StructArrayField

const HDF5Parent = Union{HDF5.File, HDF5.Group}
const HDF5BitsOrBool = Union{HDF5.BitsType,Bool}

mutable struct MatlabHDF5File <: HDF5.H5DataStore
    plain::HDF5.File
    toclose::Bool
    writeheader::Bool
    refcounter::Int
    compress::Bool
    subsystem::Subsystem

    function MatlabHDF5File(plain, toclose::Bool=true, writeheader::Bool=false, refcounter::Int=0, compress::Bool=false)
        f = new(plain, toclose, writeheader, refcounter, compress, Subsystem())
        if toclose
            finalizer(close, f)
        end
        f
    end
end

function Base.show(io::IO, f::MatlabHDF5File)
    print(io, "MatlabHDF5File(")
    print(io, f.plain, ", ")
    print(io, f.toclose, ", ")
    print(io, f.writeheader, ", ")
    print(io, f.refcounter, ", ")
    print(io, f.compress, ")")
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
            if Base.ENDIAN_BOM == 0x04030201
                magic[127] = 0x49
                magic[128] = 0x4d
            else
                magic[127] = 0x4d
                magic[128] = 0x49
            end
            rawfid = open(f.plain.filename, "r+")
            write(rawfid, magic)
            close(rawfid)
        end
        f.toclose = false
    end
    nothing
end

function matopen(filename::AbstractString, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool, compress::Bool, endian_indicator::Bool; table::Type=MatlabTable, convert_opaque::Bool=true)
    local f
    if ff && !wr
        error("Cannot append to a read-only file")
    end
    if !cr && !isfile(filename)
        error("File ", filename, " cannot be found")
    end
    fapl = HDF5.FileAccessProperties()
    fapl.fclose_degree = :strong
    if cr && (tr || !isfile(filename))
        # We're truncating, so we don't have to check the format of an existing file
        # Set the user block to 512 bytes, to save room for the header
        fcpl = HDF5.FileCreateProperties()
        fcpl.userblock = 512
        f = HDF5.API.h5f_create(filename, HDF5.API.H5F_ACC_TRUNC, fcpl, fapl)
        writeheader = true
    else
        f = HDF5.API.h5f_open(filename, wr ? HDF5.API.H5F_ACC_RDWR : HDF5.API.H5F_ACC_RDONLY, fapl)
        writeheader = false
    end
    close(fapl)
    fid = MatlabHDF5File(HDF5.File(f, filename), true, writeheader, 0, compress)
    pathrefs = "/#refs#"
    if haskey(fid.plain, pathrefs)
        g = fid.plain[pathrefs]
        fid.refcounter = length(g)-1
        close(g)
    end
    subsys_refs = "#subsystem#"
    if rd && haskey(fid.plain, subsys_refs)
        fid.subsystem.table_type = table
        fid.subsystem.convert_opaque = convert_opaque
        subsys_data = m_read(fid.plain[subsys_refs], fid.subsystem)
        MAT_subsys.load_subsys!(fid.subsystem, subsys_data, endian_indicator)
    elseif wr
        MAT_subsys.init_save!(fid.subsystem)
    end
    fid
end

### Matlab file format specification ###

const name_type_attr_matlab = "MATLAB_class"
const empty_attr_matlab = "MATLAB_empty"
const sparse_attr_matlab = "MATLAB_sparse"
const int_decode_attr_matlab = "MATLAB_int_decode"
const object_type_attr_matlab = "MATLAB_object_decode"
const object_decode_attr_matlab = "MATLAB_object_decode"
const struct_field_attr_matlab = "MATLAB_fields"

### Reading
function read_complex(dtype::HDF5.Datatype, dset::HDF5.Dataset, ::Type{T}) where T
    if !check_datatype_complex(dtype)
        close(dtype)
        error("Unrecognized compound data type when reading ", HDF5.name(dset))
    end
    return read(dset, Complex{T})
end

function read_cell(dset::HDF5.Dataset, subsys::Subsystem)
    refs = read(dset, Reference)
    out = Array{Any}(undef, size(refs))
    f = HDF5.file(dset)
    for i = 1:length(refs)
        dset = f[refs[i]]
        try
            out[i] = m_read(dset, subsys)
        finally
            close(dset)
        end
    end
    return out
end

function m_read(dset::HDF5.Dataset, subsys::Subsystem)
    if haskey(dset, empty_attr_matlab)
        # Empty arrays encode the dimensions as the dataset
        dims = convert(Vector{Int}, read(dset))
        mattype = read_attribute(dset, name_type_attr_matlab)
        if mattype == "char"
            return ""
        elseif mattype == "struct"
            # Not sure if this check is necessary but it is checked in
            # `m_read(g::HDF5.Group)`
            if haskey(dset, struct_field_attr_matlab)
                field_names = [join(n) for n in read_attribute(dset, struct_field_attr_matlab)]
                return MatlabStructArray(field_names, tuple(dims...))
            else
                return Dict{String,Any}()
            end
        else
            T = mattype == "canonical empty" ? Union{} : str2type_matlab[mattype]
            return Array{T}(undef, dims...)
        end
    end

    objecttype = haskey(dset, object_type_attr_matlab) ? read_attribute(dset, object_type_attr_matlab) : nothing
    mattype = haskey(dset, name_type_attr_matlab) ? read_attribute(dset, name_type_attr_matlab) : "struct_array_field"

    if mattype == "cell" && objecttype === nothing
        # Cell arrays, represented as an array of refs
        return read_cell(dset, subsys)
    elseif objecttype !== nothing
        if objecttype != 3
            @warn "MATLAB Object Type $mattype is currently not supported."
            return missing
        end
        if mattype == "FileWrapper__"
            return read_cell(dset, subsys)
        end
        if haskey(dset, struct_field_attr_matlab)
            @warn "Enumeration Instances are not supported currently."
            return missing
        end
    elseif mattype == "struct_array_field"
        # This will be converted into MatlabStructArray in `m_read(g::HDF5.Group)`
        return StructArrayField(read_cell(dset, subsys))
    elseif !haskey(str2type_matlab,mattype)
        @warn "MATLAB $mattype values are currently not supported."
        return missing
    end

    # Regular arrays of values
    # Convert to Julia type
    if objecttype === nothing
        T = str2type_matlab[mattype]
    else
        T = UInt32 # FIXME: Default for MATLAB objects?
    end

    # Check for a COMPOUND data set, and if so handle complex numbers specially
    dtype = datatype(dset)
    try
        class_id = HDF5.API.h5t_get_class(dtype.id)
        d = class_id == HDF5.API.H5T_COMPOUND ? read_complex(dtype, dset, T) : read(dset, T)
        if objecttype !== nothing
            return MAT_subsys.load_mcos_object(d, "MCOS", subsys)
        else
            return length(d) == 1 ? d[1] : d
        end
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

function read_sparse_matrix(g::HDF5.Group, mattype::String)
    local data
    fn = keys(g)
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
            class_id = HDF5.API.h5t_get_class(dtype.id)
            try
                data = class_id == HDF5.API.H5T_COMPOUND ? read_complex(dtype, dset, T) : read(dset, T)
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
end

function read_struct_as_dict(g::HDF5.Group, subsys::Subsystem)
    if haskey(g, struct_field_attr_matlab)
        fn = [join(f) for f in read_attribute(g, struct_field_attr_matlab)]
    else
        fn = keys(g)
    end
    s = Dict{String, Any}()
    for i = 1:length(fn)
        dset = g[fn[i]]
        try
            s[fn[i]] = m_read(dset, subsys)
        finally
            close(dset)
        end
    end
    return s
end

# reading a struct, struct array, or sparse matrix
function m_read(g::HDF5.Group, subsys::Subsystem)
    if HDF5.name(g) == "/#subsystem#"
        mattype = "#subsystem#"
    else
        mattype = read_attribute(g, name_type_attr_matlab)
    end
    is_object = false
    if mattype != "struct"
        attr = attributes(g)
        # Check if this is a sparse matrix.
        if haskey(attr, sparse_attr_matlab)
            return read_sparse_matrix(g, mattype)
        elseif mattype == "function_handle"
            if haskey(attr, object_decode_attr_matlab) && read_attribute(g, object_decode_attr_matlab)==1
                is_object = true
            end
        else
            if haskey(attr, object_decode_attr_matlab) && read_attribute(g, object_decode_attr_matlab)==2
                # I think this means it's an old object class similar to mXOBJECT_CLASS in MAT_v5
                is_object = true
            elseif mattype != "#subsystem#"
                @warn "Unknown non-struct group of type $mattype detected; attempting to read as struct"
            end
        end
    end
    if is_object
        class = mattype
    else
        class = ""
    end
    s = read_struct_as_dict(g, subsys)
    out = convert_struct_array(s, class)
    return out
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
        val = m_read(obj, f.subsystem)
    finally
        close(obj)
    end
    val
end

"""
    keys(matfile_handle) -> Vector{String}

Return a list of variables in an opened Matlab file.

See `matopen`.
"""
Base.keys(f::MatlabHDF5File) = filter!(x -> x!="#refs#" && x!="#subsystem#", keys(f.plain))

"""
    haskey(matfile_handle, varname) -> Bool

Return true if a variable is present in an opened Matlab file.

See `matopen`.
"""
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
                               deflate = 3, chunk = HDF5.heuristic_chunk(adata))
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
                              deflate = 3, chunk = HDF5.heuristic_chunk(adata))
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

function _normalize_arr(x)
    if x isa Array
        x
    elseif x isa AbstractArray
        collect(x)
    else
        x
    end
end

# Write a scalar or array
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, data::Union{T, Complex{T}, AbstractArray{T}, AbstractArray{Complex{T}}}, ) where {T<:HDF5BitsOrBool}
    data = _normalize_arr(data)
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

# Char
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, c::AbstractChar)
    m_write(mfile, parent, name, string(c))
end

# Tuple
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, t::Tuple)
    m_write(mfile, parent, name, [x for x in t])
end

# Symbol
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, s::Symbol)
    m_write(mfile, parent, name, string(s))
end

# Write cell arrays
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, data::AbstractArray{T}, object_decode::UInt32=UInt32(0)) where T
    data = _normalize_arr(data)
    refs = _write_references(mfile, parent, data)
    # Write the references as the chosen variable
    cset, ctype = create_dataset(parent, name, refs)
    try
        write_dataset(cset, ctype, refs)
        if object_decode == UInt32(3)
            write_attribute(cset, object_decode_attr_matlab, object_decode)
            write_attribute(cset, name_type_attr_matlab, "FileWrapper__")
        else
            write_attribute(cset, name_type_attr_matlab, "cell")
        end
    finally
        close(ctype)
        close(cset)
    end
end

function _write_references(mfile::MatlabHDF5File, parent::HDF5Parent, data::AbstractArray)
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
                close(a)
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
    return refs
end

function _write_field_reference(mfile::MatlabHDF5File, parent::HDF5Parent, k::Vector{String})
    pathrefs = "/#refs#"
    fid = HDF5.file(parent)
    local g
    local ref
    if !haskey(fid, pathrefs)
        g = create_group(fid, pathrefs)
    else
        g = fid[pathrefs]
    end

    try
        mfile.refcounter +=1
        itemname = string(mfile.refcounter)
        cset, ctype = create_dataset(g, itemname, HDF5.VLen(k))
        write_dataset(cset, ctype, HDF5.VLen(k))
        tmp = g[itemname]
        ref = Reference(tmp, pathrefs*"/"*itemname)
        close(tmp)
    finally
        close(g)
    end
    return ref
end

function _write_struct_fields(mfile::MatlabHDF5File, parent::Union{HDF5.Group, HDF5.Dataset}, fieldnames::Vector{String})
    total_chars = sum(length, fieldnames)
    if total_chars < 4096
        write_attribute(parent, struct_field_attr_matlab, HDF5.VLen(fieldnames))
    else
        # Write Reference instead
        ref = _write_field_reference(mfile, parent, fieldnames)
        write_attribute(parent, struct_field_attr_matlab, ref)
    end
end

# Struct array: Array of Dict => MATLAB struct array
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, arr::AbstractArray{<:AbstractDict})
    m_write(mfile, parent, name, MatlabStructArray(arr))
end

# MATLAB struct array
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, arr::MatlabStructArray)
    first_value = first(arr.values)
    if isempty(first_value)
        # write an empty struct array
        adata = [size(first_value)...]
        dset, dtype = create_dataset(parent, name, adata)
        try
            write_attribute(dset, empty_attr_matlab, 0x01)
            write_attribute(dset, name_type_attr_matlab, "struct")
            _write_struct_fields(mfile, dset, arr.names)
            write_dataset(dset, dtype, adata)
        finally
            close(dtype); close(dset)
        end
    else
        g = create_group(parent, name)
        try
            if isempty(arr.class)
                write_attribute(g, name_type_attr_matlab, "struct")
                _write_struct_fields(mfile, g, arr.names)
            else
                write_attribute(g, name_type_attr_matlab, arr.class)
                write_attribute(g, object_decode_attr_matlab, UInt32(2))
                _write_struct_fields(mfile, g, arr.names)
            end
            for (fieldname, field_values) in arr
                refs = _write_references(mfile, parent, field_values)
                dset, dtype = create_dataset(g, fieldname, refs)
                try
                    write_dataset(dset, dtype, refs)
                finally
                    close(dtype); close(dset)
                end
            end
        finally
            close(g)
        end
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

function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, arr::AbstractArray{MatlabClassObject})
    m_write(mfile, parent, name, MatlabStructArray(arr))
end

function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, obj::MatlabClassObject)
    g = create_group(parent, name)
    try
        write_attribute(g, name_type_attr_matlab, obj.class)
        write_attribute(g, object_decode_attr_matlab, UInt32(2))
        all_keys = collect(keys(obj))
        _write_struct_fields(mfile, g, all_keys)
        for (ki, vi) in zip(all_keys, values(obj))
            m_write(mfile, g, ki, vi)
        end
    finally
        close(g)
    end
end

# Write empty (zero-dimensional) structs with no fields
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, s::EmptyStruct)
    dset, dtype = create_dataset(parent, name, s.dims)
    try
        write_attribute(dset, empty_attr_matlab, 0x01)
        write_attribute(dset, name_type_attr_matlab, "struct")
        write_dataset(dset, dtype, s.dims)
    finally
        close(dtype); close(dset)
    end
end

# Write a struct from arrays of keys and values
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, k::Vector{String}, v::Vector)
    if length(k) == 0
        # empty struct
        adata = UInt64[1, 1]
        dset, dtype = create_dataset(parent, name, adata)
        try
            write_attribute(dset, empty_attr_matlab, 0x01)
            write_attribute(dset, name_type_attr_matlab, "struct")
            write_dataset(dset, dtype, adata)
        finally
            close(dtype); close(dset)
        end
        return
    end

    g = create_group(parent, name)
    try
        write_attribute(g, name_type_attr_matlab, "struct")
        for i = 1:length(k)
            m_write(mfile, g, k[i], v[i])
        end
        _write_struct_fields(mfile, g, k)
    finally
        close(g)
    end
end

# Write Associative as a struct
m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, s::AbstractDict) =
    m_write(mfile, parent, name, check_struct_keys(collect(keys(s))), collect(values(s)))

# Write named tuple as a struct
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, nt::NamedTuple)
    m_write(mfile, parent, name, [string(x) for x in keys(nt)], collect(nt))
end

# Write generic CompositeKind as a struct
function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, s)
    if isbits(s)
        error("This is the write function for CompositeKind, but the input doesn't fit")
    elseif Tables.istable(s)
        error("writing tables is not yet supported")
    end
    T = typeof(s)
    m_write(mfile, parent, name, check_struct_keys([string(x) for x in fieldnames(T)]), [getfield(s, x) for x in fieldnames(T)])
end

function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, dat::ScalarOrArray{T}) where T<:Dates.AbstractTime
    error("writing of type $T is not yet supported")
end

function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, dat::ScalarOrArray{T}) where T<:Union{Dates.DateTime, Dates.Millisecond}
    m_write(mfile, parent, name, MatlabOpaque(dat))
end

function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, obj::MatlabOpaque)
    if obj.class == "FileWrapper__"
        m_write(mfile, parent, name, obj["__filewrapper__"], UInt32(3))
        return
    end

    metadata = MAT_subsys.set_mcos_object_metadata(mfile.subsystem, obj)
    dset, dtype = create_dataset(parent, name, metadata)
    try
        write_dataset(dset, dtype, metadata)
        write_attribute(dset, name_type_attr_matlab, obj.class)
        write_attribute(dset, object_type_attr_matlab, UInt32(3))
    finally
        close(dset)
        close(dtype)
    end
end

function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, obj::AbstractArray{MatlabOpaque})
    metadata = MAT_subsys.set_mcos_object_metadata(mfile.subsystem, obj)
    dset, dtype = create_dataset(parent, name, metadata)
    try
        # TODO: Handle empty array case
        write_dataset(dset, dtype, metadata)
        write_attribute(dset, name_type_attr_matlab, first(obj).class)
        write_attribute(dset, object_type_attr_matlab, UInt32(3))
    finally
        close(dset)
        close(dtype)
    end
end

function m_write(mfile::MatlabHDF5File, parent::HDF5Parent, name::String, arr::PooledArray)
    error("writing of PooledArray types as categorical is not yet supported")
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

function write_subsys(mfile::MatlabHDF5File, subsys_data::Dict{String,Any})
    name = "#subsystem#"
    m_write(mfile, mfile.plain, name, subsys_data)
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
    memtype = create_datatype(HDF5.API.H5T_COMPOUND, 2*sizeof(T))
    HDF5.API.h5t_insert(memtype, "real", 0, HDF5.hdf5_type_id(T))
    HDF5.API.h5t_insert(memtype, "imag", sizeof(T), HDF5.hdf5_type_id(T))
    return memtype
end

function check_datatype_complex(dtype::HDF5.Datatype)
    n = HDF5.API.h5t_get_nmembers(dtype.id)
    if n != 2
        return false
    end
    if HDF5.API.h5t_get_member_name(dtype.id, 0) != "real" ||
       HDF5.API.h5t_get_member_name(dtype.id, 1) != "imag"
        return false
    end
    true
end

end
