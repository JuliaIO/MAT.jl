# MAT_subsys.jl
# Tools for processing MAT-file subsystem data in Julia
#
# Copyright (C) 2025    Nithin Lakshmisha
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

module MAT_subsys

import ..MAT_types: MatlabStructArray, MatlabOpaque, convert_opaque

export Subsystem

const FWRAP_VERSION = 4

mutable struct Subsystem
    object_cache::Dict{UInt32, MatlabOpaque}
    num_names::UInt32
    mcos_names::Vector{String}
    class_id_metadata::Vector{UInt32}
    object_id_metadata::Vector{UInt32}
    saveobj_prop_metadata::Vector{UInt32}
    obj_prop_metadata::Vector{UInt32}
    dynprop_metadata::Vector{UInt32}
    _u6_metadata::Vector{UInt32}
    _u7_metadata::Vector{UInt32}
    prop_vals_saved::Vector{Any}
    _c3::Any
    _c2::Any
    prop_vals_defaults::Any
    handle_data::Any
    java_data::Any
    table_type::Type

    Subsystem() = new(
        Dict{UInt32, MatlabOpaque}(),
        UInt32(0),
        String[],
        UInt32[],
        UInt32[],
        UInt32[],
        UInt32[],
        UInt32[],
        UInt32[],
        UInt32[],
        Any[],
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        Nothing
    )
end

function get_object!(subsys::Subsystem, oid::UInt32, classname::String)
    if haskey(subsys.object_cache, oid)
        # object is already cached, just retrieve it
        obj = subsys.object_cache[oid]
    else # it's a new object
        prop_dict = Dict{String,Any}()
        obj = MatlabOpaque(prop_dict, classname)
        # cache the new object
        subsys.object_cache[oid] = obj
        # caching must be done before a next call to `get_properties` to avoid any infinite recursion
        merge!(prop_dict, get_properties(subsys, oid))
    end
    return obj
end

function load_subsys!(subsystem_data::Dict{String,Any}, swap_bytes::Bool)
    subsys = Subsystem()
    load_subsys!(subsys, subsystem_data, swap_bytes)
end

function load_subsys!(subsys::Subsystem, subsystem_data::Dict{String,Any}, swap_bytes::Bool)
    subsys.handle_data = get(subsystem_data, "handle", nothing)
    subsys.java_data = get(subsystem_data, "java", nothing)
    mcos_data = get(subsystem_data, "MCOS", nothing)
    if mcos_data === nothing
        return
    end

    if mcos_data isa Tuple
        # Backward compatibility with MAT_v5
        mcos_data = mcos_data[2]
    end
    fwrap_metadata = vec(mcos_data[1, 1])

    # FIXME: Is this the best way to read?
    # Integers are written as uint8 (with swap), interpret as uint32
    version = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[1:4]) : fwrap_metadata[1:4])[1]
    if version <= 1 || version > FWRAP_VERSION
        error("Cannot read subsystem: Unsupported FileWrapper version: $version")
    end

    subsys.num_names = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[5:8]) : fwrap_metadata[5:8])[1]
    region_offsets = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[9:40]) : fwrap_metadata[9:40])

    # Class and Property Names stored as list of null-terminated strings
    start = 41
    pos = start
    name_count = 0
    while name_count < subsys.num_names
        if fwrap_metadata[pos] == 0x00
            push!(subsys.mcos_names, String(fwrap_metadata[start:pos-1]))
            name_count += 1
            start = pos + 1
            if name_count == subsys.num_names
                break
            end
        end
        pos += 1
    end

    subsys.class_id_metadata = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[region_offsets[1]+1:region_offsets[2]]) : fwrap_metadata[region_offsets[1]+1:region_offsets[2]])
    subsys.saveobj_prop_metadata = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[region_offsets[2]+1:region_offsets[3]]) : fwrap_metadata[region_offsets[2]+1:region_offsets[3]])
    subsys.object_id_metadata = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[region_offsets[3]+1:region_offsets[4]]) : fwrap_metadata[region_offsets[3]+1:region_offsets[4]])
    subsys.obj_prop_metadata = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[region_offsets[4]+1:region_offsets[5]]) : fwrap_metadata[region_offsets[4]+1:region_offsets[5]])
    subsys.dynprop_metadata = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[region_offsets[5]+1:region_offsets[6]]) : fwrap_metadata[region_offsets[5]+1:region_offsets[6]])

    if region_offsets[7] != 0
        subsys._u6_metadata = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[region_offsets[6]+1:region_offsets[7]]) : fwrap_metadata[region_offsets[6]+1:region_offsets[7]])
    end

    if region_offsets[8] != 0
        subsys._u7_metadata = reinterpret(UInt32, swap_bytes ? reverse(fwrap_metadata[region_offsets[7]+1:region_offsets[8]]) : fwrap_metadata[region_offsets[7]+1:region_offsets[8]])
    end

    if version == 2
        subsys.prop_vals_saved = mcos_data[3:end-1, 1]
    elseif version == 3
        subsys.prop_vals_saved = mcos_data[3:end-2, 1]
        subsys._c2 = mcos_data[end-1, 1]
    else
        subsys.prop_vals_saved = mcos_data[3:end-3, 1]
        subsys._c3 = mcos_data[end-2, 1]
    end

    subsys.prop_vals_defaults = mcos_data[end, 1]
    for el in subsys.prop_vals_defaults
        update_nested_props!(el, subsys) # just in case
    end

    return subsys
end

function get_classname(subsys::Subsystem, class_id::UInt32)
    namespace_idx = subsys.class_id_metadata[class_id*4+1]
    classname_idx = subsys.class_id_metadata[class_id*4+2]

    namespace = if namespace_idx == 0
        ""
    else
        subsys.mcos_names[namespace_idx] * "."
    end

    classname = namespace * subsys.mcos_names[classname_idx]
    return classname
end

function get_object_metadata(subsys::Subsystem, object_id::UInt32)
    return subsys.object_id_metadata[object_id*6+1:object_id*6+6]
end

function get_default_properties(subsys::Subsystem, class_id::UInt32)
    return Dict{String,Any}(subsys.prop_vals_defaults[class_id+1, 1])
end

function get_property_idxs(subsys::Subsystem, obj_type_id::UInt32, saveobj_ret_type::Bool)
    prop_field_idxs = saveobj_ret_type ? subsys.saveobj_prop_metadata : subsys.obj_prop_metadata
    nfields = 3
    offset = 1
    while obj_type_id > 0
        nprops = prop_field_idxs[offset]
        offset += 1 + (nfields * nprops)
        offset += (offset + 1) % 2  # Padding
        obj_type_id -= 1
    end
    nprops = prop_field_idxs[offset]
    offset += 1
    return prop_field_idxs[offset:offset+nprops*nfields-1]
end

update_nested_props!(prop_value, subsys::Subsystem) = prop_value

function update_nested_props!(prop_value::Union{AbstractDict, MatlabStructArray}, subsys::Subsystem)
    # Handle nested objects in structs
    for (key, value) in prop_value
        prop_value[key] = update_nested_props!(value, subsys)
    end
    return prop_value
end

function update_nested_props!(prop_value::Array{Any}, subsys::Subsystem)
    # Handle nested objects in a Cell
    for i in eachindex(prop_value)
        prop_value[i] = update_nested_props!(prop_value[i], subsys)
    end
    return prop_value
end

function update_nested_props!(prop_value::Array{UInt32}, subsys::Subsystem)
    # Hacky way to find and update nested objects
    # Nested objects are stored as a uint32 Matrix with a unique signature
    # MATLAB probably uses some kind of placeholders to decode
    # But this should work here

    if first(prop_value) == 0xdd000000
        # MATLAB identifies any uint32 array with first value 0xdd000000 as an MCOS object
        return load_mcos_object(prop_value, "MCOS", subsys)
    else
        return prop_value
    end
end

function get_saved_properties(subsys::Subsystem, obj_type_id::UInt32, saveobj_ret_type::Bool)
    save_prop_map = Dict{String,Any}()
    prop_field_idxs = get_property_idxs(subsys, obj_type_id, saveobj_ret_type)
    nprops = length(prop_field_idxs) ÷ 3
    for i in 0:nprops-1
        prop_name = subsys.mcos_names[prop_field_idxs[i*3+1]]
        prop_type = prop_field_idxs[i*3+2]
        if prop_type == 0
            prop_value = subsys.mcos_names[prop_field_idxs[i*3+3]]
        elseif prop_type == 1
            prop_value = subsys.prop_vals_saved[prop_field_idxs[i*3+3]+1]
        elseif prop_type == 2
            prop_value = prop_field_idxs[i*3+3]
        else
            error("Unknown property type ID: $prop_type encountered during deserialization")
        end
        save_prop_map[prop_name] = update_nested_props!(prop_value, subsys)
    end
    return save_prop_map
end

function get_properties(subsys::Subsystem, object_id::UInt32)
    if object_id == 0
        return Dict{String,Any}()
    end

    class_id, _, _, saveobj_id, normobj_id, _ = get_object_metadata(subsys, object_id)
    if saveobj_id != 0
        saveobj_ret_type = true
        obj_type_id = saveobj_id
    else
        saveobj_ret_type = false
        obj_type_id = normobj_id
    end

    defaults = get_default_properties(subsys, class_id)
    prop_map = merge(defaults, get_saved_properties(subsys, obj_type_id, saveobj_ret_type))
    # TODO: Add dynamic properties
    return prop_map
end

function load_mcos_object(metadata::Any, type_name::String, subsys::Subsystem)
    @warn "Expected MCOS metadata to be an Array{UInt32}, got $(typeof(metadata)). Returning metadata."
    return metadata
end

function load_mcos_object(metadata::Dict, type_name::String, subsys::Subsystem)
    @warn "Loading enumeration instances are not supported. Returning Metadata"
    return metadata
end

function load_mcos_object(metadata::Array{UInt32}, type_name::String, subsys::Subsystem)
    if type_name != "MCOS"
        @warn "Loading Type:$type_name is not implemented. Returning metadata."
        return metadata
    end

    if metadata[1, 1] != 0xDD000000
        @warn "MCOS object metadata is corrupted. Returning raw data."
        return metadata
    end

    ndims = metadata[2, 1]
    dims = metadata[3:2+ndims, 1]
    nobjects = prod(dims)
    object_ids = metadata[3+ndims:2+ndims+nobjects, 1]

    class_id = metadata[end, 1]
    classname = get_classname(subsys, class_id)

    if nobjects == 1
        oid = object_ids[1]
        obj = get_object!(subsys, oid, classname)
        return convert_opaque(obj; table=subsys.table_type)
    else
        object_arr = Array{Any}(undef, convert(Vector{Int}, dims)...)
        for i = 1:length(object_arr)
            oid = object_ids[i]
            obj = get_object!(subsys, oid, classname)
            object_arr[i] = convert_opaque(obj; table=subsys.table_type)
        end
        return object_arr
    end
end

end