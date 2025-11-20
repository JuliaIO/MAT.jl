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

# For reference
# https://github.com/foreverallama/matio/blob/main/docs/subsystem_data_format.md

module MAT_subsys

import ..MAT_types: MatlabStructArray, MatlabOpaque, convert_opaque

export Subsystem

const FWRAP_VERSION = 4
const MCOS_IDENTIFIER = 0xdd000000

mutable struct Subsystem
    object_cache::Dict{UInt32,MatlabOpaque}
    num_names::UInt32 # number of mcos_names
    mcos_names::Vector{String} # Class and Property Names
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
    table_type::Type # Julia type to convert Matlab tables into

    function Subsystem()
        return new(
            Dict{UInt32,MatlabOpaque}(),
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
            Nothing,
        )
    end
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
    return load_subsys!(subsys, subsystem_data, swap_bytes)
end

function load_subsys!(subsys::Subsystem, subsystem_data::Dict{String,Any}, swap_bytes::Bool)
    subsys.handle_data = get(subsystem_data, "handle", nothing)
    subsys.java_data = get(subsystem_data, "java", nothing)
    mcos_data = get(subsystem_data, "MCOS", nothing)
    if mcos_data === nothing
        return nothing
    end

    if mcos_data isa Tuple
        # Backward compatibility with MAT_v5
        mcos_data = mcos_data[2]
    end
    fwrap_metadata::Vector{UInt8} = vec(mcos_data[1, 1])

    version = swapped_reinterpret(fwrap_metadata[1:4], swap_bytes)[1]
    if version <= 1 || version > FWRAP_VERSION
        error("Cannot read subsystem: Unsupported FileWrapper version: $version")
    end

    subsys.num_names = swapped_reinterpret(fwrap_metadata[5:8], swap_bytes)[1]
    load_mcos_names!(subsys, fwrap_metadata)

    load_mcos_regions!(subsys, fwrap_metadata, swap_bytes)

    if version == 2
        subsys.prop_vals_saved = mcos_data[3:(end - 1), 1]
    elseif version == 3
        subsys.prop_vals_saved = mcos_data[3:(end - 2), 1]
        subsys._c2 = mcos_data[end - 1, 1]
    else
        subsys.prop_vals_saved = mcos_data[3:(end - 3), 1]
        subsys._c3 = mcos_data[end - 2, 1]
    end

    subsys.prop_vals_defaults = mcos_data[end, 1]
    return subsys
end

# Class and Property Names are stored as list of null-terminated strings
function load_mcos_names!(subsys::Subsystem, fwrap_metadata::AbstractArray{UInt8})
    start = 41
    pos = start
    name_count = 0
    while name_count < subsys.num_names
        if fwrap_metadata[pos] == 0x00
            push!(subsys.mcos_names, String(fwrap_metadata[start:(pos - 1)]))
            name_count += 1
            start = pos + 1
            if name_count == subsys.num_names
                break
            end
        end
        pos += 1
    end
end

function load_mcos_regions!(
    subsys::Subsystem, fwrap_metadata::AbstractArray{UInt8}, swap_bytes::Bool
)
    region_offsets = swapped_reinterpret(fwrap_metadata[9:40], swap_bytes)

    subsys.class_id_metadata = swapped_reinterpret(
        get_region(fwrap_metadata, region_offsets, 1), swap_bytes
    )
    subsys.saveobj_prop_metadata = swapped_reinterpret(
        get_region(fwrap_metadata, region_offsets, 2), swap_bytes
    )
    subsys.object_id_metadata = swapped_reinterpret(
        get_region(fwrap_metadata, region_offsets, 3), swap_bytes
    )
    subsys.obj_prop_metadata = swapped_reinterpret(
        get_region(fwrap_metadata, region_offsets, 4), swap_bytes
    )
    subsys.dynprop_metadata = swapped_reinterpret(
        get_region(fwrap_metadata, region_offsets, 5), swap_bytes
    )

    if region_offsets[7] != 0
        subsys._u6_metadata = swapped_reinterpret(
            get_region(fwrap_metadata, region_offsets, 6), swap_bytes
        )
    end

    if region_offsets[8] != 0
        subsys._u7_metadata = swapped_reinterpret(
            get_region(fwrap_metadata, region_offsets, 7), swap_bytes
        )
    end
end

function get_region(
    fwrap_metadata::Vector{UInt8}, region_offsets::AbstractVector{UInt32}, region::Integer
)
    return fwrap_metadata[(region_offsets[region] + 1):region_offsets[region + 1]]
end

function swapped_reinterpret(T::Type, A::AbstractArray{UInt8}, swap_bytes::Bool)
    return reinterpret(T, swap_bytes ? reverse(A) : A)
end
# integers are written as uint8 (with swap), interpret as uint32
function swapped_reinterpret(A::AbstractArray{UInt8}, swap_bytes::Bool)
    return swapped_reinterpret(UInt32, A, swap_bytes)
end

function get_classname(subsys::Subsystem, class_id::UInt32)
    namespace_idx = subsys.class_id_metadata[class_id * 4 + 1]
    classname_idx = subsys.class_id_metadata[class_id * 4 + 2]

    namespace = if namespace_idx == 0
        ""
    else
        subsys.mcos_names[namespace_idx] * "."
    end

    classname = namespace * subsys.mcos_names[classname_idx]
    return classname
end

function get_object_metadata(subsys::Subsystem, object_id::UInt32)
    return subsys.object_id_metadata[(object_id * 6 + 1):(object_id * 6 + 6)]
end

function get_default_properties(subsys::Subsystem, class_id::UInt32)
    default_props = Dict{String,Any}(subsys.prop_vals_defaults[class_id + 1, 1])
    for (key, value) in default_props
        default_props[key] = update_nested_props!(value, subsys)
    end
    return default_props
end

function get_property_idxs(subsys::Subsystem, obj_type_id::UInt32, saveobj_ret_type::Bool)
    prop_field_idxs =
        saveobj_ret_type ? subsys.saveobj_prop_metadata : subsys.obj_prop_metadata
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
    return prop_field_idxs[offset:(offset + nprops * nfields - 1)]
end

update_nested_props!(prop_value, subsys::Subsystem) = prop_value

function update_nested_props!(
    prop_value::Union{AbstractDict,MatlabStructArray}, subsys::Subsystem
)
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

    if first(prop_value) == MCOS_IDENTIFIER
        # MATLAB identifies any uint32 array with first value 0xdd000000 as an MCOS object
        return load_mcos_object(prop_value, "MCOS", subsys)
    else
        return prop_value
    end
end

function get_saved_properties(
    subsys::Subsystem, obj_type_id::UInt32, saveobj_ret_type::Bool
)
    save_prop_map = Dict{String,Any}()
    prop_field_idxs = get_property_idxs(subsys, obj_type_id, saveobj_ret_type)
    nprops = length(prop_field_idxs) ÷ 3
    for i in 0:(nprops - 1)
        prop_name = subsys.mcos_names[prop_field_idxs[i * 3 + 1]]
        prop_type = prop_field_idxs[i * 3 + 2]
        if prop_type == 0
            prop_value = subsys.mcos_names[prop_field_idxs[i * 3 + 3]]
        elseif prop_type == 1
            prop_value = subsys.prop_vals_saved[prop_field_idxs[i * 3 + 3] + 1]
        elseif prop_type == 2
            prop_value = prop_field_idxs[i * 3 + 3]
        else
            @warn "Unknown property type ID: $prop_type for property $prop_name encountered during deserialization"
            prop_value = prop_field_idxs[i * 3 + 3]
        end
        save_prop_map[prop_name] = update_nested_props!(prop_value, subsys)
    end
    return save_prop_map
end

function get_dynamic_properties(subsys::Subsystem, dep_id::UInt32)
    offset = 1
    while dep_id > 0
        nprops = subsys.dynprop_metadata[offset]
        offset += 1 + nprops
        offset += (offset + 1) % 2  # Padding
        dep_id -= 1
    end

    ndynprops = subsys.dynprop_metadata[offset]
    offset += 1
    dyn_prop_obj_ids = subsys.dynprop_metadata[offset:(offset + ndynprops - 1)]

    if dyn_prop_obj_ids == UInt32[]
        return Dict{String,Any}()
    end
    dyn_prop_map = Dict{String,Any}()
    for (i, obj_id) in enumerate(dyn_prop_obj_ids)
        dyn_class_id = get_object_metadata(subsys, obj_id)[1]
        classname = get_classname(subsys, dyn_class_id)
        dynobj_props = Dict{String,Any}()
        dynobj = MatlabOpaque(dynobj_props, classname)
        merge!(dynobj_props, get_properties(subsys, obj_id))
        dyn_prop_map["__dynamic_property_$(i)__"] = dynobj
    end
    return dyn_prop_map
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
    dyn_props = get_dynamic_properties(subsys, object_id)
    merge!(prop_map, dyn_props)
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

    if metadata[1, 1] != MCOS_IDENTIFIER
        @warn "MCOS object metadata is corrupted. Returning raw data."
        return metadata
    end

    ndims = metadata[2, 1]
    dims = metadata[3:(2 + ndims), 1]
    nobjects = prod(dims)
    object_ids = metadata[(3 + ndims):(2 + ndims + nobjects), 1]

    class_id = metadata[end, 1]
    classname = get_classname(subsys, class_id)

    if nobjects == 1
        oid = object_ids[1]
        obj = get_object!(subsys, oid, classname)
        return convert_opaque(obj; table=subsys.table_type)
    else
        # no need to convert_opaque, matlab wraps object arrays in a single class normally
        object_arr = Array{MatlabOpaque}(undef, convert(Vector{Int}, dims)...)
        for i in 1:length(object_arr)
            oid = object_ids[i]
            obj = get_object!(subsys, oid, classname)
            object_arr[i] = obj
        end
        return object_arr
    end
end

end