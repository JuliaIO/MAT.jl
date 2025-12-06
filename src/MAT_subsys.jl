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

import ..MAT_types: MatlabStructArray, MatlabOpaque, convert_opaque, EmptyStruct

export Subsystem

const FWRAP_VERSION = 4
const MCOS_IDENTIFIER = 0xdd000000

const matlab_saveobj_ret_types = String[
    "string",
    "timetable"
]

mutable struct Subsystem
    load_object_cache::Dict{UInt32,MatlabOpaque}
    save_object_cache::IdDict{MatlabOpaque,UInt32}
    num_names::UInt32 # number of mcos_names
    mcos_names::Vector{String} # Class and Property Names
    class_id_metadata::Vector{UInt32}
    object_id_metadata::Vector{Vector{UInt32}}
    saveobj_prop_metadata::Vector{Vector{UInt32}}
    obj_prop_metadata::Vector{Vector{UInt32}}
    dynprop_metadata::Vector{UInt32}
    _u6_metadata::Vector{UInt32}
    _u7_metadata::Vector{UInt32}
    prop_vals_saved::Vector{Any}
    _c3::Any
    mcos_class_alias_metadata::Any
    prop_vals_defaults::Any
    handle_data::Any
    java_data::Any

    # automatic MatlabOpaque conversion
    convert_opaque::Bool
    table_type::Type

    # Counters for saving
    saveobj_counter::UInt32
    normalobj_counter::UInt32
    obj_id_counter::UInt32
    class_id_counter::UInt32

    function Subsystem()
        return new(
            Dict{UInt32,MatlabOpaque}(),
            IdDict{MatlabOpaque,UInt32}(),
            UInt32(0),
            String[],
            UInt32[],
            Vector{Vector{UInt32}}(),
            Vector{Vector{UInt32}}(),
            Vector{Vector{UInt32}}(),
            UInt32[],
            UInt32[],
            UInt32[],
            Any[],
            nothing,
            Int32[],
            nothing,
            nothing,
            nothing,
            true,
            Nothing,
            UInt32(0),
            UInt32(0),
            UInt32(0),
            UInt32(0),
        )
    end
end

function swapped_reinterpret(T::Type, A::AbstractArray{UInt8}, swap_bytes::Bool)
    return reinterpret(T, swap_bytes ? reverse(A) : A)
end
# integers are written as uint8 (with swap), interpret as uint32
function swapped_reinterpret(A::AbstractArray{UInt8}, swap_bytes::Bool)
    return swapped_reinterpret(UInt32, A, swap_bytes)
end

function init_save!(subsys::Subsystem)
    append!(subsys.class_id_metadata, UInt32[0, 0, 0, 0])
    append!(subsys.dynprop_metadata, UInt32[0, 0])
    append!(subsys.mcos_class_alias_metadata, Int32[0])

    # These are Vector{Vector{uint32}}
    # need mutable inner vectors to handle nested properties
    push!(subsys.object_id_metadata, UInt32[0, 0, 0, 0, 0, 0])
    push!(subsys.saveobj_prop_metadata, UInt32[0, 0])
    push!(subsys.obj_prop_metadata, UInt32[0, 0])

    subsys._c3 = Any[]
    subsys.prop_vals_defaults = Any[]

    return subsys
end

function load_subsys!(subsystem_data::Dict{String,Any}, swap_bytes::Bool)
    subsys = Subsystem()
    return load_subsys!(subsys, subsystem_data, swap_bytes)
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
    subsys.saveobj_prop_metadata = [swapped_reinterpret(
        get_region(fwrap_metadata, region_offsets, 2), swap_bytes
    )]
    subsys.object_id_metadata = [swapped_reinterpret(
        get_region(fwrap_metadata, region_offsets, 3), swap_bytes
    )]
    subsys.obj_prop_metadata = [swapped_reinterpret(
        get_region(fwrap_metadata, region_offsets, 4), swap_bytes
    )]
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

function load_subsys!(subsys::Subsystem, subsystem_data::Dict{String,Any}, swap_bytes::Bool)
    try
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
            subsys.mcos_class_alias_metadata = mcos_data[end - 1, 1]
        else
            subsys.prop_vals_saved = mcos_data[3:(end - 3), 1]
            subsys._c3 = mcos_data[end - 2, 1]
        end

        subsys.prop_vals_defaults = mcos_data[end, 1]
        return subsys
    catch
        @warn "Failed to load MAT-file subsystem data. Opaque objects will be skipped. Error: $(catch_backtrace())"
        return subsys
    end
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
    return subsys.object_id_metadata[1][(object_id * 6 + 1):(object_id * 6 + 6)]
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

function get_default_properties(subsys::Subsystem, class_id::UInt32)
    default_props = Dict{String,Any}(subsys.prop_vals_defaults[class_id + 1, 1])
    for (key, value) in default_props
        default_props[key] = update_nested_props!(value, subsys)
    end
    return default_props
end

function get_property_idxs(subsys::Subsystem, obj_type_id::UInt32, saveobj_ret_type::Bool)
    prop_field_idxs =
        saveobj_ret_type ? subsys.saveobj_prop_metadata[1] : subsys.obj_prop_metadata[1]
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

    if length(subsys.dynprop_metadata) == 0
        return Dict{String,Any}()
    end

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

function get_object!(subsys::Subsystem, oid::UInt32, classname::String)
    if haskey(subsys.load_object_cache, oid)
        # object is already cached, just retrieve it
        obj = subsys.load_object_cache[oid]
    else # it's a new object
        prop_dict = Dict{String,Any}()
        obj = MatlabOpaque(prop_dict, classname)
        # cache the new object
        subsys.load_object_cache[oid] = obj
        # caching must be done before a next call to `get_properties` to avoid any infinite recursion
        merge!(prop_dict, get_properties(subsys, oid))
    end
    return obj
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
    try
        if length(subsys.class_id_metadata) == 0
            # Subsystem was not loaded properly
            # No need to warn again
            return metadata
        end

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
            if subsys.convert_opaque
                return convert_opaque(obj; table=subsys.table_type)
            else
                return obj
            end
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

    catch e
        @warn "Failed to load MCOS object. Returning raw metadata. Error: $e"
        return metadata
    end
end

function set_mcos_name!(subsys::Subsystem, name::String)
    idx = findfirst(isequal(name), subsys.mcos_names)
    if idx !== nothing
        return UInt32(idx)
    else
        push!(subsys.mcos_names, name)
        subsys.num_names += UInt32(1)
        return subsys.num_names
    end
end

function check_valid_property_name(s::AbstractString)
    if match(r"^[a-zA-Z][a-zA-Z0-9_]*$", s) === nothing
        error("Invalid property name \"$s\": property names must start with a letter and contain only alphanumeric characters and underscore")
    end
end

function set_class_id!(subsys::Subsystem, classname::String)
    # number of existing class entries (each class has 4 UInt32 metadata entries)
    class_count = UInt32(length(subsys.class_id_metadata) ÷ 4 - 1) # skip class_id = 0 case

    # Check if class name already exists
    for cid in UInt32(1):UInt32(class_count)
        existing_name = get_classname(subsys, cid)
        if existing_name == classname
            return cid
        end
    end

    # Add new class metadata
    subsys.class_id_counter += UInt32(1)
    # split namespace and class name at last dot
    pos = findlast(==('.'), classname)
    if pos === nothing
        namespace = ""
        cname = classname
    else
        namespace = classname[1:pos-1]
        cname = classname[pos+1:end]
    end

    cname_idx = set_mcos_name!(subsys, cname)
    namespace_idx = pos === nothing ? UInt32(0) : set_mcos_name!(subsys, namespace)

    append!(subsys.class_id_metadata, UInt32[namespace_idx, cname_idx, 0, 0])
    append!(subsys.mcos_class_alias_metadata, Int32[0]) # Placeholder for no aliases

    return subsys.class_id_counter
end

save_nested_props(prop_value, subsys::Subsystem) = prop_value

function save_nested_props(
    prop_value::Union{AbstractDict,MatlabStructArray}, subsys::Subsystem
)
    # Save nested objects in structs
    for (key, value) in prop_value
        prop_value[key] = save_nested_props(value, subsys)
    end
    return prop_value
end

function save_nested_props(prop_value::Array{Any}, subsys::Subsystem)
    # Save nested objects in a Cell
    for i in eachindex(prop_value)
        prop_value[i] = save_nested_props(prop_value[i], subsys)
    end
    return prop_value
end

function save_nested_props(prop_value::Union{MatlabOpaque, Array{MatlabOpaque}}, subsys::Subsystem)
    # Nested objects are saved by the uint32 Matrix signature
    # however they don't have any mxOPAQUE_CLASS headers
    # so we search within containers here

    # FIXME: Does this overwrite prop_value from the user dict?
    # Might have to create a copy instead - Test needed
    prop_value = set_mcos_object_metadata(subsys, prop_value)
    return prop_value
end

function serialize_object_props!(subsys::Subsystem, obj::MatlabOpaque, obj_prop_metadata::Vector{UInt32})

    object_id = subsys.obj_id_counter
    nprops = length(obj)

    # property metadata format: [num_props, (prop_name_idx, prop_type, prop_value_idx)*]
    push!(obj_prop_metadata, UInt32(nprops))

    for (prop_name, prop_value) in obj
        if startswith(prop_name, "__dynamic_property_")
            @warn "Dynamic properties are not supported when writing: skipping $prop_name"
            continue
        end

        check_valid_property_name(prop_name)

        field_name_idx = set_mcos_name!(subsys, prop_name)
        prop_value = save_nested_props(prop_value, subsys)

        prop_vals = UInt32[field_name_idx, 1, 0]

        cell_idx = length(subsys.prop_vals_saved) # these are zero-indexed in matlab
        push!(subsys.prop_vals_saved, prop_value)
        prop_vals[3] = cell_idx
        append!(obj_prop_metadata, prop_vals)
    end

    if length(obj_prop_metadata) % 2 == 1
        push!(obj_prop_metadata, UInt32(0)) # padding
    end

    # placeholder for dynamic props
    append!(subsys.dynprop_metadata, UInt32[0, 0])

    ndeps = subsys.obj_id_counter - object_id
    return ndeps
end

function set_object_id(subsys::Subsystem, obj::MatlabOpaque, saveobj_ret_type=false)

    if length(obj) == 0
        # This is a deleted object
        # MATLAB keeps weak references to deleted objects for some reason
        class_id = set_class_id!(subsys, obj.class)
        obj_id = 0
        return obj_id, class_id
    end

    if haskey(subsys.save_object_cache, obj)
        class_id = set_class_id!(subsys, obj.class)
        obj_id = subsys.save_object_cache[obj]
        return obj_id, class_id
    end

    subsys.obj_id_counter += UInt32(1)
    mat_obj_id = subsys.obj_id_counter
    subsys.save_object_cache[obj] = mat_obj_id

    prop_metadata = UInt32[]
    saveobj_id = 0
    normobj_id = 0

    if saveobj_ret_type
        subsys.saveobj_counter += UInt32(1)
        saveobj_id = subsys.saveobj_counter
        push!(subsys.saveobj_prop_metadata, prop_metadata)
    else
        subsys.normalobj_counter += UInt32(1)
        normobj_id = subsys.normalobj_counter
        push!(subsys.obj_prop_metadata, prop_metadata)
    end

    obj_id_metadata = zeros(UInt32, 6)
    push!(subsys.object_id_metadata, obj_id_metadata)

    if saveobj_ret_type && !(length(obj) == 1 && haskey(obj, "any"))
        error("Object of class $(obj.classname) marked with a saveobj return type must have a single property 'any' containing the return value of its saveobj method")
    end

    ndeps = serialize_object_props!(subsys, obj, prop_metadata)

    class_id = set_class_id!(subsys, obj.class)
    obj_id_metadata[1] = class_id
    if saveobj_ret_type
        obj_id_metadata[4] = saveobj_id
    else
        obj_id_metadata[5] = normobj_id
    end

    # it works idk how
    obj_id_metadata[6] = mat_obj_id + ndeps
    for i in 1:ndeps
        obj_id = subsys.obj_id_counter - ndeps + i
        subsys.object_id_metadata[obj_id + 1][6] -= 1
    end

    return mat_obj_id, class_id
end

function create_mcos_metadata_array(dims::Tuple{Vararg{Int}}, arr_ids::Vector{UInt32}, class_id::UInt32)
    ndims = length(dims)
    metadata = UInt32[]
    push!(metadata, MCOS_IDENTIFIER)
    push!(metadata, UInt32(ndims))
    for d in dims
        push!(metadata, UInt32(d))
    end
    for oid in arr_ids
        push!(metadata, oid)
    end
    push!(metadata, class_id)

    return reshape(metadata, :, 1)
end

function set_mcos_object_metadata(subsys::Subsystem, obj::Array{MatlabOpaque})

    arr_ids = UInt32[]
    if length(obj) == 0
        # TODO: Handle 1x0, 0x0, 0x1 objects
        # placeholder for empty array case
        dims = size(obj)
        return create_mcos_metadata_array(dims, UInt32(0), UInt32(0))
    end
    classname = first(obj).class

    # this is not needed but added for consistency
    # future update can support mentioning classes with saveobj methods
    saveobj_ret_type = classname in matlab_saveobj_ret_types

    class_id = UInt32(0)
    for obj_elem in obj
        obj_id, class_id = set_object_id(subsys, obj_elem, saveobj_ret_type)
        append!(arr_ids, obj_id)
    end

    dims = size(obj)
    return create_mcos_metadata_array(dims, arr_ids, class_id)
end

function set_mcos_object_metadata(subsys::Subsystem, obj::MatlabOpaque)

    classname = obj.class
    saveobj_ret_type = classname in matlab_saveobj_ret_types

    object_id, class_id = set_object_id(subsys, obj, saveobj_ret_type)
    dims = (1, 1)
    return create_mcos_metadata_array(dims, [object_id], class_id)
end

function set_fwrap_metadata!(subsys::Subsystem)

    version_bytes = vec(reinterpret(UInt8, UInt32[FWRAP_VERSION]))
    num_names_bytes = copy(vec(reinterpret(UInt8, UInt32[subsys.num_names])))

    region_offsets = zeros(UInt32, 8)

    names_str = isempty(subsys.mcos_names) ? "" : join(subsys.mcos_names, '\0') * '\0'
    names_bytes = collect(codeunits(names_str))
    pad_len = (8 - (length(names_bytes) % 8)) % 8
    if pad_len > 0
        append!(names_bytes, zeros(UInt8, pad_len))
    end

    region1_bytes = vec(reinterpret(UInt8, subsys.class_id_metadata))
    region_offsets[1] = UInt32(40 + length(names_bytes))
    region_offsets[2] = region_offsets[1] + UInt32(length(region1_bytes))

    # This region (including id=0) does not exist if there are no saveobj type objects
    region2_bytes = UInt8[]
    if length(subsys.saveobj_prop_metadata) > 1
        for sub in subsys.saveobj_prop_metadata
            append!(region2_bytes, vec(reinterpret(UInt8, sub)))
        end
    end
    region_offsets[3] = region_offsets[2] + UInt32(length(region2_bytes))

    region3_bytes = UInt8[]
    for sub in subsys.object_id_metadata
        append!(region3_bytes, vec(reinterpret(UInt8, sub)))
    end
    region_offsets[4] = region_offsets[3] + UInt32(length(region3_bytes))

    region4_bytes = UInt8[]
    for sub in subsys.obj_prop_metadata
        append!(region4_bytes, vec(reinterpret(UInt8, sub)))
    end
    region_offsets[5] = region_offsets[4] + UInt32(length(region4_bytes))

    region5_bytes = vec(reinterpret(UInt8, subsys.dynprop_metadata))
    region_offsets[6] = region_offsets[5] + UInt32(length(region5_bytes))

    region6_bytes = UInt8[]
    region_offsets[7] = region_offsets[6] + UInt32(length(region6_bytes))

    region7_bytes = zeros(UInt8, 8)
    region_offsets[8] = region_offsets[7] + UInt32(length(region7_bytes))

    fwrap = Vector{UInt8}()
    append!(fwrap, version_bytes)
    append!(fwrap, num_names_bytes)
    append!(fwrap, vec(reinterpret(UInt8, region_offsets)))
    append!(fwrap, names_bytes)
    append!(fwrap, region1_bytes)
    append!(fwrap, region2_bytes)
    append!(fwrap, region3_bytes)
    append!(fwrap, region4_bytes)
    append!(fwrap, region5_bytes)
    append!(fwrap, region6_bytes)
    append!(fwrap, region7_bytes)

    return fwrap
end

function set_fwrap_data!(subsys::Subsystem)

    fwrap_data = Any[]
    fwrap_metadata = set_fwrap_metadata!(subsys)
    push!(fwrap_data, reshape(fwrap_metadata, :, 1))
    push!(fwrap_data, reshape(Any[], 0, 0))
    append!(fwrap_data, subsys.prop_vals_saved)

    empty_struct = EmptyStruct([1, 0])
    for i in 0:subsys.class_id_counter
        push!(subsys._c3, empty_struct)
        push!(subsys.prop_vals_defaults, empty_struct)
    end

    push!(fwrap_data, reshape(subsys._c3, :, 1))
    push!(fwrap_data, reshape(subsys.mcos_class_alias_metadata, :, 1))
    push!(fwrap_data, reshape(subsys.prop_vals_defaults, :, 1))

    fw_obj = MatlabOpaque(Dict{String,Any}("__filewrapper__" => reshape(fwrap_data, :, 1)), "FileWrapper__")
    return fw_obj
end

function set_subsystem_data!(subsys::Subsystem)
    if subsys.class_id_counter == 0
        # No MCOS objects to serialize
        return nothing
    end

    fwrap = set_fwrap_data!(subsys)
    subsys_struct = Dict{String,Any}()
    subsys_struct["MCOS"] = fwrap
    return subsys_struct
end

end
