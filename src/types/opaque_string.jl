# opaque_string.jl
# Conversion of the MATLAB `string` class

# the character count MATLAB stores for a <missing> string
const MISSING_STRING_COUNT = typemax(UInt64)

# for reference: https://github.com/foreverallama/matio/blob/main/matio/utils/converters/matstring.py
function from_string(obj::MatlabOpaque, encoding::Encoding=Encoding(Symbol("UTF-16LE")))
    data = obj["any"]
    if isnothing(data) || isempty(data)
        return String[]
    end
    if data[1, 1] != 1
        @warn "String saved from a different MAT-file version. Returning empty string"
        return ""
    end
    ndims = data[1, 2]
    shape = Int.(data[1, 3:(2 + ndims)])
    num_strings = prod(shape)
    char_counts = data[1, (3 + ndims):(2 + ndims + num_strings)]
    byte_data = data[1, (3 + ndims + num_strings):end]
    bytes = reinterpret(UInt8, byte_data)

    strings = Union{Missing,String}[]
    pos = 1

    for char_count in char_counts
        if char_count == MISSING_STRING_COUNT
            # <missing> has no characters stored
            push!(strings, missing)
            continue
        end
        byte_length = char_count * 2  # UTF-16 encoding
        extracted_bytes = bytes[pos:(pos + byte_length - 1)]
        str = StringEncodings.decode(extracted_bytes, encoding)
        push!(strings, str)
        pos += byte_length
    end
    # keep the element type String unless a string is missing
    any(ismissing, strings) || (strings = convert(Vector{String}, strings))

    if num_strings == 1
        return first(strings)
    else
        return reshape(strings, shape...)
    end
end
