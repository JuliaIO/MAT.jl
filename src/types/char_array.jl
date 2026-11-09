# char_array.jl
# Decoding of MATLAB char arrays to Julia strings

# MATLAB char array storage types:
# - miUINT8: ASCII or UTF-8. Can be directly converted to String.
# - miUINT16: one character code per value, as current MATLAB writes v6 text (codes above
#   U+00FF are lost on save). Older MATLAB versions packed UTF-8 bytes into the values,
#   two to a value; a row holding any value above 0x00FF is read that way, as before.
# - miUTF8/16/32: UTF-8/16/32 encoding. Transcode to UTF-8 if needed, then convert to String
_decode_row(row::AbstractVector{UInt8},  codec)  = String(row)
_decode_row(row::AbstractVector{UInt32}, codec) = String(transcode(UInt8, row))
function _decode_row(row::AbstractVector{UInt16}, codec)
    if codec == :utf16
        return String(transcode(UInt8, row))
    end

    # packed UTF-8 can leave a lone byte of a character in a value of its own
    packed = any(>(0x00ff), row)

    n = length(row)
    out = Vector{UInt8}(undef, 2n)
    j = 1

    @inbounds for x in row
        msb = UInt8(x >> 8)
        lsb = UInt8(x & 0xff)

        if msb == 0x00 && (lsb < 0x80 || packed)
            out[j] = lsb
            j += 1
        elseif msb == 0x00
            # the code point U+0080 to U+00FF, as two UTF-8 bytes
            out[j] = 0xc0 | (lsb >> 6)
            out[j+1] = 0x80 | (lsb & 0x3f)
            j += 2
        else
            out[j] = msb
            out[j+1] = lsb
            j += 2
        end
    end

    resize!(out, j - 1)
    return String(out)
end

function decode_char_array(arr::AbstractArray{T}, codec::Symbol) where T <: Union{UInt8, UInt16, UInt32}
    char_len = size(arr, 2)
    other_dims = (size(arr, 1), size(arr)[3:end]...)

    if char_len == 0
        # an empty char array of any shape, as the v7.3 reader returns it
        return String[]
    end

    if ndims(arr) == 2 && size(arr, 1) == 1
        # Return as single string instead of 1-element array
        row = view(arr, 1, :)
        return _decode_row(row, codec)
    end

    out = Array{String}(undef, other_dims...)
    for I in CartesianIndices(out)
        idx = Tuple(I)
        row = view(arr, idx[1], :, idx[2:end]...)
        out[I] = _decode_row(row, codec)
    end

    return out
end
