# MAT.jl
# Tools for reading MATLAB v5 files in Julia
#
# Copyright (C) 2012   Timothy E. Holy and Simon Kornblith
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

VERSION >= v"0.4.0-dev+6521" && __precompile__()

module MAT

using HDF5, Compat

include("MAT_HDF5.jl")
include("MAT_v5.jl")

using .MAT_HDF5, .MAT_v5

export matopen, matread, matwrite, names, exists, @read, @write

# Open a MATLAB file
const HDF5_HEADER = UInt8[0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]
function matopen(filename::AbstractString, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool)
    # When creating new files, create as HDF5 by default
    fs = filesize(filename)
    if cr && (tr || fs == 0)
        return MAT_HDF5.matopen(filename, rd, wr, cr, tr, ff)
    elseif fs == 0
        error("File \"$filename\" does not exist and create was not specified")
    end

    # Test whether this is a MAT file
    if fs < 128
        error("File \"$filename\" is too small to be a supported MAT file")
    end
    rawfid = open(filename, "r")

    # Check for MAT v4 file
    magic = read!(rawfid, Vector{UInt8}(4))
    for i = 1:length(magic)
        if magic[i] == 0
            close(rawfid)
            error("\"$filename\" is not a MAT file, or is an unsupported (v4) MAT file")
        end
    end

    # Check for MAT v5 file
    seek(rawfid, 124)
    version = read(rawfid, UInt16)
    endian_indicator = read(rawfid, UInt16)
    if (version == 0x0100 && endian_indicator == 0x4D49) ||
       (version == 0x0001 && endian_indicator == 0x494D)
        if wr || cr || tr || ff
            error("creating or appending to MATLAB v5 files is not supported")
        end
        return MAT_v5.matopen(rawfid, endian_indicator)
    end

    # Check for HDF5 file
    for offset = 512:512:fs-8
        seek(rawfid, offset)
        if read!(rawfid, Vector{UInt8}(8)) == HDF5_HEADER
            close(rawfid)
            return MAT_HDF5.matopen(filename, rd, wr, cr, tr, ff)
        end
    end

    close(rawfid)
    error("\"$filename\" is not a MAT file")
end

function matopen(fname::AbstractString, mode::AbstractString)
    mode == "r"  ? matopen(fname, true , false, false, false, false) :
    mode == "r+" ? matopen(fname, true , true , false, false, false) :
    mode == "w"  ? matopen(fname, false, true , true , true , false) :
#     mode == "w+" ? matopen(fname, true , true , true , true , false) :
#     mode == "a"  ? matopen(fname, false, true , true , false, true ) :
#     mode == "a+" ? matopen(fname, true , true , true , false, true ) :
    error("invalid open mode: ", mode)
end

matopen(fname::AbstractString) = matopen(fname, "r")

function matopen(f::Function, args...)
    fid = matopen(args...)
    try
        f(fid)
    finally
        close(fid)
    end
end

"""
    matopen(filename [, mode]) -> handle
    matopen(f::Function, filename [, mode]) -> f(handle)

Mode defaults to "r" for read.  It can also be "w" for write, or "r+" for
read or write without creation or truncation.

Use with `read`, `write`, `close`, `names`, and `exists`.
"""
matopen

# Read all variables from a MATLAB file
"""
    matread(filename) -> Dict

Return a dictionary of all the variables and values in a Matlab file,
opening and closing it automatically.
"""
function matread(filename::AbstractString)
    file = matopen(filename)
    local vars
    try
        vars = read(file)
    finally
        close(file)
    end
    vars
end

# Write a dict to a MATLAB file
"""
    matwrite(filename, d::Dict)

Write a dictionary containing variable names as keys and values as values
to a Matlab file, opening and closing it automatically.
"""
function matwrite{S, T}(filename::AbstractString, dict::Associative{S, T})
    file = matopen(filename, "w")
    try
        for (k, v) in dict
            local kstring
            try
                kstring = ascii(convert(String, k))
            catch x
                error("matwrite requires a Dict with ASCII keys")
            end
            write(file, kstring, v)
        end
    finally
        close(file)
    end
end
end
