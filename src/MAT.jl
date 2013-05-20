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

require("MAT/src/MAT_HDF5")
require("MAT/src/MAT_v5")
require("MAT/src/MAT_macros")
module MAT
using MAT_HDF5, MAT_v5
import Base.read, Base.write

export matopen, matread, matwrite, save, load

# Open a MATLAB file
function matopen(filename::String, rd::Bool, wr::Bool, cr::Bool, tr::Bool, ff::Bool)
    # When creating new files, create as HDF5 by default
    fs = filesize(filename)
    if cr && (tr || fs == 0)
        return MAT_HDF5.matopen(filename, rd, wr, cr, tr, ff)
    elseif fs == 0
        error("File \"$filename\" does not exist and create was not specified")
    end

    # Test whether this is a MAT file
    if fs < 19
        error("File \"$filename\" is too small to be a supported MAT file")
    end
    rawfid = open(filename, "r")
    magic = read(rawfid, Uint8, 19)
    close(rawfid)

    # Send to appropriate module
    magic == "MATLAB 7.3 MAT-file".data ? MAT_HDF5.matopen(filename, rd, wr, cr, tr, ff) :
    magic == "MATLAB 5.0 MAT-file".data ? MAT_v5.matopen(filename, rd, wr, cr, tr, ff) :
        error("\"$filename\" is not a MAT file, or is an unsupported (v4) MAT file")
end

function matopen(fname::String, mode::String)
    mode == "r"  ? matopen(fname, true , false, false, false, false) :
    mode == "r+" ? matopen(fname, true , true , false, false, false) :
    mode == "w"  ? matopen(fname, false, true , true , true , false) :
#     mode == "w+" ? matopen(fname, true , true , true , true , false) :
#     mode == "a"  ? matopen(fname, false, true , true , false, true ) :
#     mode == "a+" ? matopen(fname, true , true , true , false, true ) :
    error("invalid open mode: ", mode)
end

matopen(fname::String) = matopen(fname, "r")

# Read all variables from a MATLAB file
function matread(filename::String)
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
function matwrite{S, T}(filename::String, dict::Dict{S, T})
    file = matopen(filename, "w")
    try
        for (k, v) in dict
            local kstring
            try
                kstring = ascii(k)
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
