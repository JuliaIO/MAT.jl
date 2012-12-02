# MAT.jl
# Tools for reading MATLAB v5 files in Julia

# Copyright (C) 2012   Simon Kornblith

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
load("UTF16")

module MAT
using Base, UTF16

export matread

const miINT8 = 1
const miUINT8 = 2
const miINT16 = 3
const miUINT16 = 4
const miINT32 = 5
const miUINT32 = 6
const miSINGLE = 7
const miDOUBLE = 9
const miINT64 = 12
const miUINT64 = 13
const miMATRIX = 14
const miCOMPRESSED = 15
const miUTF8 = 16
const miUTF16 = 17
const miUTF32 = 18
const mxCELL_CLASS = 1
const mxSTRUCT_CLASS = 2
const mxOBJECT_CLASS = 3
const mxCHAR_CLASS = 4
const mxSPARSE_CLASS = 5
const mxDOUBLE_CLASS = 6
const mxSINGLE_CLASS = 7
const mxINT8_CLASS = 8
const mxUINT8_CLASS = 9
const mxINT16_CLASS = 10
const mxUINT16_CLASS = 11
const mxINT32_CLASS = 12
const mxUINT32_CLASS = 13
const mxINT64_CLASS = 14
const mxUINT64_CLASS = 15

const READ_TYPES = Type[Int8, Uint8, Int16, Uint16, Int32, Uint32, Float32, None, Float64,
	None, None, Int64, Uint64]
const CONVERT_TYPES = Type[None, None, None, String, None, Float64, Float32, Int8, Uint8,
	Int16, Uint16, Int32, Uint32, Int64, Uint64]

read_bswap{T}(f::IOStream, swap_bytes::Bool, ::Type{T}) = 
	swap_bytes ? bswap(read(f, T)) : read(f, T)
read_bswap{T}(f::IOStream, swap_bytes::Bool, ::Type{T}, dim::Union(Int, (Int...))) = 
	swap_bytes ? [bswap(x) for x=read(f, T, dim)] : read(f, T, dim)

skip_padding(f::IOStream, nbytes::Int64, hbytes::Int) = if nbytes % hbytes != 0
	skip(f, hbytes-(nbytes % hbytes))
end

function read_header(f::IOStream, swap_bytes::Bool)
	dtype = read_bswap(f, swap_bytes, Uint32)

	if (dtype & 0xFFFF0000) != 0
		# Small Data Element Format
		(dtype & 0x0000FFFF, int64(dtype >> 16), 4)
	else
		# Data Element Format
		(dtype, int64(read_bswap(f, swap_bytes, Uint32)), 8)
	end
end

function read_element{T}(f::IOStream, swap_bytes::Bool, ::Type{T})
	(dtype, nbytes, hbytes) = read_header(f, swap_bytes)
	data = read_bswap(f, swap_bytes, T, div(nbytes, sizeof(T)))
	skip_padding(f, nbytes, hbytes)
	data
end

function read_data{T}(f::IOStream, swap_bytes::Bool, ::Type{T}, dimensions::Vector{Int32})
	(dtype, nbytes, hbytes) = read_header(f, swap_bytes)
	read_type = READ_TYPES[dtype]

	read_array = any(dimensions .!= 1)
	if sizeof(read_type)*prod(dimensions) != nbytes
		error("Invalid element length")
	end
	if read_array
		data = read_bswap(f, swap_bytes, read_type, tuple(int(dimensions)...))
	else
		data = read_bswap(f, swap_bytes, read_type)
	end
	skip_padding(f, nbytes, hbytes)

	if read_type != T
		data = read_array ? convert(Array{T}, data) : convert(T, data)
	end
	data
end

function read_cell(f::IOStream, swap_bytes::Bool, dimensions::Vector{Int32})
	data = cell(int(dimensions)...)
	for i = 1:numel(data)
		(ignored_name, data[i]) = read_matrix(f, swap_bytes)
	end
	data
end

function read_struct(f::IOStream, swap_bytes::Bool, dimensions::Vector{Int32}, is_object::Bool)
	field_length = read_element(f, swap_bytes, Int32)[1]
	field_names = read_element(f, swap_bytes, Uint8)
	n_fields = div(length(field_names), field_length)
	data = Dict{String, Any}(n_fields)

	if is_object
		data["class"] = ascii(read_element(f, swap_bytes, Uint8))
	end

	for i = 1:n_fields
		sname = field_names[(i-1)*field_length+1:i*field_length]
		index = findfirst(sname, 0)
		sname = ascii(index == 0 ? sname : sname[1:index-1])
		(ignored_name, data[sname]) = read_matrix(f, swap_bytes)
	end
	data
end

function read_string(f::IOStream, swap_bytes::Bool, dimensions::Vector{Int32})
	(dtype, nbytes, hbytes) = read_header(f, swap_bytes)
	if dtype <= 2 || dtype == 16
		if dimensions[1] == 1
			data = utf8(read(f, Uint8, dimensions[2]))
		else
			data = Array(UTF16String, dimensions[1])
			for i = 1:dimensions[1]
				data[i] = utf8(read(f, Uint8, dimensions[2]))
			end
		end
	elseif dtype <= 4 || dtype == 17
		if dimensions[1] == 1
			data = UTF16String(read(f, Uint16, dimensions[2]))
		else
			data = Array(UTF16String, dimensions[1])
			for i = 1:dimensions[1]
				data[i] = UTF16String(read(f, Uint16, dimensions[2]))
			end
		end
	else
		error("Unsupported string type")
	end
	skip_padding(f, nbytes, hbytes)
	data
end

function read_matrix(f::IOStream, swap_bytes::Bool)
	(dtype, nbytes) = read_header(f, swap_bytes)
	if dtype != miMATRIX
		error("Unexpected data type")
	end

	flags = read_element(f, swap_bytes, Uint32)
	dimensions = read_element(f, swap_bytes, Int32)
	name = ascii(read_element(f, swap_bytes, Uint8))

	class = flags[1] & 0xFF
	local data
	if class == mxCELL_CLASS
		data = read_cell(f, swap_bytes, dimensions)
	elseif class == mxSTRUCT_CLASS || class == mxOBJECT_CLASS
		data = read_struct(f, swap_bytes, dimensions, class == mxOBJECT_CLASS)
	elseif class == mxSPARSE_CLASS
		error("Sparse matrices not currently supported")
	elseif class == mxCHAR_CLASS && (length(dimensions) <= 2 || all(dimensions[3:end] .== 2))
		data = read_string(f, swap_bytes, dimensions)
	else
		convert_type = CONVERT_TYPES[class]
		data = read_data(f, swap_bytes, convert_type, dimensions)
		if (flags[1] & 0x0A00) != 0
			data += im*read_data(f, swap_bytes, convert_type, dimensions)
		end
	end

	return (name, data)
end

function matread(file::ASCIIString)
	f = open(file, "r")
	header = read(f, Uint8, 116)
	skip(f, 8)
	version = read(f, Uint16)

	endian_indicator = read(f, Uint16)
	local swap_bytes
	if endian_indicator == 0x4D49
		swap_bytes = false
	elseif endian_indicator == 0x494D
		swap_bytes = true
	else
		error("Invalid endian indicator")
	end

	if swap_bytes
		version = bswap(version)
	end
	if version != 0x0100
		error("Unsupported MATLAB file version")
	end

	seek_end(f)
	eof = position(f)
	seek(f, 128)

	vars = Dict{ASCIIString, Any}()
	while position(f) != eof
		(name, data) = read_matrix(f, swap_bytes)
		vars[name] = data
	end
	
	close(f)
	vars
end
end