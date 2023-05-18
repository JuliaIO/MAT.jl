module MAT_v4

import Base: read, write, close

mutable struct Matlabv4File
    ios::IOStream
    varnames::Dict{String, Int64}

    Matlabv4File(ios) = new(ios)
end

const V4_ELTYPE = [Float64, Float32, Int32, Int16, UInt16, UInt8]

matopen(ios::IOStream) = Matlabv4File(ios)

function unpack_header_type(type::Int32)
    T, P, O, M = digits(type; pad=4)
    iszero(O) || error("file is not a v4 MAT file (magic digit not 0)")
    iszero(M) || error("only little endian v4 MAT currently supported")
    iszero(T) || error("only full numeric matrices supported for v4 MAT files")

    P in 0:5 || error("invalid eltype digit in v4 MAT file: $P")
    eltype = V4_ELTYPE[P+1]

    return eltype
end

function read_one_mat!(mat::Matlabv4File)
    type, mrows, ncols, imagf, namlen = read!(mat.ios, Vector{Int32}(undef, 5))
    eltype = unpack_header_type(type)

    iszero(imagf) || error("no imaginary matrix support")

    name = String(read!(mat.ios, Vector{UInt8}(undef, namlen-1)))

    # one null byte before start of matrix data
    Base.read(mat.ios, UInt8)

    @show mrows, ncols
    value = read!(mat.ios, Matrix{eltype}(undef, mrows, ncols))

    return name => value
end

function read(mat::Matlabv4File)
    seekstart(mat.ios)

    results = Dict{String,Any}()
    while !eof(mat.ios)
        name, value = read_one_mat!(mat)
        results[name] = value
    end
    
    return results
end

close(mat::Matlabv4File) = close(mat.ios)

end # module
