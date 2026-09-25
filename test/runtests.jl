using SparseArrays, LinearAlgebra
using Test, MAT

@testset "MAT" begin
    include("types.jl")
    include("read.jl")
    include("readwrite4.jl")
    include("write.jl")
    include("timezones.jl")   # last: loading TimeZones.jl changes how zoned datetimes read
end
