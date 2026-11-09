# table.jl
# MatlabTable: conversion of the MATLAB `table` and `timetable` classes

"""
    MatlabTable(names, columns)

Representation of a MATLAB table, with column names in `names` and the corresponding
column arrays in `columns`. Implements the Tables.jl column interface.

For vector columns, convert to a DataFrames.jl table with `DataFrame(table)`; the
destination controls whether columns are copied. To convert MATLAB tables while
reading a file, pass `table=DataFrame` to `matread`.

MATLAB columns with several values per row are retained as matrices. Adapt those
columns before converting to a destination that only accepts vectors.
"""
struct MatlabTable
    names::Vector{Symbol}
    columns::Vector{Any}  # concrete, so code using it is type stable; columns differ in type
end
Tables.istable(::Type{MatlabTable}) = true
Tables.columnaccess(::Type{MatlabTable}) = true
Tables.columns(t::MatlabTable) = t
Tables.columnnames(t::MatlabTable) = t.names
Tables.getcolumn(t::MatlabTable, i::Int) = t.columns[i]
Tables.getcolumn(t::MatlabTable, nm::Symbol) = t[nm]
function find_index(m::MatlabTable, s::Symbol)
    idx = findfirst(isequal(s), m.names)
    if isnothing(idx)
        throw(KeyError(s))
    end
    return idx
end
function Base.getindex(m::MatlabTable, s::Symbol)
    idx = find_index(m, s)
    return getindex(m.columns, idx)
end
Base.getindex(m::MatlabTable, s::AbstractString) = getindex(m, Symbol(s))
MatlabTable(t::Tables.CopiedColumns{MatlabTable}) = Tables.source(t)

function from_table(obj::MatlabOpaque, ::Type{T}=MatlabTable) where {T}
    names = vec(Symbol.(obj["varnames"]))
    cols = vec([try_vec(c) for c in obj["data"]])
    t = MatlabTable(names, cols)
    return T(Tables.CopiedColumns(t))
end
# option to not convert and get the MatlabOpaque as table
from_table(obj::MatlabOpaque, ::Type{Nothing}) = obj

# A timetable is a table whose first column holds the row times, named after the row
# dimension ("Time" unless renamed), then the variables. Row times are DateTime or
# Millisecond; they are stored one per row, or, for a regular timetable, as a start time
# and a sample rate (or time step), from which they are generated here. Row times that
# cannot be interpreted leave the timetable as the MatlabOpaque it was read as.
function from_timetable(obj::MatlabOpaque, ::Type{T}=MatlabTable) where {T}
    tt = haskey(obj, "any") ? obj["any"] : obj
    times = timetable_rowtimes(tt["rowTimes"], Int(tt["numRows"]))
    if isnothing(times)
        @warn "timetable row times of this kind are not converted; returning the MatlabOpaque"
        return obj
    end
    names = Symbol[Symbol(tt["dimNames"][1]); Symbol.(vec(tt["varNames"]))]
    cols = vcat(Any[times], Any[try_vec(c) for c in vec(tt["data"])])
    t = MatlabTable(names, cols)
    return T(Tables.CopiedColumns(t))
end
from_timetable(obj::MatlabOpaque, ::Type{Nothing}) = obj

# row times stored one per row, already converted to DateTime or Millisecond
timetable_rowtimes(times::AbstractArray, n::Int) = vec(times)
# a regular timetable: a start time and a sample rate in Hz; MATLAB also stores the rate
# when the time step was given instead, so the rate serves for both
function timetable_rowtimes(regular::AbstractDict, n::Int)
    origin, rate = get(regular, "origin", nothing), get(regular, "sampleRate", nothing)
    # the start is a DateTime, a ZonedDateTime (TimeZones extension) or a duration
    if !(origin isa Union{Dates.AbstractDateTime,Dates.Period} && rate isa Real && isfinite(rate) && rate > 0)
        return nothing
    end
    return [origin + Millisecond(round(Int, 1000 * k / rate)) for k in 0:(n-1)]
end
timetable_rowtimes(times, n::Int) = nothing

try_vec(c::Vector) = c
try_vec(c) = [c]
function try_vec(c::AbstractArray)
    return (size(c, 2) == 1) ? vec(c) : c
end
