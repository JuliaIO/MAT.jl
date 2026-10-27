# opaque_duration.jl
# Conversion of the MATLAB `duration` class

function from_duration(obj::MatlabOpaque)
    dat = obj["millis"]
    #fmt = obj["fmt"] # TODO: format, e.g. 'd' to Day
    if isnothing(dat) || isempty(dat)
        return Millisecond[]
    end
    return map_or_not(ms_to_millisecond, dat)
end

# milliseconds are stored as Float64: rounded to whole ones, as datetimes are, and a
# missing (NaN) or infinite duration is missing
ms_to_millisecond(ms::Real) = isfinite(ms) ? Millisecond(round(Int, ms)) : missing

to_matlab_millis(d::Millisecond) = Float64(Dates.value(d))

function MatlabOpaque(d::ScalarOrArray{Millisecond})
    return MatlabOpaque(
        Dict{String,Any}(
            "millis" => map_or_not(to_matlab_millis, d),
        ),
        "duration"
    )
end
