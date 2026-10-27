# opaque_datetime.jl
# Conversion of the MATLAB `datetime` class

"""
    to_zoned(utc::DateTime, tz::String)

A zoned datetime from its UTC instant and MATLAB's name for the zone (an IANA name such
as "Europe/London", "UTC", or a fixed offset such as "+05:30"). Defined by the TimeZones
extension (Julia 1.9 and later, when TimeZones.jl is loaded), which returns a
`ZonedDateTime`; without it, zoned datetimes are read as their UTC instant.
"""
function to_zoned end

function from_datetime(obj::MatlabOpaque)
    dat = obj["data"]
    if isnothing(dat) || isempty(dat)
        return DateTime[]
    end
    #isdate = obj["isDateOnly"] # optional: convert to Date instead of DateTime?
    tz = haskey(obj, "tz") ? obj["tz"] : ""
    if !(tz isa AbstractString) || isempty(tz)
        return map_or_not(ms_to_datetime, dat)  # no zone: the stored wall-clock time
    end
    # MATLAB stores a zoned datetime as its UTC instant, except in "UTCLeapSeconds", where
    # the count includes leap seconds, which DateTime cannot represent
    if tz == "UTCLeapSeconds"
        @warn "datetime with timezone \"UTCLeapSeconds\" is not converted (leap seconds cannot be represented); returning the MatlabOpaque"
        return obj
    end
    utc = map_or_not(ms_to_datetime, dat)
    if !hasmethod(to_zoned, Tuple{DateTime,String})
        @warn "datetime with timezone \"$tz\" is returned as its UTC instant; load TimeZones.jl to read it as a ZonedDateTime"
        return utc
    end
    try
        return map_or_not(t -> ismissing(t) ? missing : to_zoned(t, String(tz)), utc)
    catch err
        @warn "timezone \"$tz\" was not recognised ($(sprint(showerror, err))); returning the UTC instant"
        return utc
    end
end

# The imaginary part is the low half of a double-double: smaller than half an ulp of the
# real part, so real + imag == real in Float64 and nothing is lost at DateTime resolution.
ms_to_datetime(ms::Complex) = ms_to_datetime(real(ms))
function ms_to_datetime(ms::Real)
    isfinite(ms) || return missing  # NaT is stored as NaN; DateTime has no ±Inf either
    s, ms_rem = fldmod(ms, 1_000)  # whole seconds and remainder milliseconds
    return DateTime(1970, 1, 1) + Second(s) + Millisecond(round(Int, ms_rem))
end

function to_matlab_data(d::DateTime)
    ms = Dates.value(d)

    # matlab counts w.r.t. 1970
    matlab_offset = Dates.value(DateTime(1970, 1, 1))

    # note: can be negative, that's fine
    matlab_ms = ms - matlab_offset
    return Complex(Float64(matlab_ms))
end

function MatlabOpaque(d::ScalarOrArray{DateTime})
    return MatlabOpaque(
        Dict{String,Any}(
            "tz"   => "",
            "data" => map_or_not(to_matlab_data, d),
            "fmt"  => "",
        ),
        "datetime"
    )
end
