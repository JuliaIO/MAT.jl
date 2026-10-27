# Zoned datetimes with TimeZones.jl loaded: ZonedDateTime, through the package extension
# (Julia 1.9 and later). read.jl reads the same file before TimeZones.jl is loaded.
using TimeZones, Dates

@testset "time zones with TimeZones.jl" begin
    for format in ("v7", "v7.3")
        # timezone.mat is written in MATLAB by timezone_gen.m
        filepath = joinpath(dirname(@__FILE__), format, "timezone.mat")
        if !isfile(filepath) || VERSION < v"1.9"
            @test_skip isfile(filepath) && VERSION >= v"1.9"
            continue
        end
        vars = @test_logs (:warn, r"UTCLeapSeconds") match_mode=:any matread(filepath)
        london, newyork = tz"Europe/London", tz"America/New_York"

        @test vars["dt_london_summer"] isa ZonedDateTime
        @test vars["dt_london_summer"] == ZonedDateTime(2022, 7, 20, 12, london)   # BST
        @test timezone(vars["dt_london_summer"]) == london
        @test vars["dt_london_winter"] == ZonedDateTime(2022, 1, 20, 12, london)   # GMT
        @test vars["dt_utc"] == ZonedDateTime(2022, 7, 20, 12, tz"UTC")
        @test timezone(vars["dt_utc"]) == tz"UTC"
        @test timezone(vars["dt_offset"]) == FixedTimeZone("+05:30")               # a fixed offset
        @test DateTime(vars["dt_offset"]) == DateTime(2022, 7, 20, 12)
        @test vars["dt_newyork"] == [ZonedDateTime(2022, 1, 20, 12, newyork) ZonedDateTime(2022, 7, 20, 12, newyork)]
        @test all(timezone.(vars["dt_newyork"]) .== newyork)

        @test vars["dt_unzoned"] == DateTime(2022, 7, 20, 12)   # no zone: a DateTime, as before
        @test vars["dt_leap"] isa MatlabOpaque                  # leap seconds counted: left as read

        t = vars["tt_zoned"]                                    # zoned row times, one per row
        @test t[:Time] == ZonedDateTime(2022, 7, 20, 12, london) .+ Millisecond.([0, 1500])
        @test all(timezone.(t[:Time]) .== london)
        t = vars["tt_zoned_step"]                               # regular, from a zoned start
        @test t[:Time] == ZonedDateTime(2022, 7, 20, 12, london) .+ Millisecond.([0, 500, 1000])
        @test all(timezone.(t[:Time]) .== london)
    end
end
