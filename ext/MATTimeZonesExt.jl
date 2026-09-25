module MATTimeZonesExt

# Zoned MATLAB datetimes as ZonedDateTime. MATLAB stores a zoned datetime as its UTC
# instant plus the zone's name: an IANA name ("Europe/London"), "UTC", or a fixed offset
# such as "+05:30", all of which TimeZone() parses.

using MAT, Dates, TimeZones

MAT.MAT_types.to_zoned(utc::DateTime, tz::String) = ZonedDateTime(utc, TimeZone(tz); from_utc = true)

end
