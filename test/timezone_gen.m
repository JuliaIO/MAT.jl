% Generates test/v7/timezone.mat and test/v7.3/timezone.mat: datetimes with and without
% time zones, and timetables whose row times have one.
% Run from the test directory: matlab -batch "timezone_gen"

dt_unzoned = datetime(2022, 7, 20, 12, 0, 0);
dt_utc = datetime(2022, 7, 20, 12, 0, 0, 'TimeZone', 'UTC');
dt_london_summer = datetime(2022, 7, 20, 12, 0, 0, 'TimeZone', 'Europe/London');   % BST, UTC+1
dt_london_winter = datetime(2022, 1, 20, 12, 0, 0, 'TimeZone', 'Europe/London');   % GMT
dt_offset = datetime(2022, 7, 20, 12, 0, 0, 'TimeZone', '+05:30');
dt_newyork = datetime(2022, [1 7], 20, 12, 0, 0, 'TimeZone', 'America/New_York'); % EST, EDT
dt_leap = datetime(2016, 12, 31, 23, 59, 60, 'TimeZone', 'UTCLeapSeconds');

Time = datetime(2022, 7, 20, 12, 0, 0, 'TimeZone', 'Europe/London') + seconds([0; 1.5]);
tt_zoned = timetable(Time, [1; 2], 'VariableNames', {'x'});
tt_zoned_step = timetable([1; 2; 3], 'TimeStep', seconds(0.5), ...
    'StartTime', datetime(2022, 7, 20, 12, 0, 0, 'TimeZone', 'Europe/London'), 'VariableNames', {'y'});

vars = {'dt_unzoned', 'dt_utc', 'dt_london_summer', 'dt_london_winter', 'dt_offset', ...
        'dt_newyork', 'dt_leap', 'tt_zoned', 'tt_zoned_step'};
save(fullfile('v7', 'timezone.mat'), vars{:}, '-v7');
save(fullfile('v7.3', 'timezone.mat'), vars{:}, '-v7.3');
