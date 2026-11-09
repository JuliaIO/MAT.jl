% Generates test/v7/timetable.mat and test/v7.3/timetable.mat for the timetable tests
% in read.jl. Each variable is one way MATLAB stores a timetable's row times.
% Run from the test directory: matlab -batch "timetable_gen"

% explicit datetime row times, with sub-millisecond parts, and a string variable
Time = datetime(2022, 7, 20, 2, 27, 27) + milliseconds([428.25; 495.75; 563.5]);
Current = [1.5e-11; 2.5e-11; -3e-12];
Channel = ["AB"; "A"; "B"];
tt_datetime = timetable(Time, Current, Channel);

% explicit duration row times
tt_duration = timetable(seconds([0; 0.5; 1.25]), [1; 2; 3], 'VariableNames', {'x'});

% regular: a sample rate (row times start at 0 s)
x = [10; 20; 30; 40];
tt_rate = timetable(x, 'SampleRate', 1000);

% regular: a time step from a datetime start
y = [1; 2; 3];
tt_step = timetable(y, 'TimeStep', seconds(0.25), 'StartTime', datetime(2022, 7, 20, 2, 0, 0));

% no rows and no variables, as an acquisition leaves an unused timetable
tt_empty = timetable(datetime.empty(0, 1));

% a variable with two columns
tt_matrix = timetable(seconds([1; 2]), [1 2; 3 4], 'VariableNames', {'m'});

% the row-times dimension renamed
tt_dimname = tt_duration;
tt_dimname.Properties.DimensionNames{1} = 'Timestamp';

% regular in calendar months: no fixed sample rate, so not converted
z = [1; 2; 3];
tt_calendar = timetable(z, 'TimeStep', calmonths(1), 'StartTime', datetime(2022, 1, 1));

vars = {'tt_datetime', 'tt_duration', 'tt_rate', 'tt_step', 'tt_empty', 'tt_matrix', 'tt_dimname', 'tt_calendar'};
save(fullfile('v7', 'timetable.mat'), vars{:}, '-v7');
save(fullfile('v7.3', 'timetable.mat'), vars{:}, '-v7.3');
