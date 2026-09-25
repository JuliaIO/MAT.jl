% Generates test/v7/missing_times.mat and test/v7.3/missing_times.mat: datetime and
% duration arrays holding missing (NaT, NaN), infinite, or fractional-millisecond values.
% Run from the test directory: matlab -batch "missing_gen"

dt_nat = [datetime(2022, 7, 20, 12, 0, 0), NaT];
dt_nat_scalar = NaT;
dt_inf = [datetime(2022, 7, 20, 12, 0, 0), datetime(Inf, 'ConvertFrom', 'posixtime'), ...
          datetime(-Inf, 'ConvertFrom', 'posixtime')];
dt_nat_zoned = [datetime(2022, 7, 20, 12, 0, 0, 'TimeZone', 'Europe/London'), NaT('TimeZone', 'Europe/London')];
dur_nan = [seconds(1.5), seconds(NaN)];
dur_nan_scalar = seconds(NaN);
dur_frac = milliseconds([1.5; 2.25; 2.75]);

vars = {'dt_nat', 'dt_nat_scalar', 'dt_inf', 'dt_nat_zoned', 'dur_nan', 'dur_nan_scalar', 'dur_frac'};
save(fullfile('v7', 'missing_times.mat'), vars{:}, '-v7');
save(fullfile('v7.3', 'missing_times.mat'), vars{:}, '-v7.3');
