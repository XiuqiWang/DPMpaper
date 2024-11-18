function [peak_value, t_peak, index_peak, t_peak_error]=getIndexPeak(t,Value)
Value_smooth=Value;
[max_value, max_index] = max(Value_smooth);
t_peak=t(max_index);
peak_value = max_value;
index_peak = max_index;

% Fit a Gaussian to the data near the peak
Valuestd = std(Value_smooth);
lower_bounds = [peak_value*0.8, t_peak*0.8, Valuestd*0.8];  % replace with actual lower bounds
upper_bounds = [peak_value*1.2, t_peak*1.2, Valuestd*1.2];  % replace with actual upper bounds

% Perform the fit with bounds
fit_obj = fit(t', Value_smooth', 'gauss1', 'Lower', lower_bounds, 'Upper', upper_bounds);

% Get confidence intervals for the Gaussian fit parameters
ci = confint(fit_obj);            % 95% confidence intervals for all parameters

% The peak of the Gaussian is stored in the 'b1' parameter of the fit object
t_peak_fit = fit_obj.b1;              % x-value where the peak occurs (center of the Gaussian)

% The error on the peak is half the width of the confidence interval for 'b1'
t_peak_error = (ci(2,2) - ci(1,2)) / 2;
end