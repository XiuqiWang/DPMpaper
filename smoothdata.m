function smoothed_data=smoothdata(Q)
data = Q;
% Define window size for moving average
window_size = 30;
% Apply moving average filter
smoothed_data = movmean(data, window_size);
% Define a finer window size for the beginning part
window_size_begin = 5;
smoothedBegin = movmean(data(1:ceil(0.3*end)),window_size_begin);
smoothed_data(1:ceil(0.3*end))=smoothedBegin;
end