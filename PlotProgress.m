function PlotProgress(iteration,total,prefix,length_value)
%PlotProgress  Prints in the command window the progression of a for loop.
%   PlotProgress(iteration,total,prefix,length_value)
%
%Arguments:
%   iteration [integer]: The iteration value of a for loop i.e. "i"
%   total [integer]: The total value i.e. for i = 1:total
%   prefix [string]: The prefix string before the progress bar
%   length_value: The number of characters wide for the progress bar to be
%
% Written by R J Scales 2024-17-09 in MATLAB 2024a
% Taken inspiration from a Python equivalent from the following:
% https://stackoverflow.com/a/34325723
arguments
    iteration (1,1) double {mustBeNonempty,mustBeInteger,mustBeFinite}
    total (1,1) double {mustBeNonempty,mustBeInteger,mustBeFinite} 
    prefix (1,1) string 
    length_value (1,1) double {mustBeNonempty,mustBeInteger,mustBeFinite} 
end

if and(iteration == 1, isMATLABReleaseOlderThan('R2024a') == true) == true
    warning('PlotProgress was written using MATLAB R2024a...')
end

fill = char(9724);
percent = sprintf('%0.1f', 100 * (iteration / total));
filledLength = floor((length_value * iteration) / total);
bar = horzcat(repelem(fill, 1, filledLength), repelem('-', 1, length_value - filledLength));
fprintf("\r%s |%s| {%s}\r", prefix, bar, percent)
% Print New Line on Complete
if iteration == total
    disp('')

end

