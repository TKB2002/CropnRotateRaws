function checkVersion()
%CHECKVERSION Summary of this function goes here
%   Detailed explanation goes here
    if ~ispc == true
        warning('This code was written on a Windows PC. Differences may occur.')
    end
    if isMATLABReleaseOlderThan('R2019a') == true
        warning('This code was written using MATLAB R2019.')
    end
end

