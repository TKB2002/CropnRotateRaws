function [RotatedMatrix] = RotateAndCropMatrix(Matrix,Rotation,MinzShift,MaxzShift,zShift,uint)
% RotateAndCropMatrix uses the inputted zShift's for whole loadsteps and
% rotation infomation to make series of equal size rotated Volumes.
%   [RotatedMatrix] = RotateAndCropMatrix(Matrix,Rotation,MinzShift,MaxzShift,zShift)
%
%   Arguments:
%       Matrix double {mustBeNonempty}: Matrix to be rotated
%       Rotation Geometric transformation object {mustBeNonempty}: The second matrix
%       Offset double: The offset value
%       uint single: The bit depth of number in '8' or '16'
%   Returns:
%   RotatedMatrix double: Rotated Matrix
%       
%    
RotatedMatrix = uint16(zeros(length(Matrix(:,1,1)),length(Matrix(1,:,1)),length(Matrix(1,1,:))-MaxzShift+MinzShift));
    if uint == 16
        RotatedMatrix = uint16(zeros(length(Matrix(:,1,1)),length(Matrix(1,:,1)),length(Matrix(1,1,:))-MaxzShift+MinzShift));
        for i =1-MinzShift+zShift:1:(length(Matrix(1,1,:))-MaxzShift+zShift)            
            RotatedMatrix(:,:,i+MinzShift-zShift) = uint16(imwarp(Matrix(:,:,i), imref2d(size(Matrix(:,:,i))), Rotation, 'OutputView', imref2d(size(Matrix(:,:,i))), 'SmoothEdges', true));
        end
    elseif uint == 8
        RotatedMatrix = uint8(zeros(length(Matrix(:,1,1)),length(Matrix(1,:,1)),length(Matrix(1,1,:))-MaxzShift+MinzShift));
        for i =1-MinzShift+zShift:1:(length(Matrix(1,1,:))-MaxzShift+zShift)
            RotatedMatrix(:,:,i+MinzShift-zShift) = uint8(imwarp(Matrix(:,:,i), imref2d(size(Matrix(:,:,i))), Rotation, 'OutputView', imref2d(size(Matrix(:,:,i))), 'SmoothEdges', true));
        end

    end