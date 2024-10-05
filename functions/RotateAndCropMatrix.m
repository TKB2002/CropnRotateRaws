function [RotatedMatrix] = RotateAndCropMatrix(Matrix,Rotation,MinzShift,MaxzShift,zShift)

    for i =1-MinzShift:1:(length(Matrix(1,1,:))-MaxzShift+zShift)
        RotatedMatrix(:,:,i+MinzShift) = imwarp(Matrix(:,:,i), imref2d(size(Matrix(:,:,i))), Rotation, 'OutputView', imref2d(size(Matrix(:,:,i))), 'SmoothEdges', true);
    end