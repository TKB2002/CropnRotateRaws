function [RotatedMatrix] = CropMatrix(Matrix,MinzShift,MaxzShift,zShift)

    for i =1-MinzShift:1:(length(Matrix(1,1,:))-MaxzShift+zShift)
        RotatedMatrix(:,:,i+MinzShift) = Matrix(:,:,i);
    end