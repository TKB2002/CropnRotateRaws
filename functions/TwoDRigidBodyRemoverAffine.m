function [M,MeanResI,MenResF,SSi,SSf] = TwoDRigidBodyRemoverAffine(Matrix1,Matrix2,Offset)
% TwoDRigidBodyRemoverAffine  Find XY rigid body of a Matrix using top slice.
%   [M,MeanResI,MenResF,SSi,SSf] = TwoDRigidBodyRemoverAffine(Matrix1,Matrix2,Offset) adds A to itself.
%
%   Arguments:
%       Matrix1 double {mustBeNonempty}: The first matrix
%       Matrix2 double {mustBeNonempty}: The second matrix
%       Offset double: The offset value
%   Returns:
%       
%       
arguments
    Matrix1 double {mustBeNonempty}
    Matrix2 double {mustBeNonempty} 
    Offset double
end 

figure
if Offset >= 0
    I1 = (Matrix1(:,:,0.8*length(Matrix1(1,1,:))-Offset));
    I2 = (Matrix2(:,:,0.8*length(Matrix2(1,1,:))));
else
    I1 = (Matrix1(:,:,0.8*length(Matrix1(1,1,:))));
    I2 = (Matrix2(:,:,0.8*length(Matrix2(1,1,:))+Offset));
end
I2f = registerImagesAffine(I2,I1);
M= I2f.Transformation;
IRi = abs(I2-I1);
MeanResI = mean(IRi,"all");
disp("Mean Residual prewarp" + MeanResI)

tiledlayout(2,3)
sgtitle('Images registration')
nexttile
imshow(I1)
title('I1')
nexttile
imshow(I2)
title('I2')
nexttile
imshow(IRi)
title('Residual')


nexttile
imshow(I1)
title('I1')
nexttile
imshow(I2f.RegisteredImage)
title('I2 warped using Registration Estimator')
IRf = abs(I2f.RegisteredImage-I1);
title('Residual post warp')
MenResF = mean(IRf,"all");
disp("Mean Residual postwarp" + MenResF)
nexttile
imshow(IRf)

SSi = ssim(I2,I1);
SSf = ssim(I2f.RegisteredImage,I1);
end