function [M,MeanResI,MeanResF,SSi,SSf] = TwoDRigidBodyRemover(Matrix1,Matrix2,Offset) 
%Find XY rigid body of a Matrix using top slice.
    figure
    if Offset >= 0
        I1 = (Matrix1(:,:,0.8*length(Matrix1(1,1,:))));
        I2 = (Matrix2(:,:,0.8*length(Matrix2(1,1,:))+Offset));
    else
        I1 = (Matrix1(:,:,0.8*length(Matrix1(1,1,:))));
        I2 = (Matrix2(:,:,0.8*length(Matrix2(1,1,:))+Offset));       
    end
    I2f = registerImages(I2,I1);
    M= I2f.Transformation;
    IRi = abs(I2-I1);
    MeanResI = mean(IRi,"all");
    disp("Mean Residual prewarp" + MeanResI)
    
    tiledlayout(2,3)
    sgtitle('Images registration')
    nexttile
    imshow(I1)
    title('Image 1')
    nexttile
    imshow(I2)
    title('Image 2')
    nexttile
    imshow(IRi)
    title({ "\fontsize{8}Mean Residual prewarp:"+ MeanResI ;'\fontsize{15}Residual'})

    nexttile
    imshow(I1)
    title('Image 1')
    nexttile
    imshow(I2f.RegisteredImage)
    title('I2 warped using Registration Estimator')
    IRf = abs(I2f.RegisteredImage-I1);
    
    MeanResF = mean(IRf,"all");
    disp("Mean Residual postwarp" + MeanResF)
    nexttile
    imshow(IRf)
     title({ "\fontsize{8}Mean Residual postwarp:"+ MeanResF ;'\fontsize{15}Residual'})

   SSi = ssim(I2,I1);
   SSf = ssim(I2f.RegisteredImage,I1);
end