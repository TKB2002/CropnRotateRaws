function printResidualImages(Vol1,Vol2,VolRemapped,VolResidual,zSlice)

nexttile
imshow(Vol1(:,:,zSlice))
title('Image1')
nexttile
imshow(Vol2(:,:,zSlice))
title('Image2')
nexttile
imshow(VolRemapped(:,:,zSlice))
title('Image 2 Remapped')
nexttile
imshow(VolResidual(:,:,zSlice))
title('Residual')