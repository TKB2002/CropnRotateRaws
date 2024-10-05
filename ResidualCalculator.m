function Res = ResidualCalculator(Matrix1,Matrix2)

for i = 1:min(length(Matrix2(1,1,:)),length(Matrix1(1,1,:)))
    Res(:,:,i)= Matrix2(:,:,i) - Matrix1(:,:,i);
end