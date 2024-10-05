function [MinEnergyPoint,MeanEnergy] = zshift(Matrix1,Matrix2,MovingAverage)
    %Calculate z shift in DCDC samples

%% Standard Deviation calculations
sd1 = zeros(length(Matrix1(1,1,:))-MovingAverage,1);
sd2 = zeros(length(Matrix1(1,1,:))-MovingAverage,1);
for z = MovingAverage+1:length(Matrix1(1,1,:))-MovingAverage
    sd1(z-MovingAverage) = std2(Matrix1(:,:,(z-MovingAverage):(z+MovingAverage)));
    sd2(z-MovingAverage) = std2(Matrix2(:,:,(z-MovingAverage):(z+MovingAverage)));
end
z = 1:length(sd1);

sd1 = sd1 - mean(sd1);
sd2 = sd2 - mean(sd2);

%% Plot std graph 1
figure
tiledlayout(3,1)
nexttile
hold on
plot(z,sd1)
plot(z,sd2)
hold off
title('\bfGraph Showing the Standard Deviation per Slice')
xlabel('\bfSlice')
ylabel('\bfStandard Error')
legend('Fixed Volume','Moving Volume')


%% Sample shift meathod
PositiveEnergy = nan(round(0.8*length(sd2),0),1);
NegativeEnergy = nan(round(0.8*length(sd2),0),1);
for i = 1:0.8*length(sd2)
    PositiveEnergy(i) = mean(abs(sd1(1:(length(sd2)-i)) - sd2(1+i:length(sd2))));
end
for i = 1:0.8*length(sd2)
    NegativeEnergy(i) = mean(abs(sd1(1+i:length(sd2)) - sd2(1:(length(sd1)-i))));
end
ZeroEnergy = mean(abs(sd1(1:length(sd2)) - sd2(1:(length(sd1)))));
% %% 1 Subset meathod
% PositiveEnergy = nan(round(0.8*length(sd2),0),1);
% NegativeEnergy = nan(round(0.8*length(sd2),0),1);
% for i = 1:round(0.8*length(sd2),0)
%     PositiveEnergy(i+length(sd2)) = mean(abs(sd1(i:(i+round(length(sd2)*0.2,0))) - sd2(1:(round(length(sd2)*0.2,0)+1))));
% end
% for i = 1:round(0.8*length(sd2),0)
%     NegativeEnergy(i) = mean(abs(sd2(i:(round(length(sd1)*0.2,0)+i)) - sd1(1:(round(length(sd1)*0.2,0)+1))));
% end
% ZeroEnergy = mean(abs(sd1(1:length(sd2)) - sd2(1:(length(sd1)))));
%% Find min energy

[MinPositiveEnergy,MinPositiveEnergyPoint] = min(PositiveEnergy);
[MinNegativeEnergy,MinNegativeEnergyPoint] = min(NegativeEnergy);

if MinNegativeEnergy < MinPositiveEnergy
    MinEnergy = MinNegativeEnergy;
    MinEnergyPoint = -MinNegativeEnergyPoint;
elseif MinNegativeEnergy > MinPositiveEnergy
    MinEnergy = MinPositiveEnergy;
    MinEnergyPoint = MinPositiveEnergyPoint;
end

if MinEnergy > ZeroEnergy
    print('Sample Registered')
end

Energy =[flip(NegativeEnergy);ZeroEnergy;PositiveEnergy];
MeanEnergy = mean(Energy);
%% Plot Energy Shift diagram
nexttile
hold on
plot(-length(NegativeEnergy):length(PositiveEnergy),Energy)
title('\bfGraph showing the change in Energy with zShift')
% title({'\rm\fontsize{8} 0<zShift<100' ;'\fontsize{15}\bfGraph showing the change in Energy with zShift'})
% xlim([0 100])
xlabel('\bfzShift(Pixels)')
ylabel({'\fontsize{8}Mean Difference in Standard Deviation';'\fontsize{15}\bf{Energy}'})
scatter(MinEnergyPoint,MinEnergy,250,'x','ColorVariable',[0.9290 0.6940 0.1250])
legend('Energy','Minimum Value','Location','southeast')
if MinEnergyPoint >= length(sd2)
    MinEnergyPoint= -(MinEnergyPoint-length(sd2));
end

%% Plot Corrected zShifts.
nexttile
hold on
plot(z+ones(1,length(z))*MinEnergyPoint,sd1)
plot(z,sd2)
hold off
title('\bfGraph Showing the Standard Deviation per Slice')
xlabel('\bfSlice')
ylabel('\bfStandard Error')
legend('Fixed Volume','Moving Volume')