function [MinEnergyPoint,Error] = zshift(Matrix1,Matrix2,MovingAverage,Normaliser)
% zShift Calculates the shift in the z direction of samples by calculating and comparing slice standard
% devations 
% zshift(Matrix1,Matrix2,MovingAverage)
%
% Arguments:
% Matrix1 [Integer, nxnxn]: The Reference Matrix 
% Matrix2 [Integer, nxnxn]: The Moving Matrix
% MovingAverage[Integer]: The Moving average of slices in which to
% calculate the standard deviations from.
% Normaliser[Boolian]: To normalise standard deviations for comparison
%% Standard Deviation calculations

% Initiate
sd1 = zeros(length(Matrix1(1,1,:))-MovingAverage,1);
sd2 = zeros(length(Matrix1(1,1,:))-MovingAverage,1);

for z = MovingAverage+1:length(Matrix1(1,1,:))-MovingAverage % For every Slice
    sd1(z-MovingAverage) = std2(Matrix1(:,:,(z-MovingAverage):(z+MovingAverage))); % Calculate the s.d of the slice + Moving average in Reference Volume
    sd2(z-MovingAverage) = std2(Matrix2(:,:,(z-MovingAverage):(z+MovingAverage))); % Calculate the s.d of the slice + Moving average in Moving Volume
end

if Normaliser == true
    sd1 = sd1 - mean(sd1);
    sd2 = sd2 - mean(sd2);
end

%% Plot std graph 1
z = 1:length(sd1); % Initalise z (x-axis) for graphing

figure % Initialise figure
tiledlayout(3,1)

% Plot the 2 standard devilations overlapped for plot 1
nexttile
hold on
plot(z,sd1)
plot(z,sd2) 
hold off

title('\bfGraph Showing the Standard Deviation per Slice')
xlabel('\bfSlice')
ylabel('\bfStandard Error')
legend('Fixed Volume','Moving Volume')


%% Find shift - Sample shift meathod
%Initalise 
PositiveEnergy = nan(round(0.8*length(sd2),0),1);
NegativeEnergy = nan(round(0.8*length(sd2),0),1);

for i = 1:0.8*length(sd2)
    %Move Volume 2 up 1, take 1 off the end of Volume 1 to keep the size the same
    PositiveEnergy(i) = mean(abs(sd1(1:(length(sd2)-i)) - sd2(1+i:length(sd2)))); %Find total difference in standard devations for each shift
end
for i = 1:0.8*length(sd2)
    %Move Volume 1 up 1, take 1 off the end of Volume 2 to keep the size the same
    NegativeEnergy(i) = mean(abs(sd1(1+i:length(sd2)) - sd2(1:(length(sd1)-i)))); %Find total difference in standard devations for each shift
end
ZeroEnergy = mean(abs(sd1(1:length(sd2)) - sd2(1:(length(sd1))))); % Find standard devation difference in initial condition



%% Find min energy

% Find best match from both shift directions.
[MinPositiveEnergy,MinPositiveEnergyPoint] = min(PositiveEnergy);
[MinNegativeEnergy,MinNegativeEnergyPoint] = min(NegativeEnergy);

% Find which direction shift had best match and select the best match to 
if MinNegativeEnergy < MinPositiveEnergy 
    MinEnergy = MinNegativeEnergy;
    MinEnergyPoint = -MinNegativeEnergyPoint;
    ShiftPos = false;
elseif MinNegativeEnergy > MinPositiveEnergy
    MinEnergy = MinPositiveEnergy;
    MinEnergyPoint = MinPositiveEnergyPoint;
    ShiftPos = true;
end

% Make sure the sample isnt already registered as best as possible
if MinEnergy > ZeroEnergy
    print('Sample Registered')
    MinEnergyPoint = 0;
end

% Create varible with all Energy's for plotting
Energy =[flip(NegativeEnergy);ZeroEnergy;PositiveEnergy];

%% Create Variable to show zShift quality

% Take the difference in Shift Quality and devide by the worst shift for
% an error parameter
Error = ((max(Energy)-MinEnergy)/max(Energy));


%% Plot Energy Shift diagram
nexttile
hold on
plot(-length(NegativeEnergy):length(PositiveEnergy),Energy)
title('\bfGraph showing the change in Energy with zShift')
title({'\rm\fontsize{8} 0<zShift<100' ;'\fontsize{15}\bfGraph showing the change in Energy with zShift'})
xlim([0 100])
xlabel('\bfzShift(Pixels)')
ylabel({'\fontsize{8}Mean Difference in Standard Deviation';'\fontsize{15}\bf{Energy}'})
scatter(MinEnergyPoint,MinEnergy,250,'x','ColorVariable',[0.9290 0.6940 0.1250])
legend('Energy','Minimum Value','Location','southeast')



%% Plot Corrected zShifts.
nexttile
hold on
plot(z,sd1)
plot(z+ones(1,length(z))*MinEnergyPoint,sd2)
hold off
title('\bfGraph Showing the Standard Deviation per Slice')
xlabel('\bfSlice')
ylabel('\bfStandard Error')
legend('Fixed Volume','Moving Volume')