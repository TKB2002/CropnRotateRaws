%% Crop and rotate samples.
% Quick Automatic registration of sample's saved as .raw or equivalent
% Requires a boundary/large feature for z registration
% Does not work on samples with xz or xy rotations
% Word Documentation provided with MATLAB doc

% Code can be altered for different sample geometries

% By Tushar Bhudia (March 2024)
% tushar.bhudia@materials.ox.ac.uk
% tusharbhudia@gmail.com
% Version 1.0

% Many Thanks to James Marrow for advice, Thomas Zillhart's dataset for the
% idea, Marcus Williamson for the idea to use MATLAB rather then ImageJ and
% Robin Scales for help professionalising the code

%% Checking Computer
if ~ispc == true
    warning('This code was written on a Windows PC. Differences may occur.')
end
if isMATLABReleaseOlderThan('R2019a') == true
    warning('This code was written using MATLAB R2019.')
end

%% Import Functions
% This adds all the functions required automatically if they are in the 
% same folder as done in the word Doc
% Code taken from:
% https://uk.mathworks.com/matlabcentral/answers/247180-how-can-i-add-all-
% subfolders-to-my-matlab-path#answer_194998

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


%%  USER INPUTS

% Step 2
% Import the dimentions of your .raw file.
xdim =1000; % Units in Pixels
ydim =1000; % Units in Pixels
zdim =1000; % Units in Pixels

i8 = false;
i16 = true;
% Import the Amount of tomographs per loadstep, (Version 1.0 only works
% with 1 Scan per LoadStep.
ScansPerLoadStep = 1; % Amount of tomographs imaging the crack, per loadstep

% Import the zslice where your sample begins, everything else outside of
% this will be cropped.
Volume1StartOfSample = 0; %Z slice where the sample starts

% Step 2.1
% Does your .raw include filler around the material.
Cropping = false; % (false = No Cropping/ no filler, true= cropping/ has filler)

% Step 2.1.1
% If your sample needs cropping specify the pixels in it needs cropping
Croppedxdim = 800:1600; % Units in Pixels
Croppedydim = 800:1600; % Units in Pixels
Croppedzdim = 800:1600; % Units in Pixels

% Croppedxdim = 300:1700; % Units in Pixels
% Croppedydim = 300:1700; % Units in Pixels

% Step 2.2
% Does your sample include a boarder
Boarder = false; % (true = has a boarder, false = no boarder);


% Step 2.3
% Is their zShift in your sample?
Shift = true; % (true = their is a shift, false = no shift)

zShiftBinning = 1;

% Step 2.4
% Does your sample require xy rotation.
Rotation = true;

% Step 2.4.1
% Specify the rotation meathod
Method = false; % (true = affine warping, false = rigid warping

% Step 2.5
% Specify if your sample is flipped
Flipped = true;

% Step 3
% Calculate Residual?
ResidualQ = true;

%% Select files to import
% Can be replaced with file directories for further automation

% Select files to import
clear Files Path
disp('Select .raw files')
[Files,Path] = uigetfile({'*.raw';'*.raww'},'Load in All .raw files','MultiSelect','on'); % Displays UI for user to import files


%% Select File to Save data in

% Select folder to save Registered samples into
SavePath = uigetdir(Path,'Save registered .raw files'); % Displays UI for user to select save folder

while strcmp(Path,SavePath) == true
    warning('You must chose a different path to not overwrite the data!') % Displays warning if path chosen = source path
    SavePath = uigetdir(Path,'Save registered .raw files'); % Makes user change path until, path is different from source
end


%% Loading in the data and formatting it

% Clearing necessary variable names
clear OpenedFiles Volume Matrix

% Initalising

LoadStepsCount =length(Files)/ScansPerLoadStep; % Initalise simplified variable representing the amount of Load Steps


OpenedFiles = cell(ScansPerLoadStep,LoadStepsCount); % Initalise a cell containing file infomation
Volume = cell(ScansPerLoadStep,LoadStepsCount); % Initalise a cell containing the .raw data in vector form
Matrix = cell(ScansPerLoadStep,LoadStepsCount); % Initalise a cell containing the .raw data in matrix form

fprintf('Importing data');
tic_start = tic; %Start Timer

for j = 1:LoadStepsCount
    for i = 1:ScansPerLoadStep

        OpenedFiles{i,j} = fopen(strcat(Path,Files{(j-1)*ScansPerLoadStep+i})); % Get the file infomation
        if i16 == true
            Volume{i,j} = fread(OpenedFiles{i,j},xdim*ydim*zdim,'uint16'); % Read the .raw data into RAM
            Matrix{i,j} = uint16(reshape(Volume{i,j},xdim,ydim,[])); % Format the .raw data into Matrix form

        elseif i8 == true
            Volume{i,j} = fread(OpenedFiles{i,j},xdim*ydim*zdim,'uint8'); % Read the .raw data into RAM
            Matrix{i,j} = uint8(reshape(Volume{i,j},xdim,ydim,[])); % Format the .raw data into Matrix form
        end
        %Cropping Sample
        if Cropping == true %If the user specified to crop the sample in step 2.1
            Matrix{i,j} = Matrix{i,j}(Croppedxdim,Croppedydim,Croppedzdim); %Replace the Matrix data with data of a specific part.
        end
    end
end

clear Volume % Clear the redundant .raw data in vector form.

% This step is currently using double the ram it should be because it loads
% every loadstep twice when it should only load 1 loadstep at a time twice.
% Work is going into optimisng this

ImportTime = toc(tic_start); % Save the time required for the code to import and format data
fprintf('Data Imported %.2f\n', ImportTime); % Display time taken for data being imported

%% Removing sample boarder
if Boarder == true
    for i = 1:LoadStepsCount
        RemovedBoarderGreyLevelValue = mode(Matrix{1,i});
        Matrix{1,i}(Matrix{1,i}==RemovedBoarderGreyLevelValue) = NaN;
    end
end
%% Calculating the zShift of each sample
% We calculate sd of 1st Matrix every time we do this, we can optimise
% 2 meathods WIP

% Clearing necessary variable names
clear zShift
clear Offset

ProcessTimes = nan(LoadStepsCount,1);

if Shift == true
    % Initalising
    zShift = nan(LoadStepsCount,1); % Initalise a 1D array of zShifts relative to the first sample
    zShiftError = nan(LoadStepsCount,1);

    zShift(1) = 0;
    zShiftError(1) = 0;
    if Flipped == true
        zFlipShift = nan(LoadStepsCount,1); % Initalise a 1D array of zShifts relative to the first sample
        zFlipShiftError = nan(LoadStepsCount,1);

        zFlipShift(1) = 0;
        zFlipShiftError(1) = 0;
        zFlippedMatrix = cell(ScansPerLoadStep,LoadStepsCount);
    end

    for i = 2:LoadStepsCount
        tic
        fprintf('z Shift of Sample %d \n',i)
        [zShift(i),zShiftError(i)] = zshift(Matrix{1,1},Matrix{1,i},zShiftBinning); % Calculate zShift with novel function, standard deviation energy reduction meathod

        if Flipped == true
            zFlippedMatrix{i} = flip(Matrix{i},3);
            tic
            fprintf('z Shift of zFlipped Sample %d \n',i)
            [zFlipShift(i),zFlipShiftError(i)] = zshift(Matrix{1,1},zFlippedMatrix{1,i},zShiftBinning); % Calculate zShift with novel function, standard deviation energy reduction meathod
        end
        ProcessTimes(i) = toc; % Calculate time needed to run
        PlotProgress(i-1,LoadStepsCount-1,'zShift Progress',50)
    end
else
    zShift = zeros(LoadStepsCount);
end

%% Calculating xy transformations of registered z slices. Affine Method (User Input, Method = true)
if and(Rotation,Method) == true
    % Clearing necessary variable names
    clear RB Affine11 Affine12 Affine21 Affine22 Affine31 Affine32 RB MeanResI MeanResF

    %Initalising
    FailedXY = NaN(LoadStepsCount,1);
    RB = cell(LoadStepsCount,1);
    MeanResI = nan(LoadStepsCount,1);
    MeanResF = nan(LoadStepsCount,1);
    Affine11 = nan(LoadStepsCount,1);
    Affine12 = nan(LoadStepsCount,1);
    Affine21 = nan(LoadStepsCount,1);
    Affine22 = nan(LoadStepsCount,1);
    Affine31 = nan(LoadStepsCount,1);
    Affine32 = nan(LoadStepsCount,1);
    RotationPassed = nan(LoadStepsCount,1);
    SSIMi = nan(LoadStepsCount,1);
    SSIMf = nan(LoadStepsCount,1);

    % Setting cos terms of affine matrix on the reference matrix to 1 as cos(0) = 1
    Affine11(1) = 1;
    Affine22(1) = 1;
    disp(i) % Displaying the xy transformation calculating
    for i = 2:LoadStepsCount
        try
            [RB{i},MeanResI(i),MeanResF(i),SSIMi(i),SSIMf(i)] = TwoDRigidBodyRemoverAffine(Matrix{ScansPerLoadStep,1},Matrix{ScansPerLoadStep,i},zShift(i));
            Affine11(i) = RB{i}.T(1,1);
            Affine12(i) = RB{i}.T(1,2);
            Affine21(i) = RB{i}.T(2,1);
            Affine22(i) = RB{i}.T(2,2);
            Affine31(i) = RB{i}.T(3,1);
            Affine32(i) = RB{i}.T(3,2);
            RotationPassed(i) = true;
        catch
            FailedXY(i) = true;
            fprintf('xy rotation failed on sample %i\n', i )
            RotationPassed(i) = false;
            continue
        end
    end
end

%% Calculating xy transformations of registered z slices. Rigid Method (User Input, Method = false)
if and(Rotation,~Method) == true
    % Clearing necessary variable names
    clear RB RotationAngle xShift yShift FailedXY

    % Initalising

    RotationAngle= nan(LoadStepsCount,1);
    xShift = nan(LoadStepsCount,1);
    yShift= nan(LoadStepsCount,1);
    FailedXY = nan(LoadStepsCount,1);
    SSIMi = nan(LoadStepsCount,1);
    SSIMf = nan(LoadStepsCount,1);
    RotationPassed = nan(LoadStepsCount,1);

    RB = cell(LoadStepsCount,1);


    disp(i) % RJS: Why?
    for i = 2:length(Files)/ScansPerLoadStep
        try
            [RB{i},MeanResI(i),MeanResF(i),SSIMi(i),SSIMf(i)] = TwoDRigidBodyRemoverRigid(Matrix{ScansPerLoadStep,1},Matrix{ScansPerLoadStep,i},zShift(i));
            RotationAngle(i) = acos(RB{i}.T(1,1))*57.2958;
            xShift(i) = RB{i}.T(3,1);
            yShift(i) = RB{i}.T(3,2);
            RotationPassed(i) = true;
        catch
            FailedXY(i) = true;
            fprintf('xy rotation detection failed on sample %i\n', i )
            RotationPassed(i) = false;
            continue
        end
    end
end
%% ???TUSHAR TITLE OF THIS SECTION???
if Rotation == true
    figure('Name', 'Essential Shifts of Sample Per Load Step')
    yyaxis left
    hold on
    plot(RotationAngle)
    xlabel('LoadStep [1]')
    ylabel('Angle [\circ]')
    yyaxis right
    plot(xShift)
    plot(yShift)

    plot(zShift)
    legend('Rotation Angle','xshift','yshift','zshift')
end


%% Rotate Matrices

% Fix issues if RB doesnt work
clear MatrixRotated
MatrixRotated = cell(LoadStepsCount,1);
if Rotation == true
    if RotationPassed(i) == true
        for i = 2:length(Matrix)
            try
                MatrixRotated{i} = RotateAndCropMatrix(Matrix{i},RB{i},min(zShift),max(zShift),zShift(i));
            catch
                fprintf('xy rotation failed on sample %i\n', i )
                continue
            end
        end
    end
else
    for i = 2:length(Matrix)
        try
            MatrixRotated{i} = CropMatrix(Matrix{i},min(zShift),max(zShift),zShift(i));
        catch
            fprintf('Shift failed on sample %i\n', i )
            continue
        end
    end
end

MatrixRotated{1} = Matrix{1}(:,:,(1-min(zShift)):(zdim-max(zShift)));
%% Save Rotated Matrices

clear SaveFile
clear SaveFileID
SaveFile = cell(length(Matrix),1);
SaveFileID = cell(length(Matrix),1);
for i = 2:length(Matrix)
    if RotationPassed(i) == true
        SaveFile{i} = strcat(SavePath,'\Cropped+Rotated ',Files{i});
        SaveFileID{i} = fopen(SaveFile{i},'w+');
        fwrite(SaveFileID{i},reshape(MatrixRotated{i},1,[]),'uint16');
        fclose(SaveFileID{i});
    end
end

%% Calculate Residuals

if ResidualQ == true
    clear Residual
    Residual = cell(length(Matrix),1);
    if RotationPassed(i) == true
        for i = 1:length(Matrix)
            Residual{i} = ResidualCalculator(MatrixRotated{1},MatrixRotated{i});
            fprintf('xy rotation failed on sample %i\n', i )
        end
    end
end

%% Save Residuals
% Save Residual to file
if ResidualQ == true
    ResSavePath = uigetdir(Path,'Save residual .raw files');
    while strcmp(Path,ResSavePath) == true || strcmp(SavePath,ResSavePath)
        warning('You must chose a different path to not overwrite the data!')
        ResSavePath = uigetdir(Path,'Save residual .raw files');
    end

    clear ResSaveFile
    clear ResSaveFileID
    ResSaveFile = cell(length(Matrix),1);
    ResSaveFileID = cell(length(Matrix),1);
    for i = 2:length(Matrix)
        if RotationPassed(i) == true
            ResSaveFile{i} = strcat(ResSavePath,'\Residual ',Files{i});
            ResSaveFileID{i} = fopen(ResSaveFile{i},'w+');
            fwrite(ResSaveFileID{i},reshape(Residual{i},1,[]),'uint16');
            fclose(ResSaveFileID{i});
        end
    end
end

fprintf('Code finished!\nTotal run time %.2f\n', toc(tic_start)); % Display the total time taken
