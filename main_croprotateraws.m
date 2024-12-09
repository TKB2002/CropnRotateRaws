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

% RJS: Done to also add in files in the current directory
addpath(genpath(cd));

%% Checking Computer
% Checks whether the computer and MATLAB version are the same as that this
% code was developed for. Will only print out warnings. Added by RJS.
checkVersion;



%%  USER INPUTS

% Step 2
% Import the dimentions of your .raw file.
xdim =1000; % Units in Pixels
ydim =1000; % Units in Pixels
zdim =1000; % Units in Pixels

i8 = false;
i16 = true;

% Import the zslice where your sample begins, everything else outside of
% this will be cropped.
Volume1StartOfSample = 0; %Z slice where the sample starts

% Step 2.1
% Does your .raw include filler around the material(i.e air or sample holder).
Cropping = false; % (false = No Cropping/ no filler, true= cropping/ has filler)

% Step 2.1.1
% If your sample needs cropping, specify the pixels to be included
Croppedxdim = 800:1600; % Units in Pixels
Croppedydim = 800:1600; % Units in Pixels
Croppedzdim = 800:1600; % Units in Pixels

% Step 2.2
% Does your sample include a tomograph boarder
% Uses mode to remove boarders, increasing the sensitivity of z-Shift
% calculations
Boarder = false; % (true = has a boarder, false = no boarder);

% Step 2.3
% Is their zShift in your sample?
Shift = true; % (true = their is a shift, false = no shift)

Normilzer = true;
zShiftBinning = 1;

% Step 2.4
% Does your sample require xy rotation.
Rotation = true;

% Step 2.4.1
% Specify the rotation meathod
Method = false; % (true = affine warping, false = rigid warping)

% Step 3
% Calculate and save Residual?
ResidualQ = true;

% Step 4
% Clear Inital Volume from RAM after rotation
Clear = true;

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
    warning('You must chose a different path to not overwrite the data!') % Displays warning if path chosen = Source path
    SavePath = uigetdir(Path,'Save registered .raw files'); % Makes user change path until, path is different from source
end

if ResidualQ == true % Of user selected option for Residual Calulcation is on
    ResSavePath = uigetdir(Path,'Save residual .raw files');
    while strcmp(Path,ResSavePath) == true || strcmp(SavePath,ResSavePath)
        warning('You must chose a different path to not overwrite the data!')
        ResSavePath = uigetdir(Path,'Save residual .raw files');
    end
end
%% Loading in the data and formatting it

% Clearing necessary variable names
clear OpenedFiles Volume Matrix

% Initalising
LoadStepsCount =length(Files); % Initalise simplified variable representing the amount of Load Steps

OpenedFiles = cell(LoadStepsCount,1); % Initalise a cell containing file infomation
Volume = cell(LoadStepsCount,1); % Initalise a cell containing the .raw data in vector form
Matrix = cell(LoadStepsCount,1); % Initalise a cell containing the .raw data in matrix form

fprintf('Importing data');
tic_start = tic; % Start Timer

for j = 1:LoadStepsCount
    OpenedFiles{j} = fopen(strcat(Path,Files{j})); % Get the file infomation
    if i16 == true
        Volume{j} = fread(OpenedFiles{j},xdim*ydim*zdim,'uint16'); % Read the .raw data into RAM
        Matrix{j} = uint16(reshape(Volume{j},xdim,ydim,[])); % Format the .raw data into Matrix form
        Volume{j} = [];
    elseif i8 == true
        Volume{j} = fread(OpenedFiles{j},xdim*ydim*zdim,'uint8'); % Read the .raw data into RAM
        Matrix{j} = uint8(reshape(Volume{j},xdim,ydim,[])); % Format the .raw data into Matrix form
        Volume{j} = [];
    end

    %Cropping Sample
    if Cropping == true %If the user specified to crop the sample in step 2.1
        Matrix{j} = Matrix{j}(Croppedxdim,Croppedydim,Croppedzdim); % Replace the Matrix data with data of a specific part.
    end
end


clear Volume % Clear the redundant .raw data in vector form.

% This step is currently using double the ram it should be because it loads
% every loadstep twice when it should only load 1 loadstep at a time twice.
% Work is going into optimisng this

ImportTime = toc(tic_start); % Save the time required for the code to import and format data
fprintf('Data Imported %.2f\n', ImportTime); % Display time taken for data being imported

%% Removing tomograph boarder
if Boarder == true
    for i = 1:LoadStepsCount
        RemovedBoarderGreyLevelValue = mode(Matrix{i}); % Find most common grey level value
        Matrix{i}(Matrix{i}==RemovedBoarderGreyLevelValue) = NaN; % Set most common grey level value to NaN
    end
end
%% Calculating the zShift of each sample
% 2 meathods WIP

% Clearing necessary variable names
clear zShift
clear Offset

ProcessTimes = nan(LoadStepsCount,1);

if Shift == true
    % Initalising
    zShift = nan(LoadStepsCount,1); % Initalise a 1D array of zShifts relative to the first sample
    zShiftError = nan(LoadStepsCount,1);

    %Setting the first value to the reference
    zShift(1) = 0;
    zShiftError(1) = 0;

    for i = 2:LoadStepsCount
        tic
        fprintf('z Shift of Sample %d \n',i)
        [zShift(i),zShiftError(i)] = zshift(Matrix{1},Matrix{i},zShiftBinning,Normilzer); % Calculate zShift with novel function, standard deviation energy reduction meathod
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

    for i = 2:LoadStepsCount
        try
            % Comparing a z-shifted slice of sample i to its equivilent in
            % sample 1.
            % Exporting the rotation properties (RB) and quality of fit
            % (2D-Residual and SSIM)
            [RB{i},MeanResI(i),MeanResF(i),SSIMi(i),SSIMf(i)] = TwoDRigidBodyRemoverAffine(Matrix{1},Matrix{i},zShift(i)); 
            
            % Breaking down the Rotation Matrix
            Affine11(i) = RB{i}.T(1,1);
            Affine12(i) = RB{i}.T(1,2);
            Affine21(i) = RB{i}.T(2,1);
            Affine22(i) = RB{i}.T(2,2);
            Affine31(i) = RB{i}.T(3,1);
            Affine32(i) = RB{i}.T(3,2);
            RotationPassed(i) = true;
        catch % If rotation function spits out an error
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
    
    RotaationAngle(1) = 0;
    xShift(1) = 0;
    yShift(1) = 0;

    for i = 2:length(Files)
        try

            % Comparing a z-shifted slice of sample i to its equivilent in
            % sample 1.
            % Exporting the rotation properties (RB) and quality of fit
            % (2D-Residual and SSIM)
            [RB{i},MeanResI(i),MeanResF(i),SSIMi(i),SSIMf(i)] = TwoDRigidBodyRemoverRigid(Matrix{1},Matrix{i},zShift(i));
            
            % Breaking down the Rotation Matrix
            xShift(i) = RB{i}.T(3,1);
            yShift(i) = RB{i}.T(3,2);
            RotationAngle(i) = acos(RB{i}.T(1,1))*57.2958; % Converting Radians to degrees
            RotationPassed(i) = true;
        catch % If rotation function spits out an error
            FailedXY(i) = true;
            fprintf('xy rotation detection failed on sample %i\n', i )
            RotationPassed(i) = false;
            continue
        end
    end
end
%% Plotting calculated shifts
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
if Clear == false
    % Clearing necessary variable names
    clear MatrixRotated
    
    % Initalising
    MatrixRotated = cell(LoadStepsCount,1);
    
    
    if Rotation == true % Where user has selected their is rotation
        if RotationPassed(i) == true % Where rotation has worked
            for i = 2:length(Matrix)
                try
                    MatrixRotated{i} = RotateAndCropMatrix(Matrix{i},RB{i},min(zShift),max(zShift),zShift(i),16); % Using rotation properties to rotate Matrix
                    PlotProgress(i-1,LoadStepsCount-1,'Rotation Progress',50)
                catch
                    fprintf('xy rotation failed on sample %i\n', i )
                    continue
                end
            end
        end
    else
        for i = 2:length(Matrix)
            try
                MatrixRotated{i} = CropMatrix(Matrix{i,1},min(zShift),max(zShift),zShift(i)); % z-cropping sample so the start of the sample matches up
            catch
                fprintf('Shift failed on sample %i\n', i )
                continue
            end
        end
    end
    
    MatrixRotated{1} = Matrix{1}(:,:,(1-min(zShift)):(zdim-max(zShift))); % Cropping loadstep 1
end
%% Save Rotated Matrices
if Clear == false
    % Clearing necessary variable names
    clear SaveFile
    clear SaveFileID
    
    % Initalising
    SaveFile = cell(length(Matrix),1);
    SaveFileID = cell(length(Matrix),1);
    
    for i = 2:length(Matrix)
        if RotationPassed(i) == true % Where rotation has worked
            SaveFile{i} = strcat(SavePath,'\Cropped+Rotated ',Files{i}); % Saving files at the User Selected File Directory with characters 'Cropped+Rotated' added
            SaveFileID{i} = fopen(SaveFile{i},'w+');
            if i8 == true
                fwrite(SaveFileID{i},reshape(uint8(MatrixRotated{i}),1,[]),'uint8');
            end
            if i16 == true
                fwrite(SaveFileID{i},reshape(uint16(MatrixRotated{i}),1,[]),'uint16');
            end
            fclose(SaveFileID{i});
            PlotProgress(i-1,LoadStepsCount-1,'Saving Progress',50)
        end
    end
end
%% Calculate Residuals
if Clear == false
    if ResidualQ == true % Of user selected option for Residual Calulcation is on
    
        % Clearing necessary variable names
        clear Residual
    
        % Initalising
        Residual = cell(length(Matrix),1);
    
        for i = 1:length(Matrix)
            if RotationPassed(i) == true     
                Residual{i} = ResidualCalculator(MatrixRotated{1},MatrixRotated{i});
                PlotProgress(i-1,LoadStepsCount-1,'Resdiual Calculation Progress',50)
            end
        end
    end
end
%% Save Residuals
% Save Residual to file
if Clear == false
    if ResidualQ == true
    
        % Clearing necessary variable names    
        clear ResSaveFile
        clear ResSaveFileID
    
        % Initalising
        ResSaveFile = cell(length(Matrix),1);
        ResSaveFileID = cell(length(Matrix),1);
    
        for i = 2:length(Matrix)
            if RotationPassed(i) == true
    
                ResSaveFile{i} = strcat(ResSavePath,'\Residual ',Files{i});  % Saving files at the User Selected File Directory with characters 'Residual' added
                ResSaveFileID{i} = fopen(ResSaveFile{i},'w+');
                if i8 == true
                    fwrite(ResSaveFileID{i},reshape(Residual{i},1,[]),'uint8');
                end
                if i16 == true
                    fwrite(ResSaveFileID{i},reshape(Residual{i},1,[]),'uint16');
                end
                fclose(ResSaveFileID{i});
                PlotProgress(i-1,LoadStepsCount-1,'Residual Saving Progress',50)
            end
        end
    end
    
    fprintf('Code finished!\nTotal run time %.2f\n', toc(tic_start)); % Display the total time taken
end

%% Post processing for Memory Efficency
if Clear == true
    % Clearing necessary variable names
    clear MatrixRotated
    clear SaveFile
    clear SaveFileID
    clear Residual
    clear ResSaveFile
    clear ResSaveFileID    
    % Initalising
    MatrixRotated = cell(LoadStepsCount,1);
    SaveFile = cell(length(Matrix),1);
    SaveFileID = cell(length(Matrix),1); 
    Residual = cell(length(Matrix),1);
    ResSaveFile = cell(length(Matrix),1);
    ResSaveFileID = cell(length(Matrix),1);

    % Initalising
   
    MatrixRotated{1} = Matrix{1}(:,:,(1-min(zShift)):(zdim-max(zShift))); % z Cropping loadstep 1
    SaveFile{1} = strcat(SavePath,'\Cropped+Rotated ',Files{1}); % Saving files at the User Selected File Directory with characters 'Cropped+Rotated' added
    SaveFileID{1} = fopen(SaveFile{1},'w+');
    if i8 == true
        fwrite(SaveFileID{1},reshape(uint8(MatrixRotated{1}),1,[]),'uint8');
    end
    if i16 == true
        fwrite(SaveFileID{1},reshape(uint16(MatrixRotated{1}),1,[]),'uint16');
    end
    fclose(SaveFileID{1});
    

           

    if Rotation == true % Where user has selected their is rotation
        if RotationPassed(i) == true % Where rotation has worked
            for i = 2:length(Matrix)
                % Rotate Matrix
                try
                    MatrixRotated{i} = RotateAndCropMatrix(Matrix{i},RB{i},min(zShift),max(zShift),zShift(i),16); % Using rotation properties to rotate Matrix
                    
                catch
                    fprintf('xy rotation failed on sample %i\n', i )
                    continue
                end
                % Save Rotated Files
                SaveFile{i} = strcat(SavePath,'\Cropped+Rotated ',Files{i}); % Saving files at the User Selected File Directory with characters 'Cropped+Rotated' added
                SaveFileID{i} = fopen(SaveFile{i},'w+');
                if i8 == true
                    fwrite(SaveFileID{i},reshape(uint8(MatrixRotated{i}),1,[]),'uint8');
                end
                if i16 == true
                    fwrite(SaveFileID{i},reshape(uint16(MatrixRotated{i}),1,[]),'uint16');
                end
                fclose(SaveFileID{i});

                %Calculate Residuals    
                Residual{i} = ResidualCalculator(MatrixRotated{1},MatrixRotated{i});

                %Save Residual Files
                ResSaveFile{i} = strcat(ResSavePath,'\Residual ',Files{i});  % Saving files at the User Selected File Directory with characters 'Residual' added
                ResSaveFileID{i} = fopen(ResSaveFile{i},'w+');
                if i8 == true
                    fwrite(ResSaveFileID{i},reshape(Residual{i},1,[]),'uint8');
                end
                if i16 == true
                    fwrite(ResSaveFileID{i},reshape(Residual{i},1,[]),'uint16');
                end
                fclose(ResSaveFileID{i});
                PlotProgress(i-1,LoadStepsCount-1,'Rotation Progress',50)
            end      
        end      
    else
        for i = 2:length(Matrix)
            try
                MatrixRotated{i} = CropMatrix(Matrix{i,1},min(zShift),max(zShift),zShift(i)); % z-cropping sample so the start of the sample matches up
            catch
                fprintf('Shift failed on sample %i\n', i )
                continue
            end
            % Save Rotated Files
            SaveFile{i} = strcat(SavePath,'\Cropped+Rotated ',Files{i}); % Saving files at the User Selected File Directory with characters 'Cropped+Rotated' added
            SaveFileID{i} = fopen(SaveFile{i},'w+');
            if i8 == true
                fwrite(SaveFileID{i},reshape(uint8(MatrixRotated{i}),1,[]),'uint8');
            end
            if i16 == true
                fwrite(SaveFileID{i},reshape(uint16(MatrixRotated{i}),1,[]),'uint16');
            end
            fclose(SaveFileID{i});

            %Calculate Residuals    
            Residual{i} = ResidualCalculator(MatrixRotated{1},MatrixRotated{i});

            %Save Residual Files
            ResSaveFile{i} = strcat(ResSavePath,'\Residual ',Files{i});  % Saving files at the User Selected File Directory with characters 'Residual' added
            ResSaveFileID{i} = fopen(ResSaveFile{i},'w+');
            if i8 == true
                fwrite(ResSaveFileID{i},reshape(Residual{i},1,[]),'uint8');
            end
            if i16 == true
                fwrite(ResSaveFileID{i},reshape(Residual{i},1,[]),'uint16');
            end
            fclose(ResSaveFileID{i});
            PlotProgress(i-1,LoadStepsCount-1,'Rotation Progress',50)
        end
     end
 end

    
    