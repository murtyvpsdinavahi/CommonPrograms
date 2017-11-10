% This program extracts EEG data from the raw (.eeg, .vhdr, .vmrk) files
% and stores the data in specified folders. The format is the same as for Blackrock data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program needs the following matlab files
% makeDirectory.m
% appendIfNotPresent.m
% removeIfPresent.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray,
% March 2015

function electrodesStored = getEEGDataBrainProducts(subjectName,expDate,protocolName,dataFolder,gridType,...
    goodStimTimes,timeStartFromBaseLine,deltaT,notchLineNoise,reRefFlag,refElec,decimateFlag,decimationFactor)

if ~exist('notchLineNoise','var'); notchLineNoise = 0; end;
if ~exist('reRefFlag','var'); reRefFlag = 0; end;
if ~exist('decimateFlag','var'); decimateFlag = 0; end;
if ~exist('decimationFactor','var'); decimationFactor = 1; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = [subjectName expDate protocolName '.vhdr'];

dataLog{1,2} = subjectName;
dataLog{2,2} = gridType;
dataLog{3,2} = expDate;
dataLog{4,2} = protocolName;
dataLog{14,2} = dataFolder;

[~,folderName]=getFolderDetails(dataLog);

makeDirectory(folderName);
folderIn = fullfile(dataFolder,[subjectName expDate]);
folderExtract = fullfile(folderName,'extractedData');
makeDirectory(folderExtract);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use EEGLAB plugin "bva-io" to read the file
eegInfo = pop_loadbv(folderIn,fileName,[],[]);

cAnalog = eegInfo.nbchan;
Fs = eegInfo.srate/decimationFactor; % added by MD; default is 1
analogInputNums = 1:cAnalog;
disp(['Total number of Analog channels recorded: ' num2str(cAnalog)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% EEG Decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysisOnsetTimes = goodStimTimes + timeStartFromBaseLine;
% times = (eegInfo.times/1000);


if (cAnalog>0)
    
    % Set appropriate time Range
    numSamples = deltaT*Fs;
    timeVals = timeStartFromBaseLine+ (1/Fs:1/Fs:deltaT);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare folders
    folderOut = fullfile(folderName,'segmentedData');
    makeDirectory(folderOut); % main directory to store EEG Data
    
    % Make Diectory for storing LFP data
    outputFolder = fullfile(folderOut,'LFP'); % Still kept in the folder LFP to be compatible with Blackrock data
    makeDirectory(outputFolder);
    
    % Now segment and store data in the outputFolder directory
    totalStim = length(analysisOnsetTimes);
    goodStimPos = zeros(1,totalStim);
    
    % Decimation [MD]
    if decimateFlag
        disp(['Decimating data from ',num2str(Fs*decimationFactor),' Hz to ',num2str(Fs),' Hz...']);
        for iD = 1:cAnalog
            eegInfo.data(iD,:) = decimate(eegInfo.data(iD,:),decimationFactor);
        end
    end
    
    times = 0:1/Fs:(size(eegInfo.data,2)-1)/Fs;
    
    for i=1:totalStim
        goodStimPos(i) = find(times>analysisOnsetTimes(i),1);
    end
    
    hW = waitbar(0,'Analysing electrodes...');
    
    
    
    % Rereferencing [MD]
    if reRefFlag
        disp(['Re-referencing data to elec',num2str(refElec),'...']);
        refData = eegInfo.data(refElec,:);
        for iR = 1:cAnalog
            eegInfo.data(iR,:) = eegInfo.data(iR,:) - refData;
        end
        eegInfo.data(refElec,:) = (-1).*refData;
    end
    
    % Notch filters [MD]; taken from Vinay's method in blackrock extraction
    if notchLineNoise
        for iN = 1:cAnalog
            N = length(eegInfo.data(iN,:)); % length of the full length analog data

            F = Fs/N*(0:N-1); % Span of frequencies to be considered

            clear notchRange
            notchRange{1} = find(F>48 & F<52); % 50 Hz line noise            
            notchRange{2} = find(F>98 & F<102); % 1st harmonic
            notchRange{3} = find(F>148 & F<152); % 2nd harmonic

            notchRange{4} = find(F>(Fs-52) & F<(Fs-48));
            notchRange{5} = find(F>(Fs-102) & F<(Fs-98));
            notchRange{6} = find(F>(Fs-152) & F<(Fs-148));

            % Take FFT and notch it!
            analogDataFullfft = fft(eegInfo.data(iN,:));
            analogDataFullfft(notchRange{1}) = 0;
            analogDataFullfft(notchRange{2}) = 0;
            analogDataFullfft(notchRange{3}) = 0;
            analogDataFullfft(notchRange{4}) = 0;
            analogDataFullfft(notchRange{5}) = 0;
            analogDataFullfft(notchRange{6}) = 0;

            % Take inverse fft to get the notched analog data
            analogDataFullNotched(iN,:) = ifft(analogDataFullfft);
        end
        eegInfo.data = analogDataFullNotched;
    end

    for i=1:cAnalog
%         disp(['elec' num2str(analogInputNums(i))]);
        waitbar(i/cAnalog,hW,['Analysing electrode ' num2str(i) ' of ' num2str(cAnalog) ' electrodes']);
        
        clear analogData
        analogData = zeros(totalStim,numSamples);
        for j=1:totalStim
            analogData(j,:) = eegInfo.data(analogInputNums(i),goodStimPos(j)+1:goodStimPos(j)+numSamples);
        end
        analogInfo = eegInfo.chanlocs(analogInputNums(i)); %#ok<*NASGU>
        save(fullfile(outputFolder,['elec' num2str(analogInputNums(i)) '.mat']),'analogData','analogInfo');
    end
    close(hW);
    
    % Write LFP information. For backward compatibility, we also save
    % analogChannelsStored which is the list of electrode data
    electrodesStored = analogInputNums;
    analogChannelsStored = electrodesStored;
    save(fullfile(outputFolder,'lfpInfo.mat'),'analogChannelsStored','electrodesStored','analogInputNums','goodStimPos','timeVals');
end