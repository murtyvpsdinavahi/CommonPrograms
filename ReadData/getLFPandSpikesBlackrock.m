% This program extracts the LFP and Spikes from the .nev file and stores
% the data in specified folders 

% Inputs
% fileName: name of the .nev file
% analogElectrodesToStore: List of analog electrodes to store.
% folderIn: folder where filename is located. Default: C:\Supratim\rawData
% folderOut: folder where the output data will be stored. Default: folderIn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program needs the following matlab files
% makeDirectory.m
% appendIfNotPresent.m
% removeIfPresent.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supratim Ray,
% September 2008

% June 2014: In previous versions we only recorded data from the electrodes
% but not input/auxilliary channels. The code is modified to allow reading
% both. "analogChannel" now is sepatared into either an "electrode" or an
% "analogInput".

% Modified by Murty V P S Dinavahi (MD) to extract Electrode positions and Sampling Frequencies for use in
% EEG
% Average referencing: added by MD 18/10/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [electrodeNums,elecSampleRate,AinpSampleRate] = getLFPandSpikesBlackrock(dataLog,folderSourceString,analogElectrodesToStore,neuralChannelsToStore,...
    timeStartFromBaseLine,deltaT,Fs,hFile,getLFP,getSpikes,LLFlag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('hFile','var');        hFile = [];                            end
if ~exist('getLFP','var');       getLFP=1;                              end
if ~exist('getSpikes','var');    getSpikes=0;                           end % Default changed to 0 by MD 30-11-2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,folderName]=getFolderDetails(dataLog);

subjectName=dataLog{1,2};
expDate=dataLog{3,2};
protocolName=dataLog{4,2};
gridType=dataLog{2,2};

fileName = [subjectName expDate protocolName '.nev'];
makeDirectory(folderName);
folderIn = fullfile(folderSourceString,'data','rawData',[subjectName expDate]);
folderExtract = fullfile(folderName,'extractedData');
makeDirectory(folderExtract);

load(fullfile(folderExtract,'goodStimNums.mat'));

% added by MD
elecSampleRate = 0;
AinpSampleRate = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Append appropriate tags
fileName = appendIfNotPresent(fileName,'.nev');

if isempty(hFile)
    % Load the appropriate DLL from pwd (MD 30-11-2015)
    dllName = fullfile(pwd,'NeuroShare','nsNEVLibrary64.dll');
%     dllName = 'C:\Users\LabComputer6\Documents\MATLAB\NeuroShare\nsNEVLibrary64.dll';
    [nsresult] = ns_SetLibrary(dllName);
    if (nsresult ~= 0)
        error('DLL was not found! Please add NeuroShare folder to the present working directory.');
    end
end

% Load data file and display some info about the file open data file

if ~isempty(hFile)
    % Get file information
    [nsresult, fileInfo] = ns_GetFileInfo(hFile);
    % Gives you entityCount, timeStampResolution and timeSpan
    if (nsresult ~= 0)
        error('Data file information did not load!');
    end
else
    disp('Getting file info...')
    [nsresult, hFile] = ns_OpenFile(fullfile(folderIn,fileName));
    if (nsresult ~= 0)
        error('Data file did not open!');
    else
        disp('Done getting file info...');
    end

    % Get file information
    [nsresult, fileInfo] = ns_GetFileInfo(hFile);
    % Gives you entityCount, timeStampResolution and timeSpan
    if (nsresult ~= 0)
        error('Data file information did not load!');
    end
end

% Build catalogue of entities
[~, entityInfo] = ns_GetEntityInfo(hFile, 1:fileInfo.EntityCount);

% List of EntityIDs needed to retrieve the information and data
%eventList = find([entityInfo.EntityType] == 1);
analogList = find([entityInfo.EntityType] == 2);
segmentList = find([entityInfo.EntityType] == 3);
neuralList = find([entityInfo.EntityType] == 4);

% How many of a particular entity do we have
% cEvent = length(eventList);
cAnalog = length(analogList);
cSegment = length(segmentList);
cNeural = length(neuralList);

% From the entityInfo.Labels, identify the channels that have been recorded
% eventLabels = strvcat(entityInfo(eventList).EntityLabel);
analogLabels = char(entityInfo(analogList).EntityLabel);
segmentLabels = char(entityInfo(segmentList).EntityLabel);
neuralLabels = char(entityInfo(neuralList).EntityLabel);

%%%%%%% Separate AnalogChannels into Electrode and AnalogInput data %%%%%%%
electrodeCount = 0;
electrodeNums = []; % added by MD 28-10-14 for no electrode condition
ainpCount = 0;

for i=1:cAnalog
    if strcmp(analogLabels(i,1:4),'elec')
        electrodeCount = electrodeCount+1;
        electrodeNums(electrodeCount) = str2num(analogLabels(i,5:end)); %#ok<*AGROW,*ST2NM>
        electrodeListIDs(electrodeCount,:) = analogList(i);
        
    elseif strcmp(analogLabels(i,1:4),'chan')
        electrodeCount = electrodeCount+1;
        electrodeNums(electrodeCount) = str2num(analogLabels(i,5:end)); %#ok<*AGROW,*ST2NM>
        electrodeListIDs(electrodeCount,:) = analogList(i);
        
    elseif strcmp(analogLabels(i,1:4),'ainp')
        ainpCount = ainpCount+1;
        analogInputNums(ainpCount) = str2num(analogLabels(i,5:end));
        analogInputListIDs(ainpCount) = analogList(i);
    end
end

segmentChannelNums = str2num(segmentLabels(:,5:end));
neuralChannelNums  = str2num(neuralLabels(:,5:end));

% Display these numbers
disp(['Total number of Analog channels recorded: ' num2str(cAnalog) ', electrodes: ' num2str(electrodeCount) ', Inp: ' num2str(ainpCount)]);
disp(['Total number of Segments recorded: ' num2str(cSegment)]);
disp(['Total number of Neurons recorded: ' num2str(cNeural)]);

% Get the desired number of samples
numberOfItems = [entityInfo.ItemCount];
analysisOnsetTimes = goodStimTimes + timeStartFromBaseLine;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% LFP Decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if getLFP && (cAnalog>0)
    
    % Set appropriate time Range
    numSamples = deltaT*Fs;
    timeVals = timeStartFromBaseLine+ (1/Fs:1/Fs:deltaT); %#ok<NASGU>

    % The optional input analogChannelToStore contains the list of channels that
    % should be recorded
    
    if electrodeCount ~= 0
        if isempty(analogElectrodesToStore)
            disp('Analog electrode list not given, taking all available electrodes.');
            electrodeListIDsStored = electrodeListIDs;
            electrodesStored = electrodeNums;
        elseif strcmp(analogElectrodesToStore,'ainp') % Added by MD 19-04-2015
            disp('Taking all available analog inputs...');
            electrodeCount = 0;
        else
            disp(['LFP from the specified ' num2str(length(analogElectrodesToStore)) ' electrodes and ' num2str(ainpCount) ' Ainp channels will be stored.']);

            count=1;
            for i=1:length(analogElectrodesToStore)
                % Check that the selected channel actually exist
                findChannelLoc = find(electrodeNums == analogElectrodesToStore(i));
                if isempty(findChannelLoc)
                    disp(['Analog Channel ' num2str(analogElectrodesToStore(i)) ' does not exist.']);
                else
                    electrodeListIDsStored(count) = electrodeListIDs(findChannelLoc);
                    electrodesStored(count) = electrodeNums(findChannelLoc);
                    count=count+1;
                end
            end
        end
        cElectrodeListIDsStored = length(electrodeListIDsStored);
        disp([num2str(cElectrodeListIDsStored) ' electrodes obtained']);
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare folders
    folderOut = fullfile(folderName,'segmentedData');
    makeDirectory(folderOut); % main directory to store both LFP and spikes

    % Make Diectory for storing LFP data
    outputFolder = fullfile(folderOut,'LFP');
    makeDirectory(outputFolder);

    % Now segment and store data in the outputFolder directory
    totalStim = length(analysisOnsetTimes);
    goodStimPos = zeros(1,totalStim);
    for i=1:totalStim
        [~,goodStimPos(i)] = ns_GetIndexByTime(hFile,analogList(1),analysisOnsetTimes(i),-1);
    end

    %%%%%%%%%%%%%%%%%%%%%%% Get data from electrodes %%%%%%%%%%%%%%%%%%%%%%
    if electrodeCount ~= 0
        if cElectrodeListIDsStored > 0
            sumData = zeros(totalStim,numSamples); % For average referencing: added by MD 18/10/2014
            hWbElec = waitbar(0,['Getting data from elec' num2str(electrodesStored(1))]); % Added by MD 18-04-2015
            for i=1:cElectrodeListIDsStored
                waitbar((i-1)/cElectrodeListIDsStored,hWbElec,['Getting data from elec' num2str(electrodesStored(i))]);
                clear analogInfo
                [~, analogInfo] = ns_GetAnalogInfo(hFile, electrodeListIDsStored(i)); %#ok<NASGU>

                clear analogData
                analogData = zeros(totalStim,numSamples);
                for j=1:totalStim
                    [~, ~, analogData(j,:)] = ns_GetAnalogData(hFile, ...
                        electrodeListIDsStored(i), goodStimPos(j)+1, numSamples);
%                     disp(['stim no. ' num2str(j)]);
                end
                sumData = sumData + analogData; % For average referencing: added by MD 18/10/2014
                save(fullfile(outputFolder,['elec' num2str(electrodesStored(i))]),'analogData','analogInfo');
            end
            averageData = sumData/cElectrodeListIDsStored; % For average referencing: added by MD 18/10/2014
            save(fullfile(outputFolder,'lfpAverage.mat'),'averageData'); % For average referencing: added by MD 18/10/2014
            elecSampleRate = analogInfo.SampleRate;
            close(hWbElec);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%% Get analog input data %%%%%%%%%%%%%%%%%%%%%%%%%
    if ainpCount>0
        hWbAinp = waitbar(0,['Getting data from ainp' num2str(analogInputNums(1))]); % Added by MD 18-04-2015
        for i=1:ainpCount
            waitbar((i-1)/ainpCount,hWbAinp,['Getting data from ainp' num2str(analogInputNums(i))])
%             disp(['ainp' num2str(analogInputNums(i))]);
            
            clear analogInfo
            [~, analogInfo] = ns_GetAnalogInfo(hFile, analogInputListIDs(i)); %#ok<NASGU>
            
            clear analogData
            analogData = zeros(totalStim,numSamples);
            for j=1:totalStim
                [~, ~, analogData(j,:)] = ns_GetAnalogData(hFile, ...
                    analogInputListIDs(i), goodStimPos(j)+1, numSamples);
%                 disp(['stim no. ' num2str(j)]);
            end
            if strcmp(analogElectrodesToStore,'ainp') % Added by MD 19-04-2015
                save(fullfile(outputFolder,['unallignedAinp' num2str(analogInputNums(i))]),'analogData','analogInfo');
            else
                save(fullfile(outputFolder,['ainp' num2str(analogInputNums(i))]),'analogData','analogInfo');
            end
        end
        AinpSampleRate = analogInfo.SampleRate;
        close(hWbAinp);
    else
        disp ('No Ainp channels found');
        analogInputNums = 0;
    end

    % Write LFP information. For backward compatibility, we also save
    % analogChannelsStored which is the list of electrode data
    if electrodeCount == 0
        electrodesStored = 0;
    end
    analogChannelsStored = electrodesStored; %#ok<NASGU>
    save(fullfile(outputFolder,'lfpInfo.mat'),'analogChannelsStored','electrodesStored','analogInputNums','goodStimPos','timeVals');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spike Decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if electrodeCount ~= 0
    if getSpikes
        % The optional input digitalChannelsToStore contains the list of channels that
        % should be recorded
        clear neuralListIDs neuralChannelsStored
        if isempty(neuralChannelsToStore)
            disp('Neural electrode list not given, taking all available electrodes.');
            neuralChannelsToStore = neuralChannelNums;
        end

        disp(['Spikes from the specified ' num2str(length(neuralChannelsToStore)) ' electrodes will be stored.']);
        count=1;
        for i=1:length(neuralChannelsToStore)
            % Check that the selected channel actually exist
            findChannelLocs = find(neuralChannelNums == neuralChannelsToStore(i));
            if isempty(findChannelLocs)
                disp(['Neural Channel ' num2str(neuralChannelsToStore(i)) ' does not exist.']);
            else
                for j=1:length(findChannelLocs)
                    neuralListIDs(count) = neuralList(findChannelLocs(j));
                    neuralChannelsStored(count) = neuralChannelNums(findChannelLocs(j));
                    count=count+1;
                end
            end
        end
        
        if count == 1 % count is 1 if no neural channel is stored in the loop above
            cNeuralListIDs = 0;
            neuralChannelsStored = [];
        else
            cNeuralListIDs = length(neuralListIDs);
        end
        disp([num2str(cNeuralListIDs) ' neurons obtained']);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Prepare folders
        folderOut = fullfile(folderName,'segmentedData');
        makeDirectory(folderOut); % main directory to store both LFP and spikes
    
        % Make Diectory for storing Spike data
        outputFolder = fullfile(folderOut,'Spikes');
        makeDirectory(outputFolder);

        % Now segment and store data in the outputFolder directory
        totalStim = length(analysisOnsetTimes);

        SourceUnitID = zeros(1,cNeuralListIDs);
        for i=1:cNeuralListIDs
            clear neuralInfo
            [~, neuralInfo] = ns_GetNeuralInfo(hFile, neuralListIDs(i));
            SourceUnitID(i) = neuralInfo.SourceUnitID;
            disp([neuralChannelsStored(i) SourceUnitID(i)]);

            clear neuralData
            [~,neuralData]  = ns_GetNeuralData(hFile, neuralListIDs(i),1,numberOfItems(neuralListIDs(i)));

            clear spikeData
            spikeData = cell(1,totalStim);
            for j=1:totalStim
                spikeData{j} = neuralData(intersect(find(neuralData>=analysisOnsetTimes(j))...
                    ,find(neuralData<analysisOnsetTimes(j)+deltaT)));
                if ~isempty(spikeData{j})
                    spikeData{j} = spikeData{j} - goodStimTimes(j);
                end
            end
            save(fullfile(outputFolder,['elec' num2str(neuralChannelsStored(i)) ... 
                '_SID' num2str(SourceUnitID(i)) '.mat']),'spikeData','neuralInfo');
        end

        % Write Spike information
        save(fullfile(outputFolder,'spikeInfo.mat'),'neuralChannelsStored','SourceUnitID');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% Segment Decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if getSpikes % Get segments if spikes are being recorded
        segmentChannelsToStore = neuralChannelsToStore;

        % The optional input segmentChannelsToStore contains the list of channels that
        % should be recorded
        clear segmentListIDs segmentChannelsStored
        if isempty(segmentChannelsToStore)
            disp('Segment list not given, taking all available electrodes.');
            segmentChannelsToStore = unique(segmentChannelNums);
        end

        disp(['Segments from the specified ' num2str(length(segmentChannelsToStore)) ' electrodes will be stored.']);
        count=1;
        for i=1:length(segmentChannelsToStore)
            % Check that the selected channel actually exist
            findChannelLocs = find(segmentChannelNums == segmentChannelsToStore(i));
            if isempty(findChannelLocs)
                disp(['Segements from Channel ' num2str(segmentChannelsToStore(i)) ' does not exist.']);
            else
                for j=1:length(findChannelLocs)
                    segmentListIDs(count) = segmentList(findChannelLocs(j));
                    segmentChannelsStored(count) = segmentChannelNums(findChannelLocs(j));
                    count=count+1;
                end
            end
        end

        cSegmentListIDs = length(segmentListIDs);
        disp([num2str(cSegmentListIDs) ' segments obtained']);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare folders
        folderOut = fullfile(folderName,'segmentedData');
        makeDirectory(folderOut); % main directory to store both LFP and spikes
        
        outputFolder = fullfile(folderOut,'Segments');
        makeDirectory(outputFolder);

        % Now segment and store data in the outputFolder directory

        for i=1:cSegmentListIDs
            clear segmentInfo
            [~, segmentInfo] = ns_GetSegmentInfo(hFile, segmentListIDs(i)); %#ok<NASGU>
            %disp([segmentChannelsStored(i)]);

            clear timeStamp segmentData sampleCount unitID
            [~,timeStamp,segmentData, sampleCount, unitID]  = ns_GetSegmentData(hFile, segmentListIDs(i),1:min(1000000,numberOfItems(segmentListIDs(i)))); %#ok<NASGU,ASGLU>

            save(fullfile(outputFolder,['elec' num2str(segmentChannelsStored(i)) '.mat']),'segmentData','segmentInfo','timeStamp','unitID','sampleCount');
        end

        % Write Spike information
        numItems = numberOfItems(segmentListIDs); %#ok<NASGU>
        save(fullfile(outputFolder,'segmentInfo.mat'),'segmentChannelsStored','numItems');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close file
if ~strcmp(analogElectrodesToStore,'ainp') % Added by MD 19-04-2015
    ns_CloseFile(hFile);
    clear mexprog;
end
end