% This program contained the following functions from the old format
% 1. saveLLDataSRC and saveLLDataGRF
% 2. getStimResultsLLSRC and getStimResultsLLGRF
% 3. saveEyeAndBehaviorDataSRC and saveEyeAndBehaviorDataGRF
% 4. getEyePositionAndBehavioralDataSRC and getEyePositionAndBehavioralDataGRF

% 5. saveEyeDataInDegSRC, saveEyeDataInDeg (GRF suffix added),
% saveEyeDataStimPosGRF and getEyeDataStimPosGRF

% Now it only contains 1, the remining ones (2-5) are in a separate program
% (saveEyePositionAndBehaviorData)

% 10 March 2015: More information is extracted from LL data file because
% the digital codes are no longer sufficient for EEG data

% 19 March 2015: Adding endTime (trialEnd).

function LLFileExistsFlag = saveLLData(subjectName,expDate,protocolName,folderSourceString,gridType)

datFileName = fullfile(folderSourceString,'data','rawData',[subjectName expDate],[subjectName expDate protocolName '.dat']);
if ~exist(datFileName,'file')
    disp('Lablib data file does not exist');
    LLFileExistsFlag = 0;
else
    disp('Working on Lablib data file ...');
    LLFileExistsFlag = 1;
%     folderName    = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
%     folderExtract = fullfile(folderName,'extractedData');
%     makeDirectory(folderExtract);
    
    if strncmpi(protocolName,'SRC',3) % SRC
        [LL,targetInfo,psyInfo,reactInfo] = getStimResultsLLSRC(subjectName,expDate,protocolName,folderSourceString); %#ok<*ASGLU,*NASGU>
        save(fullfile(folderExtract,'LL.mat'),'LL','targetInfo','psyInfo','reactInfo');
    elseif strncmpi(protocolName,'GAV',3) % GAV (Added by MD 16-03-2015)
        LL = getStimResultsLLGAV(subjectName,gridType,expDate,protocolName,folderSourceString);
%         save(fullfile(folderExtract,'LL.mat'),'LL');
    else 
        LL = getStimResultsLLGRF(subjectName,expDate,protocolName,folderSourceString);
        save(fullfile(folderExtract,'LL.mat'),'LL');
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LL,targetInfo,psyInfo,reactInfo] = getStimResultsLLSRC(subjectName,expDate,protocolName,folderSourceString)

frameISIinMS=10;

datFileName = fullfile(folderSourceString,'data','rawData',[subjectName expDate],[subjectName expDate protocolName '.dat']);

% Get Lablib data
header = readLLFile('i',datFileName);

if isfield(header,'eccentricityDeg')
    LL.azimuthDeg   = round(100*header.eccentricityDeg.data*cos(header.polarAngleDeg.data/180*pi))/100;
    LL.elevationDeg = round(100*header.eccentricityDeg.data*sin(header.polarAngleDeg.data/180*pi))/100;
else
    LL.azimuthDeg   = header.azimuthDeg.data;
    LL.elevationDeg = header.elevationDeg.data;
end

LL.sigmaDeg = header.sigmaDeg.data;
LL.spatialFreqCPD = header.spatialFreqCPD.data;
LL.orientationDeg = header.stimOrientationDeg.data;

% Stimulus properties
numTrials = header.numberOfTrials;
targetContrastIndex = [];
targetTemporalFreqIndex = [];

catchTrial = [];
instructTrial = [];
attendLoc = [];

eotCode=[];
startTime=[];

countTI=1;
countPI=1;
countRI=1;

for i=1:numTrials
    %disp(i);
    clear trials
    trials = readLLFile('t',i);
    
    thisEOT = [trials.trialEnd.data];
    eotCode = cat(2,eotCode,thisEOT);
    
    if isfield(trials,'trialStart')
        startTime = cat(2,startTime,[trials.trialStart.timeMS]);
    end
    
    %     if isfield(trials,'stimulusOn')
    %         allStimulusIndex = [allStimulusIndex [trials.stimulusOn.data]'];
    %     end
    %
    %     if isfield(trials,'stimulusOnTime')
    %         allStimulusOnTimes = [allStimulusOnTimes [trials.stimulusOnTime.timeMS]'];
    %     end
    
    if isfield(trials,'trial')
        thisCatchTrial    = [trials.trial.data.catchTrial];
        thisInstructTrial = [trials.trial.data.instructTrial];
        thisAttendLoc     = [trials.trial.data.attendLoc];
        thisTargetContrastIndex     = [trials.trial.data.targetContrastIndex];
        thisTargetTemporalFreqIndex = [trials.trial.data.targetTemporalFreqIndex];
        
        catchTrial    = cat(2,catchTrial,thisCatchTrial);
        instructTrial = cat(2,instructTrial,thisInstructTrial);
        attendLoc     = cat(2,attendLoc,thisAttendLoc);
        
        targetContrastIndex    = cat(2,targetContrastIndex,thisTargetContrastIndex);
        targetTemporalFreqIndex = cat(2,targetTemporalFreqIndex,thisTargetTemporalFreqIndex);
        
        
        % Adding the getTrialInfoLLFile (from the MAC) details here
        
        
        % Nothing to update on catch trials or instruction trials
        if (thisCatchTrial) || (thisInstructTrial) %#ok<*BDSCI,BDLGI>
            % Do nothing
        else
            
            % All the information needed to plot the targetIndex versus
            % performace plot are stored in the variable targetInfo
            
            targetInfo(countTI).trialNum          = i; %#ok<*AGROW>
            targetInfo(countTI).contrastIndex     = thisTargetContrastIndex;
            targetInfo(countTI).temporalFreqIndex = thisTargetTemporalFreqIndex;
            targetInfo(countTI).targetIndex       = trials.trial.data.targetIndex;
            targetInfo(countTI).attendLoc         = thisAttendLoc;
            targetInfo(countTI).eotCode           = thisEOT;
            countTI=countTI+1;
            
            % In addition, if the EOTCODE is correct, wrong or failed, make
            % another variable 'psyInfo' that stores information to plot the
            % psychometric functions
            
            if thisEOT <=2 % Correct, wrong or failed
                psyInfo(countPI).trialNum          = i;
                psyInfo(countPI).contrastIndex     = thisTargetContrastIndex;
                psyInfo(countPI).temporalFreqIndex = thisTargetTemporalFreqIndex;
                psyInfo(countPI).attendLoc         = thisAttendLoc;
                psyInfo(countPI).eotCode           = thisEOT;
                countPI=countPI+1;
            end
            
            % Finally, get the reaction times for correct trials and store them
            % in 'reactInfo'
            
            if thisEOT==0
                reactInfo(countRI).trialNum          = i;
                reactInfo(countRI).contrastIndex     = thisTargetContrastIndex;
                reactInfo(countRI).temporalFreqIndex = thisTargetTemporalFreqIndex;
                reactInfo(countRI).attendLoc         = thisAttendLoc;
                reactInfo(countRI).reactTime         = trials.saccade.timeMS-trials.visualStimsOn.timeMS ...
                                                       -trials.stimulusOn.data(trials.trial.data.targetIndex+1)*frameISIinMS;
                countRI=countRI+1;
            end
        end
    end
end

LL.eotCode = eotCode;
LL.instructTrial = instructTrial;
LL.catchTrial = catchTrial;
LL.attendLoc = attendLoc;
LL.targetContrastIndex = targetContrastIndex;
LL.targetTemporalFreqIndex = targetTemporalFreqIndex;

LL.startTime = startTime/1000; % in seconds
end
function LL = getStimResultsLLGRF(subjectName,expDate,protocolName,folderSourceString)

datFileName = fullfile(folderSourceString,'data','rawData',[subjectName expDate],[subjectName expDate protocolName '.dat']);

% Get Lablib data
header = readLLFile('i',datFileName);

% Stimulus properties
numTrials = header.numberOfTrials;

allStimulusIndex = [];
allStimulusOnTimes = [];
gaborIndex = [];
stimType = [];
azimuthDeg = [];
elevationDeg = [];
sigmaDeg = [];
radiusDeg = [];
spatialFreqCPD =[];
orientationDeg = [];
contrastPC = [];
temporalFreqHz = [];

eotCode=[];
myEotCode=[];
endTime=[];
startTime=[];
instructTrial=[];
catchTrial=[];
trialCertify=[];

for i=1:numTrials
    disp(['trial ' num2str(i) ' of ' num2str(numTrials)]);
    clear trials
    trials = readLLFile('t',i);
    
    if isfield(trials,'trial')
        instructTrial = [instructTrial [trials.trial.data.instructTrial]];
        catchTrial    = [catchTrial [trials.trial.data.catchTrial]];
    end
    
    if isfield(trials,'trialCertify')
        trialCertify = [trialCertify [trials.trialCertify.data]];
    end
    
    if isfield(trials,'trialEnd')
        eotCode = [eotCode [trials.trialEnd.data]];
        endTime = [endTime [trials.trialEnd.timeMS]];
    end
    
    if isfield(trials,'myTrialEnd')
        myEotCode = [myEotCode [trials.myTrialEnd.data]];
    end
    
    if isfield(trials,'trialStart')
        startTime = [startTime [trials.trialStart.timeMS]];
    end
    
    if isfield(trials,'stimulusOn')
        allStimulusIndex = [allStimulusIndex [trials.stimulusOn.data]'];
    end
    
    if isfield(trials,'stimulusOnTime')
        allStimulusOnTimes = [allStimulusOnTimes [trials.stimulusOnTime.timeMS]'];
    end
    
    if isfield(trials,'stimDesc')
        gaborIndex = [gaborIndex [trials.stimDesc.data.gaborIndex]];
        stimType = [stimType [trials.stimDesc.data.stimType]];
        azimuthDeg = [azimuthDeg [trials.stimDesc.data.azimuthDeg]];
        elevationDeg = [elevationDeg [trials.stimDesc.data.elevationDeg]];
        sigmaDeg = [sigmaDeg [trials.stimDesc.data.sigmaDeg]];
        if isfield(trials.stimDesc.data,'radiusDeg')
            radiusExists=1;
            radiusDeg = [radiusDeg [trials.stimDesc.data.radiusDeg]];
        else
            radiusExists=0;
        end
        spatialFreqCPD = [spatialFreqCPD [trials.stimDesc.data.spatialFreqCPD]];
        orientationDeg = [orientationDeg [trials.stimDesc.data.directionDeg]];
        contrastPC = [contrastPC [trials.stimDesc.data.contrastPC]];
        temporalFreqHz = [temporalFreqHz [trials.stimDesc.data.temporalFreqHz]];
    end
end

% Sort stim properties by stimType
numGabors = length(unique(gaborIndex));
for i=1:numGabors
    gaborIndexFromStimulusOn{i} = find(allStimulusIndex==i-1);
    gaborIndexFromStimDesc{i} = find(gaborIndex==i-1);
end

if isequal(gaborIndexFromStimDesc,gaborIndexFromStimulusOn)
    for i=1:numGabors
        LL.(['time' num2str(i-1)]) = allStimulusOnTimes(gaborIndexFromStimulusOn{i});
        LL.(['stimType' num2str(i-1)]) = stimType(gaborIndexFromStimulusOn{i});
        LL.(['azimuthDeg' num2str(i-1)]) = azimuthDeg(gaborIndexFromStimulusOn{i});
        LL.(['elevationDeg' num2str(i-1)]) = elevationDeg(gaborIndexFromStimulusOn{i});
        LL.(['sigmaDeg' num2str(i-1)]) = sigmaDeg(gaborIndexFromStimulusOn{i});
        if radiusExists
            LL.(['radiusDeg' num2str(i-1)]) = radiusDeg(gaborIndexFromStimulusOn{i});
        end
        LL.(['spatialFreqCPD' num2str(i-1)]) = spatialFreqCPD(gaborIndexFromStimulusOn{i});
        LL.(['orientationDeg' num2str(i-1)]) = orientationDeg(gaborIndexFromStimulusOn{i});
        LL.(['contrastPC' num2str(i-1)]) = contrastPC(gaborIndexFromStimulusOn{i});
        LL.(['temporalFreqHz' num2str(i-1)]) = temporalFreqHz(gaborIndexFromStimulusOn{i});
    end
else
    error('Gabor indices from stimuluOn and stimDesc do not match!!');
end

LL.eotCode = eotCode;
LL.myEotCode = myEotCode;
LL.startTime = startTime/1000; % in seconds
LL.endTime = endTime/1000; % in seconds
LL.instructTrial = instructTrial;
LL.catchTrial = catchTrial;
LL.trialCertify = trialCertify;

end