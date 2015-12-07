% measure: 'LFP', 'Spikes', 'CSD' or 'Energy'
% timeRanges: a cell array of times - in seconds
% 'LFP' and 'Spikes' will run by default. For using the 'CSD' option, the
% CSDs must first be computed. Similarly, for using the 'Energy' option,
% meanEnergy must already be saved.

% Changes
% 1. option to save the output in a different folder specified by folderOut.
% 2. timeRanges have been replaced by timeRangeList. Each entry can have one
% or two time intervals (second is baseline). For the latter, the program
% subtracts the metric of interest (RMS,Max or Power) in the second
% interval from the first.  

% 3. More variables are now saved in rfValues .mat file such that RF
% estimation can be done without needing additional information.

% 18/6/15: Now also saving the absolute of the min value. Also addnig an
% option to put a filter

function filterStr = getValuesForRFEstimation(subjectName,expDate,protocolName,folderSourceString,gridType,measure,timeRangeList,removeAvgRef,electrodeList,folderOutString,applyFilterFlag)

if ~exist('timeRangeList','var');       timeRangeList=[];                   end
if ~exist('removeAvgRef','var');        removeAvgRef=0;                     end
if ~exist('electrodeList','var');       electrodeList=[];                   end
if ~exist('folderOutString','var');     folderOutString=folderSourceString; end
if ~exist('applyFilterFlag','var');     applyFilterFlag = 0;                end

stimPosGreaterThanOne=0; % 0 - all, 1 - all except the first stimulus in each trial, 2 - only the first stimulus of each trial

if isempty(timeRangeList)
    timeRangeList{1}={[40 100]/1000};
    timeRangeList{2}={[0 200]/1000};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% foldernames
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

folderOut2 = fullfile(folderOutString,'data',subjectName,gridType,expDate,protocolName);
folderOut1 = fullfile(folderOut2,'RFMeasures');
makeDirectory(folderOut1);
folderOut = fullfile(folderOut1,measure);
makeDirectory(folderOut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if removeAvgRef
    disp('Removing average reference');
    load(fullfile(folderSegment,'LFP','avgRef'));
    avgRef = analogData;
    fileTag = 'AvgRefRemoved';
else
    fileTag = '';
end

% Load stimulus parameters, bad trials
load(fullfile(folderExtract,'parameterCombinations.mat'));
aLength = length(aValsUnique);
eLength = length(eValsUnique);

numTimePeriods = length(timeRangeList);

load(fullfile(folderSegment,'badTrials.mat'));

if strcmp(measure,'LFP')    % Case 1 - LFP analysis
    % Get Time Ranges
    load(fullfile(folderSegment,'LFP','lfpInfo.mat'));
    Fs = round(1/(timeVals(2)-timeVals(1)));
    [timePos,timePosBL] = getTimePos(timeRangeList,timeVals);
    
    stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosGreaterThanOne);
    
    if isempty(electrodeList)
        electrodeList=analogChannelsStored;
    end
    
    for i=1:length(electrodeList)
        channelNumber = electrodeList(i);
        
        % Get LFP data
        clear signal analogData meanLFPData numStimuli
        load(fullfile(folderSegment,'LFP',['elec' num2str(channelNumber) '.mat']));

        if removeAvgRef
            analogData = analogData-avgRef;
        end
        
        if applyFilterFlag
            [analogData,filterStr] = applyFilter(analogData,Fs);
        else
            filterStr='';
        end
        
        for a=1:aLength
            for e=1:eLength
                
                clear goodPos
                goodPos = intersect(parameterCombinations{a,e,1,end,end},stimPos); %#ok<*USENS>
                goodPos = setdiff(goodPos,badTrials);
                
                if isempty(goodPos)
                    rfValsRMS(e,a,channelNumber,1:numTimePeriods)=0; %#ok<*AGROW>
                    rfValsMax(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsPower(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsMin(e,a,channelNumber,1:numTimePeriods)=0;
                    numStimuli(e,a) = 0;
                else
                    clear erp erpBL erpST
                    erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>
                    numStimuli(e,a) = length(goodPos);
                    meanLFPData(e,a,:) = erp; %#ok<*NASGU>
                    
                    for j=1:numTimePeriods
                        clear erpSegment rmsVal
                        erpSegment = erp(timePos{j});
                        rmsVal = rms(erpSegment);
                        maxVal = max(abs(erpSegment));
                        minVal = abs(min(erpSegment));
                        
                        erpSegmentBL = erp(timePosBL{j});
                        if isempty(erpSegmentBL)
                            rmsValBL=0;
                            maxValBL=0;
                            minValBL=0;
                        else
                            rmsValBL = rms(erpSegmentBL);
                            maxValBL = max(abs(erpSegmentBL));
                            minValBL = abs(min(erpSegmentBL));
                        end
                        
                        rfValsRMS(e,a,channelNumber,j) = rmsVal - rmsValBL;
                        rfValsMax(e,a,channelNumber,j) = maxVal - maxValBL;
                        rfValsPower(e,a,channelNumber,j) = rmsVal.^2 - rmsValBL.^2;
                        rfValsMin(e,a,channelNumber,j) = minVal - minValBL;
                    end
                end
            end
        end
        
        % Save Mean LFP data
        save(fullfile(folderOut,['meanLFPDataChan' num2str(channelNumber) fileTag filterStr '.mat']),'meanLFPData','timeVals','numStimuli','aValsUnique','eValsUnique');
    end
    
    % Save - numStimuli should be the same for all channels
    save(fullfile(folderOut,['rfValues' fileTag filterStr '.mat']),'rfValsRMS','rfValsMax','rfValsPower','rfValsMin','numStimuli','timeRangeList','electrodeList','aValsUnique','eValsUnique');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmp(measure,'Spikes')    % Case 2 - spike analysis
    
    filterStr = ''; % No filter
    load(fullfile(folderSegment,'LFP','lfpInfo.mat'));  % to get timeVals
    load(fullfile(folderSegment,'Spikes','spikeInfo.mat'));
    stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosGreaterThanOne);
    
    if isempty(electrodeList)
        electrodeList=neuralChannelsStored;
    end
    
    for i=1:length(electrodeList)
        channelNumber = electrodeList(i);
        SID = 0; % SourceUnitID(channelNumber);
        
        % Get Spike data
        clear spikeData neuralInfo
        load(fullfile(folderSegment,'Spikes',['elec' num2str(channelNumber) '_SID' num2str(SID) '.mat']));
        
        for a=1:aLength
            for e=1:eLength
                
                clear goodPos
                goodPos = intersect(parameterCombinations{a,e,1,end,end},stimPos);
                goodPos = setdiff(goodPos,badTrials);
                
                if isempty(goodPos)
                    rfValsRMS(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsMax(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsMean(e,a,channelNumber,1:numTimePeriods)=0;
                    
                    numStimuli(e,a) = 0;
                else
                    clear firingRate
                    [firingRate,timeValsFR] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(length(timeVals))]);
                    numStimuli(e,a) = length(goodPos);
                    meanSpikeData(e,a,:) = firingRate;
                    
                    for j=1:numTimePeriods
                        timePos = intersect(find(timeValsFR>=timeRangeList{j}{1}(1)),find(timeValsFR<timeRangeList{j}{1}(2)));
                        
                        if length(timeRangeList{j})==1
                            timePosBL = [];
                        else
                            timePosBL = intersect(find(timeValsFR>=timeRangeList{j}{2}(1)),find(timeValsFR<timeRangeList{j}{2}(2)));
                        end
                        
                        clear frSegment rmsVal maxVal meanVal frSegmentBL rmsValBL maxValBL meanValBL
                        frSegment = firingRate(timePos);
                        rmsVal = rms(frSegment);
                        maxVal = max(frSegment);
                        meanVal = mean(frSegment);
                        
                        frSegmentBL = firingRate(timePosBL);
                        if isempty(frSegmentBL)
                            rmsValBL = 0;
                            maxValBL = 0;
                            meanValBL = 0;
                        else
                            rmsValBL = rms(frSegmentBL);
                            maxValBL = max(frSegmentBL);
                            meanValBL = mean(frSegmentBL);
                        end
                        
                        rfValsRMS(e,a,channelNumber,j) = rmsVal - rmsValBL;
                        rfValsMax(e,a,channelNumber,j) = maxVal - maxValBL;
                        rfValsMean(e,a,channelNumber,j) = meanVal - meanValBL;
                    end
                end
            end
        end
        
        % Save Mean LFP data
        save(fullfile(folderOut,['meanSpikeDataChan' num2str(channelNumber) '_SID' num2str(SID) '.mat']),'meanSpikeData','timeValsFR','numStimuli','aValsUnique','eValsUnique');
    end
    
    % Save - numStimuli should be the same for all channels
    save(fullfile(folderOut,'rfValues.mat'),'rfValsRMS','rfValsMax','rfValsMean','numStimuli','timeRangeList','electrodeList','aValsUnique','eValsUnique');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmpi(measure,'CSD')    % Case 3 - CSD analysis (very similar to LFP analysis)
    
    % Get Time Ranges
    load(fullfile(folderSegment,'LFP','lfpInfo.mat'));
    [timePos,timePosBL] = getTimePos(timeRangeList,timeVals);
    
    stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosGreaterThanOne);
    
    if isempty(electrodeList)
        electrodeList=analogChannelsStored;
    end
    
    for i=1:length(electrodeList)
        channelNumber = electrodeList(i);
        
        % Get CSD data
        clear signal csdData meanCSDData numStimuli
        load(fullfile(folderSegment,'CSD',['elec' num2str(channelNumber) '.mat']));
        if removeAvgRef
            csdData = csdData-avgRef;
        end
        
        for a=1:aLength
            for e=1:eLength
                
                clear goodPos
                goodPos = intersect(parameterCombinations{a,e,1,end,end},stimPos); %#ok<*USENS>
                goodPos = setdiff(goodPos,badTrials);
                
                if isempty(goodPos)
                    rfValsRMS(e,a,channelNumber,1:numTimePeriods)=0; %#ok<*AGROW>
                    rfValsMax(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsPower(e,a,channelNumber,1:numTimePeriods)=0;
                    
                    numStimuli(e,a) = 0;
                else
                    clear erp erpBL erpST
                    erp = mean(csdData(goodPos,:),1); %#ok<*NODEF>
                    numStimuli(e,a) = length(goodPos);
                    meanCSDData(e,a,:) = erp; %#ok<*NASGU>
                    
                    for j=1:numTimePeriods
                        
                        clear erpSegment rmsVal
                        erpSegment = erp(timePos{j});
                        rmsVal = rms(erpSegment);
                        maxVal = max(abs(erpSegment));
                        
                        erpSegmentBL = erp(timePosBL{j});
                        if isempty(erpSegmentBL)
                            rmsValBL=0;
                            maxValBL=0;
                        else
                            rmsValBL = rms(erpSegmentBL);
                            maxValBL = max(abs(erpSegmentBL));
                        end
                        
                        rfValsRMS(e,a,channelNumber,j) = rmsVal - rmsValBL;
                        rfValsMax(e,a,channelNumber,j) = maxVal - maxValBL;
                        rfValsPower(e,a,channelNumber,j) = rmsVal.^2 - rmsValBL.^2;
                    end
                end
            end
        end
        
        % Save Mean LFP data
        save(fullfile(folderOut,['meanCSDDataChan' num2str(channelNumber) fileTag '.mat']),'meanCSDData','timeVals','numStimuli','aValsUnique','eValsUnique');
    end
    
    % Save - numStimuli should be the same for all channels
    save(fullfile(folderOut,['rfValues' fileTag '.mat']),'rfValsRMS','rfValsMax','rfValsPower','numStimuli','timeRangeList','electrodeList','aValsUnique','eValsUnique');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmpi(measure,'Energy')    % Case 4 - Energy analysis
    
    folderMP = fullfile(folderName,'mpAnalysis');
    
    % Get Time Ranges
    clear('energyValues','downsampledFreqVals','downsampledTimeVals','numStimuli');
    load(fullfile(folderMP,['elec' num2str(electrodeList(1))],'energyMatlab',['mEnergy_a' num2str(1) 'e' num2str(1) '.mat']));
    
    [timePos,timePosBL] = getTimePos(timeRangeList,downsampledTimeVals);
    numFreqPos = length(downsampledFreqVals);
    
    %stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosGreaterThanOne);

    for i=1:length(electrodeList)
        channelNumber = electrodeList(i);
        
        clear('rfValsRMS','rfValsMax','rfValsPower','numStimuliAll');
        
        for a=1:aLength
            for e=1:eLength
                
                % Get Energy data
                clear('energyValues','downsampledFreqVals','downsampledTimeVals','numStimuli');
                load(fullfile(folderMP,['elec' num2str(channelNumber)],'energyMatlab',['mEnergy_a' num2str(a) 'e' num2str(e) '.mat']));

                %clear goodPos
                %goodPos = intersect(parameterCombinations{a,e,1,end,end},stimPos); %#ok<*USENS>
                %goodPos = setdiff(goodPos,badTrials);
                
                if numStimuli==0
                    rfValsRMS(e,a,1:numTimePeriods,1:numFreqPos)=0;
                    rfValsMax(e,a,1:numTimePeriods,1:numFreqPos)=0;
                    rfValsPower(e,a,1:numTimePeriods,1:numFreqPos)=0;
                    
                    numStimuliAll(e,a) = 0;
                else
                    numStimuliAll(e,a) = numStimuli;
                    
                    for j=1:numTimePeriods
                        clear meanEnergy meanEnergyBL maxEnergy maxEnergyBL
                        meanEnergy = mean(energyValues(:,timePos{j}),2);
                        maxVal     = sqrt(max(energyValues(:,timePos{j}),[],2)); %sqrt of max energy

                        if isempty(timePosBL{j})
                            rfValsRMS(e,a,j,:) = sqrt(meanEnergy);
                            rfValsMax(e,a,j,:) = maxVal;
                            rfValsPower(e,a,j,:) = meanEnergy;
                        else
                            meanEnergyBL = mean(energyValues(:,timePosBL{j}),2);
                            maxValBL     = sqrt(max(energyValues(:,timePosBL{j}),[],2)); %sqrt of max energy
                        
                            rfValsRMS(e,a,j,:) = sqrt(meanEnergy) - sqrt(meanEnergyBL);
                            rfValsMax(e,a,j,:) = maxVal - maxValBL;
                            rfValsPower(e,a,j,:) = meanEnergy - meanEnergyBL;
                        end
                    end
                end
            end
        end
        
        % Save Mean LFP data
        clear numStimuli
        numStimuli=numStimuliAll;
        save(fullfile(folderOut,['rfValues' num2str(channelNumber) fileTag '.mat']),'rfValsRMS','rfValsMax','rfValsPower','numStimuli','downsampledFreqVals','timeRanges','aValsUnique','eValsUnique');
    end
    
    save(fullfile(folderOut,['rfValues' fileTag '.mat']),'numStimuli','timeRangeList','electrodeList','aValsUnique','eValsUnique');
end
end
function stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosOption)

folderExtract = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData');
load(fullfile(folderExtract,'goodStimNums.mat'));
load(fullfile(folderExtract,'stimResults.mat'));

goodStimPos = stimResults.stimPosition(goodStimNums);

if exist(fullfile(folderExtract,'validStimAfterTarget.mat'),'file')
    load(fullfile(folderExtract,'validStimAfterTarget.mat'));
    if ~isempty(validStimuliAfterTarget)
        disp(['Removing ' num2str(length(validStimuliAfterTarget)) ' stimuli after target']);
    end
    goodStimPos(validStimuliAfterTarget)=-1;  % These will be not be included in either stimPos==1 or stimPos>1
end

if stimPosOption==1
    stimPos = find(goodStimPos>1);
elseif stimPosOption==2
    stimPos = find(goodStimPos==1);
else
    stimPos = find(goodStimPos>0);
end
end
function [timePos,timePosBL] = getTimePos(timeRangeList,timeVals)

numTimePeriods = length(timeRangeList);
timePos = cell(1,numTimePeriods);
timePosBL = cell(1,numTimePeriods);

for i=1:numTimePeriods
    timeRangetmp = timeRangeList{i};
    timePos{i} = intersect(find(timeVals>=timeRangetmp{1}(1)),find(timeVals<timeRangetmp{1}(2)));
    
    if length(timeRangetmp)==1
        timePosBL{i}=[];
    else
        timePosBL{i} = intersect(find(timeVals>=timeRangetmp{2}(1)),find(timeVals<timeRangetmp{2}(2)));
    end
end
end