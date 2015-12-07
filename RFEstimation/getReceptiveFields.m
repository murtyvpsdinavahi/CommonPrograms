function getReceptiveFields(subjectName,expDate,protocolName,folderSourceString,gridType,measure,removeAvgRef,poolingOptionsList,filterStr,channelNumbers)

if ~exist('removeAvgRef','var');         removeAvgRef=0;                end
if ~exist('poolingOptionsList','var');   poolingOptionsList=1:3;        end
if ~exist('filterStr','var');            filterStr='';                  end
if ~exist('channelNumbers','var');       channelNumbers=[];             end

% foldername
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderOut = fullfile(folderName,'RFMeasures',measure);

%load(fullfile(folderExtract,'parameterCombinations.mat'));

if removeAvgRef
    fileTag = 'AvgRefRemoved';
else
    fileTag = '';
end

load(fullfile(folderOut,['rfValues' fileTag filterStr '.mat']));
if strcmpi(measure,'LFP') || strcmpi(measure,'CSD') || strcmpi(measure,'Spikes')
    numTimeRanges = size(rfValsRMS,4); %#ok<*NODEF>
elseif strncmpi(measure,'Energy',6)
    load(fullfile(folderOut,['rfValues' num2str(electrodeList(1)) fileTag '.mat']));
    numTimeRanges = size(rfValsRMS,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(channelNumbers)
    if strcmpi(measure,'LFP') || strcmpi(measure,'CSD')
        load(fullfile(folderSegment,'LFP','lfpInfo.mat'));
        channelNumbers = analogChannelsStored;
    elseif strcmpi(measure,'Spikes')
        load(fullfile(folderSegment,'Spikes','spikeInfo.mat'));
        channelNumbers = neuralChannelsStored;
    elseif strcmpi(measure,'Energy')
        channelNumbers = goodElectrodes;
    end
end
% channelNumbers=electrodeList;

if strcmpi(measure,'LFP') || strcmpi(measure,'CSD') || strcmpi(measure,'Spikes')
    for ii=1:length(poolingOptionsList)
        poolingOption = poolingOptionsList(ii);
        clear paramsRMS paramsMax paramsPower paramsMean paramsMin
        
        for i=1:length(channelNumbers)
            channelNumber = channelNumbers(i);
            
            for j=1:numTimeRanges
                clear rfValsRMStmp rfValsMaxtmp rfValsPowertmp rfValsMeantmp rfValsMintmp
                rfValsRMStmp = squeeze(rfValsRMS(:,:,channelNumber,j));
                rfValsMaxtmp = squeeze(rfValsMax(:,:,channelNumber,j));
                
                paramsRMS{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsRMStmp,poolingOption); %#ok<*AGROW,*NASGU>
                paramsMax{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsMaxtmp,poolingOption);
                
                paramsRMSScaled{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsRMStmp,poolingOption,numStimuli);
                paramsMaxScaled{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsMaxtmp,poolingOption,numStimuli);
                
                if strcmp(measure,'LFP') || strcmp(measure,'CSD')
                    rfValsPowertmp = squeeze(rfValsPower(:,:,channelNumber,j));
                    paramsPower{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsPowertmp,poolingOption);
                    paramsPowerScaled{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsPowertmp,poolingOption,numStimuli);
                    
                    rfValsMintmp = squeeze(rfValsMin(:,:,channelNumber,j));
                    paramsMin{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsMintmp,poolingOption);
                    paramsMinScaled{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsMintmp,poolingOption,numStimuli);
                    
                elseif strcmp(measure,'Spikes')
                    rfValsMeantmp = squeeze(rfValsMean(:,:,channelNumber,j));
                    paramsMean{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsMeantmp,poolingOption);
                    paramsMeanScaled{channelNumber,j} = getRFcenter(aValsUnique,eValsUnique,rfValsMeantmp,poolingOption,numStimuli);
                end
            end
        end
        
        if strcmp(measure,'LFP') || strcmp(measure,'CSD')
            save(fullfile(folderOut,['rfParams' fileTag num2str(poolingOption) filterStr '.mat']),'paramsRMS','paramsMax','paramsPower','paramsMin', ...
                'paramsRMSScaled','paramsMaxScaled','paramsPowerScaled','paramsMinScaled','aValsUnique','eValsUnique');
        elseif strcmp(measure,'Spikes')
            save(fullfile(folderOut,['rfParams' num2str(poolingOption) filterStr '.mat']),'paramsRMS','paramsMax','paramsMean', ...
                'paramsRMSScaled','paramsMaxScaled','paramsMeanScaled','aValsUnique','eValsUnique');
        end
    end
    
elseif strncmpi(measure,'Energy',6)
    
    numFreqPos = length(downsampledFreqVals);
    
    for i=1:length(channelNumbers)
        channelNumber = channelNumbers(i);
        clear rfValsRMS rfValsMax rfValsPower
        load(fullfile(folderOut,['rfValues' num2str(channelNumber) fileTag '.mat']));
 
        for ii=1:length(poolingOptionsList)
            poolingOption = poolingOptionsList(ii);
            disp([channelNumber poolingOption]);
            clear paramsRMS paramsMax paramsPower paramsRMSScaled paramsMaxScaled paramsPowerScaled
            
            for j=1:numTimeRanges
                for k=1:numFreqPos
                    clear rfValsRMStmp rfValsMaxtmp rfValsPowertmp
                    rfValsRMStmp = squeeze(rfValsRMS(:,:,j,k));
                    rfValsMaxtmp = squeeze(rfValsMax(:,:,j,k));
                    rfValsPowertmp = squeeze(rfValsPower(:,:,j,k));
                    
                    paramsRMS(j,k,:) = getRFcenter(aValsUnique,eValsUnique,rfValsRMStmp,poolingOption); %#ok<*AGROW,*NASGU>
                    paramsMax(j,k,:) = getRFcenter(aValsUnique,eValsUnique,rfValsMaxtmp,poolingOption);
                    
                    paramsRMSScaled(j,k,:) = getRFcenter(aValsUnique,eValsUnique,rfValsRMStmp,poolingOption,numStimuli);
                    paramsMaxScaled(j,k,:) = getRFcenter(aValsUnique,eValsUnique,rfValsMaxtmp,poolingOption,numStimuli);
                    
                    paramsPower(j,k,:) = getRFcenter(aValsUnique,eValsUnique,rfValsPowertmp,poolingOption);
                    paramsPowerScaled(j,k,:) = getRFcenter(aValsUnique,eValsUnique,rfValsPowertmp,poolingOption,numStimuli);
                end
            end
            
            save(fullfile(folderOut,['rfParams' num2str(channelNumber) fileTag '_' num2str(poolingOption) '.mat']),'paramsRMS','paramsMax','paramsPower', ...
                'paramsRMSScaled','paramsMaxScaled','paramsPowerScaled','aValsUnique','eValsUnique','downsampledFreqVals');
        end
    end
end