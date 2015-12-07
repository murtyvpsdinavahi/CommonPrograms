% We check, for each trial, whether there is a value greater than 
% 1. threshold times the std dev. 
% 2. maxLimit
% If yes, that trial is marked as a bad trial

% 15/06/15: Adding putFilterFlag. Data sometimes shows a low frequency
% drift. We put a high-pass filter to remove this drift before applying the
% thresholds

function [allBadTrials,badTrials] = findBadTrialsWithLFP(monkeyName,expDate,protocolName,folderSourceString,gridType,checkTheseElectrodes,threshold,maxLimit,showElectrodes,minLimit,saveDataFlag,applyFilterFlag)

if ~exist('checkTheseElectrodes','var');     checkTheseElectrodes = [33 12 80 63 44];   end
if ~exist('folderSourceString','var');       folderSourceString = 'G:';                 end
if ~exist('threshold','var');                threshold = 6;                             end
if ~exist('minLimit','var');                 minLimit = -2000;                          end
if ~exist('saveDataFlag','var');             saveDataFlag = 1;                          end
if ~exist('applyFilterFlag','var');          applyFilterFlag = 0;                       end

folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);
folderSegment = fullfile(folderName,'segmentedData');

numElectrodes = length(checkTheseElectrodes);

allBadTrials = cell(1,numElectrodes);

for i=1:numElectrodes
    electrodeNum=checkTheseElectrodes(i);
    load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeNum) '.mat']));
    
    if applyFilterFlag
        analogData = applyFilter(analogData);
    end
    
    numTrials = size(analogData,1); %#ok<*NODEF>
    meanData = mean(analogData,2)';
    stdData  = std(analogData,[],2)';
    maxData  = max(analogData,[],2)';
    minData  = min(analogData,[],2)';
    
    clear tmpBadTrials tmpBadTrials2
    tmpBadTrials = unique([find(maxData > meanData + threshold * stdData) find(minData < meanData - threshold * stdData)]);
    tmpBadTrials2 = unique(find(maxData > maxLimit));
    tmpBadTrials3 = unique(find(minData < minLimit));
    allBadTrials{i} = unique([tmpBadTrials tmpBadTrials2 tmpBadTrials3]);
end

badTrials=allBadTrials{1};
for i=1:numElectrodes
    badTrials=intersect(badTrials,allBadTrials{i}); % in the previous case we took the union
end

disp(['total Trials: ' num2str(numTrials) ', bad trials: ' num2str(badTrials)]);

for i=1:numElectrodes
    if length(allBadTrials{i}) ~= length(badTrials)
        disp(['Bad trials for electrode ' num2str(checkTheseElectrodes(i)) ': ' num2str(length(allBadTrials{i}))]);
    else
        disp(['Bad trials for electrode ' num2str(checkTheseElectrodes(i)) ': all (' num2str(length(badTrials)) ')']);
    end
end

if saveDataFlag
    disp(['Saving ' num2str(length(badTrials)) ' bad trials']);
    save(fullfile(folderSegment,'badTrials.mat'),'badTrials','checkTheseElectrodes','threshold','maxLimit');
else
    disp('Bad trials will not be saved..');
end

load(fullfile(folderSegment,'LFP','lfpInfo.mat'));

lengthShowElectrodes = length(showElectrodes);
if ~isempty(showElectrodes)
    for i=1:lengthShowElectrodes
        
        if lengthShowElectrodes>1
            subplot(lengthShowElectrodes,1,i);
        else
            subplot(2,1,1);
        end
        channelNum = showElectrodes(i);

        clear signal analogData
        load(fullfile(folderSegment,'LFP',['elec' num2str(channelNum) '.mat']));
        
        if applyFilterFlag
            analogData = applyFilter(analogData);
        end
    
        if numTrials<4000
            plot(timeVals,analogData(setdiff(1:numTrials,badTrials),:),'color','k');
            hold on;
        else
            disp('More than 4000 trials...');
        end
        if ~isempty(badTrials)
            plot(timeVals,analogData(badTrials,:),'color','g');
        end
        title(['electrode ' num2str(channelNum)]);
        axis tight;
        
        if lengthShowElectrodes==1
            subplot(2,1,2);
            plot(timeVals,analogData(setdiff(1:numTrials,badTrials),:),'color','k');
            axis tight;
        end
    end
end
end