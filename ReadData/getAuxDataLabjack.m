% Labjack Data is stored using the LJStreamUD application which comes with
% the installation package. It runs only on windows. This application
% stores analog data, but can also record a 16 digit digital code as an
% analog data.

% Brain Products does not have a good solution for reading Aux channels
% such as eye position, sensor readings and so on. So we use Labjack
% instead. The first analog channel is used for reading the reward signal,
% the remaining channels are for reading sensor/eye data.

% We assume that the raw data is initially stored in
% folderSourceString\data\rawData\{subjectName}{expDate}\

% We assume that the sensor data is being recorded along with the
% BrainProducts EEG data. 

function getAuxDataLabjack(subjectName,expDate,protocolName,folderSourceString,gridType,goodStimTimes,timeStartFromBaseLine,deltaT,electrodesToStore,reallignFlag)

if ~exist('folderSourceString','var');    folderSourceString ='F:';      end
if ~exist('timeStartFromBaseLine','var'); timeStartFromBaseLine= -0.55;  end
if ~exist('deltaT','var');                deltaT = 1.024;                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataLog{1,2} = subjectName;
dataLog{2,2} = gridType;
dataLog{3,2} = expDate;
dataLog{4,2} = protocolName;
dataLog{14,2} = folderSourceString;

[~,folderName]=getFolderDetails(dataLog);

makeDirectory(folderName);
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
makeDirectory(folderSegment);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read Labjack Data
[data,t] = readLabjackData(subjectName,expDate,protocolName,folderSourceString,length(electrodesToStore));

%%%%%%%%%%%%%%%%%%%%%%%%%% Compare Reward Times %%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(fullfile(folderExtract,'LL.mat'),'LL');
% rewardTimesLL = LL.endTime(LL.eotCode==0);

load(fullfile(folderExtract,'trialResults.mat'));
for i=1:length(trialEvents) %#ok<*USENS>
    if strcmp(trialEvents{i},'TE')
        pos=i;
        break;
    end
end

endTimes = trialResults(pos).times;
endVals  = trialResults(pos).value;
rewardTimesDigital = endTimes(endVals==0);

rewardChannel=1; cutOff = -0.5;
d = diff(mod(data{rewardChannel},2));
rewardTimesLabjack = t(d<cutOff);

diffTimes = rewardTimesLabjack(:) - rewardTimesDigital(:);
offset = mean(diffTimes);

disp(['Labjack data stream ahead by mean: ' num2str(offset) ', min: ' num2str(min(diffTimes)) ', max: ' num2str(max(diffTimes)) ' seconds, max d: ' num2str(1000*(max(diffTimes) - min(diffTimes))) ' ms']); 

goodStimTimesLabjack = goodStimTimes + offset;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save LFP Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveAnalogData(data,t,electrodesToStore,folderSegment,goodStimTimesLabjack,timeStartFromBaseLine,deltaT,reallignFlag);
end

% Read Data from Labjack
function [analogData,timeVals] = readLabjackData(subjectName,expDate,protocolName,folderSourceString,numChannels)

folderIn = fullfile(folderSourceString,'data','rawData',[subjectName expDate]);

timeVals = [];
for i=1:numChannels
    analogData{i} = [];
end

fileExistFlag = 1;
fileCounter = 0;

while(fileExistFlag)
    fileName = fullfile(folderIn,[subjectName expDate protocolName '_' num2str(fileCounter) '.dat']);
    
    if exist(fileName,'file')
        disp(['Reading from ' fileName]);
        [v,t]=readLabjackDataSingleFile(fileName,numChannels);
    
        timeVals = cat(1,timeVals,t);
        for j=1:numChannels
            analogData{j} = cat(1,analogData{j},v{j});
        end
        fileCounter = fileCounter+1;
    else
        fileExistFlag = 0;
    end 
end
end
function [v,t,dataDetails]=readLabjackDataSingleFile(fileName,numChannels)

fid = fopen(fileName,'r');

% Get Parameters of interest

dataDetails.dateStr = fscanf(fid,'%s',1);
dataDetails.timeStr = fscanf(fid,'%s',2);

fscanfStr = '%f';
for i=1:numChannels
    dataDetails.channelDetails{i} = fscanf(fid,'%s',5);
    fscanfStr = [fscanfStr '%f%f'];
end

dataDetails.paramsList = fscanf(fid,'%s',2*numChannels+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = fscanf(fid,fscanfStr);

t = data(1:numChannels*2+1:end);
for i=1:numChannels
    v{i} = data(i+1:numChannels*2+1:end); %#ok<*AGROW>
end

fclose(fid);
end
function saveAnalogData(data,t,analogChannelsStored,folderSegment,goodStimTimes,timeStartFromBaseLine,deltaT,reallignFlag)

Fs = round(1./(t(2)-t(1)));

% Make Diectory for storing LFP data
outputFolder = fullfile(folderSegment,'LFP');
makeDirectory(outputFolder);

analysisOnsetTimes = goodStimTimes + timeStartFromBaseLine;
numSamples = deltaT*Fs;
timeVals = timeStartFromBaseLine+ (1/Fs:1/Fs:deltaT); %#ok<NASGU>

% Now segment and store data in the outputFolder directory
totalStim = length(analysisOnsetTimes);

for i=1:length(analogChannelsStored)
    dataThisChannel = data{i};
    
    clear analogData
    analogData = zeros(totalStim,numSamples);
    
    for j=1:totalStim
        pos = find(t<=analysisOnsetTimes(j), 1, 'last' );
        analogData(j,:) = dataThisChannel(pos+1:pos+numSamples);
    end
    if reallignFlag % Added by MD 19-04-2015
        save(fullfile(outputFolder,['unallignedAinp' num2str(analogChannelsStored(i))]),'analogData');
    else
        save(fullfile(outputFolder,['ainp' num2str(analogChannelsStored(i))]),'analogData');
    end
    
end

analogInputNum = analogChannelsStored;
% Write LFP information
if exist(fullfile(outputFolder,'lfpInfo.mat'),'file')    
    clear analogChannelsStored timeVals;
    load(fullfile(outputFolder,'lfpInfo.mat'));
    analogInputNums = analogInputNum;
    save(fullfile(outputFolder,'lfpInfo.mat'),'analogChannelsStored','electrodesStored','analogInputNums','goodStimPos','timeVals');
else
    save(fullfile(outputFolder,'lfpInfo.mat'),'analogChannelsStored','timeVals');
end

end