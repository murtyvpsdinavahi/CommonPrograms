function displayLLInforAllDaysGRF(subjectName,expDates,protocolNames,folderSourceString,gridType,removeBreaks,removeIgnores)

if ~exist('removeBreaks','var');            removeBreaks=0;             end
if ~exist('removeIgnores','var');           removeIgnores=0;            end

if ~iscell(expDates)
    [expDates,protocolNames] = convertToCell(expDates,protocolNames);
end

fileStr = getFileStr(subjectName,expDates,protocolNames,folderSourceString,gridType);

% Choose day
dayString = getDayString(subjectName,expDates,protocolNames);
hChooseDay = uicontrol('Unit','Normalized','Position',[0.85 0.9 0.1 0.1],...
    'Style','popup','String',dayString,'Callback',{@plotData_callBack});

% choose whether to use or ignore breaks and ignores
hBreaks = uicontrol('Unit','Normalized','Position',[0.75 0.95 0.1 0.05],...
    'Style','checkbox','String','remove Breaks','Callback',{@plotDataLLBar_callBack});
hIgnores = uicontrol('Unit','Normalized','Position',[0.65 0.95 0.1 0.05],...
    'Style','checkbox','String','remove Ignores','Callback',{@plotDataLLBar_callBack});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot LLBars of all days
LLBarGrid = [0.05 0.65 0.9 0.225]; 
numDays = length(expDates);
dX = min(0.035,LLBarGrid(3)/(numDays+2));

% MakeTable

tableLabelGrid = [0 0.45 LLBarGrid(1) 0.2];
tableLabels{1,1} = 'Date';
tableLabels{2,1} = 'Azimuth';
tableLabels{3,1} = 'Elevation';
tableLabels{4,1} = 'Orientation';
tableLabels{5,1} = 'Spatial F';
tableLabels{6,1} = 'Sigma';
tableLabels{7,1} = 'Day Num';
makeTable(tableLabelGrid,tableLabels);

tableGrid = [LLBarGrid(1) 0.45 numDays*dX 0.2];

tableEntries = cell(7,numDays);
for i=1:numDays
    tableEntries{1,i} = [expDates{i}(1:2) '/' expDates{i}(3:4)];
    tableEntries{2,i} = num2str(round(100*str2double(fileStr(i).azimuth))/100);
    tableEntries{3,i} = num2str(round(100*str2double(fileStr(i).elevation))/100);
    tableEntries{4,i} = num2str(fileStr(i).orientationDeg0);
    tableEntries{5,i} = num2str(fileStr(i).spatialFreqCPD0);
    tableEntries{6,i} = num2str(fileStr(i).sigmaDeg0);
    tableEntries{7,i} = num2str(i);
end
makeTable(tableGrid,tableEntries);

%%% Plotting with default values of removeBreaks and removeIgnores
for i=1:numDays
    gridPositionLLBar = [LLBarGrid(1)+(i-1)*dX LLBarGrid(2) dX LLBarGrid(4)];
    displayLLbarGRF(subjectName,expDates{i},protocolNames{i},folderSourceString,gridType,gridPositionLLBar,removeBreaks,removeIgnores);
end

gridPositionLLBar = [LLBarGrid(1)+(numDays+1)*dX LLBarGrid(2) dX LLBarGrid(4)];
displayLLbarGRF(subjectName,expDates,protocolNames,folderSourceString,gridType,gridPositionLLBar,removeBreaks,removeIgnores);


% Contrast 0 data
allContrasts=[];
meanContrasts = zeros(1,numDays);
stdContrasts = zeros(1,numDays);
for i=1:numDays
    allContrasts = cat(2,allContrasts,fileStr(i).contrastPC0);
    meanContrasts(i) = mean(fileStr(i).contrastPC0);
    stdContrasts(i)  = std(fileStr(i).contrastPC0);
end

gridPositionContrasts = [0.775 0.25 0.2 0.125];
subplot('Position',gridPositionContrasts);
errorbar(1:numDays,meanContrasts,stdContrasts);
xlabel('Day number'); ylabel('Contrast (%)');
xlim([0 numDays+1]);

gridPositionContrastHistogram = [0.5 0.25 0.2 0.125];
hContrastHist = subplot('Position',gridPositionContrastHistogram);
plotHistogram(hContrastHist,allContrasts,'Contrast (%)');

% Plot performance as a function of target Position
gridPositionPerformanceVsTargetPos = [0.05 0.05 0.2 0.175];
gridPositionPercentTargets         = [0.05 0.25 0.2 0.125];
gridPositionCorrect{1}             = [0.3 0.25 0.15 0.125]; % For target position
gridPositionCorrect{2}             = [0.3 0.05 0.15 0.15]; % For target position

displayTargetInfo(gridPositionPerformanceVsTargetPos,gridPositionPercentTargets,gridPositionCorrect,...
    subjectName,expDates,protocolNames,folderSourceString,gridType,removeBreaks,removeIgnores);

% Fixation window histogram
allFixationData=[];
for i=1:numDays
    folderNameTmp = fullfile(folderSourceString,'data',subjectName,gridType,expDates{i},protocolNames{i});
    [~,fixWindowSizeAllTrials] = getFixationWindowData(folderNameTmp);
    allFixationData=cat(2,allFixationData,fixWindowSizeAllTrials);
end

gridPositionFixationHistogram = [0.5 0.05 0.2 0.15];
hFixationHist = subplot('Position',gridPositionFixationHistogram);
plotHistogram(hFixationHist,allFixationData/2,'Fixation Window (Deg)');

% Get Eye data information (plot mean and std)
gridPositionEyeData = [0.775 0.05 0.2 0.175];
plotEyeMeansAndStd(subjectName,expDates,protocolNames,folderSourceString,gridType,gridPositionEyeData);

    function plotDataLLBar_callBack(~,~)
        
        removeBreaks = get(hBreaks,'val');
        removeIgnores = get(hIgnores,'val');
        
        for ii=1:numDays
            gridPositionLLBar = [LLBarGrid(1)+(ii-1)*dX LLBarGrid(2) dX LLBarGrid(4)];
            displayLLbarGRF(subjectName,expDates{ii},protocolNames{ii},folderSourceString,gridType,gridPositionLLBar,removeBreaks,removeIgnores);
        end

        gridPositionLLBar = [LLBarGrid(1)+(numDays+1)*dX LLBarGrid(2) dX LLBarGrid(4)];
        displayLLbarGRF(subjectName,expDates,protocolNames,folderSourceString,gridType,gridPositionLLBar,removeBreaks,removeIgnores);
        
        num = get(hChooseDay,'val');
        
        if num<=numDays
            histContrast = fileStr(num).contrastPC0;
            displayTargetInfo(gridPositionPerformanceVsTargetPos,gridPositionPercentTargets,gridPositionCorrect,...
            subjectName,expDates{num},protocolNames{num},folderSourceString,gridType,removeBreaks,removeIgnores);
            [~,fixationData] = getFixationWindowData(fullfile(folderSourceString,'data',subjectName,gridType,expDates{num},protocolNames{num}));
        else
            histContrast = allContrasts;
            displayTargetInfo(gridPositionPerformanceVsTargetPos,gridPositionPercentTargets,gridPositionCorrect,...
            subjectName,expDates,protocolNames,folderSourceString,gridType,removeBreaks,removeIgnores);
            fixationData = allFixationData;
        end
        plotHistogram(hContrastHist,histContrast,'Contrast (%)');
        plotHistogram(hFixationHist,fixationData/2,'Fixation Window (Deg)');
    end

    function plotData_callBack(~,~)

        num = get(hChooseDay,'val');
        
        if num<=numDays
            histContrast = fileStr(num).contrastPC0;
            displayTargetInfo(gridPositionPerformanceVsTargetPos,gridPositionPercentTargets,gridPositionCorrect,...
            subjectName,expDates{num},protocolNames{num},folderSourceString,gridType,removeBreaks,removeIgnores);
            [~,fixationData] = getFixationWindowData(fullfile(folderSourceString,'data',subjectName,gridType,expDates{num},protocolNames{num}));
        else
            histContrast = allContrasts;
            displayTargetInfo(gridPositionPerformanceVsTargetPos,gridPositionPercentTargets,gridPositionCorrect,...
            subjectName,expDates,protocolNames,folderSourceString,gridType,removeBreaks,removeIgnores);
            fixationData = allFixationData;
        end
        plotHistogram(hContrastHist,histContrast,'Contrast (%)');
        plotHistogram(hFixationHist,fixationData/2,'Fixation Window (Deg)');
    end
end

function displayLLbarGRF(subjectName,expDates,protocolNames,folderSourceString,gridType,gridPositionLLBar,removeBreaks,removeIgnores)

%Fonts
fontSizeLarge=16;
textPosX = 0.5; textPosY = 1.2;

if iscell(expDates)
    allEOTCodes=[];
    for i=1:length(expDates)
        clear expDate protocolName
        expDate=expDates{i};
        protocolName=protocolNames{i};
        load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData','LL.mat'));
        allEOTCodes = cat(2,allEOTCodes,LL.eotCode);
    end
else
    expDate=expDates;
    protocolName=protocolNames;
    load(fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData','LL.mat')); % Single day
    allEOTCodes = [LL.eotCode];
end

[~,~,numTrials] = displayEOTCodes(gridPositionLLBar,allEOTCodes,protocolName,1,[],removeBreaks,removeIgnores);
 text(textPosX,textPosY,num2str(numTrials), ...
     'Units','Normalized','Rotation',90,'FontSize',fontSizeLarge,'HorizontalAlignment','center','color','k');
end
function dayString=getDayString(subjectName,expDates,protocolNames)

dayString = '';
for i=1:length(expDates)
    dayString = cat(2,dayString,[num2str(i) ' ' subjectName expDates{i} protocolNames{i} '|']);
end
dayString = [dayString 'all days'];
end
function [fixWindowSizeEyeStims,fixWindowSizeGoodStims] = getFixationWindowData(folderNameMain)
% fixWindowSizeEyeStims - the fixation window size for each stimulus for
% which we analyze the energy of the eye Data.
% fixWindowSizeGoodStims - fixation window size for all good stimuli
% Get the fixationWindowSizes
load(fullfile(folderNameMain,'extractedData','BehaviorData.mat')); % returns the variable allTrials

% Get the stimResults
load(fullfile(folderNameMain,'extractedData','stimResults.mat')); % returns the variable stimResults

% Get goodStimNums
load(fullfile(folderNameMain,'extractedData','goodStimNums.mat')); % returns the variable goodStimNums

if exist(fullfile(folderNameMain,'extractedData','validStimAfterTarget.mat'),'file')
    load (fullfile(folderNameMain,'extractedData','validStimAfterTarget.mat'));
    goodStimNums(validStimuliAfterTarget)=[];  % Get rid of these stimuli
end

stimPos = stimResults.stimPosition(goodStimNums);
fixWindowSizeGoodStims = allTrials.fixWindowSize(stimResults.trialNumber(goodStimNums));

fixWindowSizeEyeStims = fixWindowSizeGoodStims(stimPos>1);
end
function plotHistogram(hHist,histContrast,xLabelStr)
axes(hHist); %#ok<*MAXES>
cla;
hist(hHist,histContrast,100);
text(0.4,0.9,[num2str(mean(histContrast)) ' +- ' num2str(std(histContrast))],'Units','Normalized');
xlabel(xLabelStr);
axis tight
end
function plotEyeMeansAndStd(subjectName,expDates,protocolNames,folderSourceString,gridType,gridPositionEyeData)
numDays = length(expDates);
plotHandles = getPlotHandles(2,1,gridPositionEyeData,0,0.025);

mX   = zeros(1,numDays); mY   = zeros(1,numDays);
stdX = zeros(1,numDays); stdY = zeros(1,numDays);

for i=1:numDays
    expDate = expDates{i};
    protocolName = protocolNames{i};
    
    folderNameMain = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    inputFolder  = fullfile(folderNameMain,'segmentedData','eyeData');
    load(fullfile(inputFolder,'eyeDataDeg.mat'));
    
    if iscell(eyeDataDegX)
        [eyeDataDegX,eyeDataDegY] = getValuesFromCellArray(eyeDataDegX,eyeDataDegY);
    end
    mX(i) = mean(eyeDataDegX(:)); stdX(i) = std(eyeDataDegX(:)); %#ok<*NODEF>
    mY(i) = mean(eyeDataDegY(:)); stdY(i) = std(eyeDataDegY(:));
    
    [~,fixationData] = getFixationWindowData(folderNameMain);
    plot(plotHandles(1,:),i,unique(fixationData)/2,'r+'); hold(plotHandles(1,:),'on');
    plot(plotHandles(1,:),i,-unique(fixationData)/2,'r+');
    
    plot(plotHandles(2,:),i,unique(fixationData)/2,'r+'); hold(plotHandles(2,:),'on');
    plot(plotHandles(2,:),i,-unique(fixationData)/2,'r+');
end

errorbar(plotHandles(1,:),mX,stdX); 
axis(plotHandles(1,:),[0 numDays+1 -1*max(fixationData) max(fixationData)]);

errorbar(plotHandles(2,:),mY,stdY); hold(plotHandles(2,:),'on');
axis(plotHandles(2,:),[0 numDays+1 -1*max(fixationData) max(fixationData)]);

end
function fileStr = getFileStr(subjectName,expDates,protocolNames,folderSourceString,gridType)

for i=1:length(expDates)
    load(fullfile(folderSourceString,'data',subjectName,gridType,expDates{i},protocolNames{i},'extractedData','LL.mat'));
    
    if sum(LL.stimType0)==0  % Task Gabor is always hidden
        fileStr(i).azimuth = '0';
        fileStr(i).elevation = '0';
        fileStr(i).contrastPC0 = 0*(LL.contrastPC0);
        fileStr(i).orientationDeg0 = 0;
        fileStr(i).sigmaDeg0 = 0;
        fileStr(i).spatialFreqCPD0 = 0;
    else
        if length(unique(LL.azimuthDeg0))>1
            disp('More than one azimuths...');
            fileStr(i).azimuth = 'many'; %#ok<*AGROW>
        else
            fileStr(i).azimuth     = num2str(unique(LL.azimuthDeg0));
        end
        
        if length(unique(LL.elevationDeg0))>1
            disp('More than one elevations...');
            fileStr(i).elevation = 'many';
        else
            fileStr(i).elevation     = num2str(unique(LL.elevationDeg0));
        end
        
        fileStr(i).contrastPC0 = LL.contrastPC0;
        fileStr(i).orientationDeg0 = min(LL.orientationDeg0); % Two values, x and x+90 (for targets)
        fileStr(i).sigmaDeg0 = unique(LL.sigmaDeg0);
        fileStr(i).spatialFreqCPD0 = unique(LL.spatialFreqCPD0);
    end
end
end
function [expDates,protocolNames] = convertToCell(expDate,protocolName)
expDates{1} = expDate;
protocolNames{1} = protocolName;
end
function [Xout,Yout] = getValuesFromCellArray(X,Y)
Xout=[];
Yout=[];

for i=1:length(X)
   Xout = cat(1,Xout,X{i}(:));
   Yout = cat(1,Yout,Y{i}(:));
end
end