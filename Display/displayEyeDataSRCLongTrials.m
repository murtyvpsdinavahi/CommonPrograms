% The program shows the eye data as a function of different stimulus
% conditions

function displayEyeDataSRCLongTrials(monkeyName,expDates,protocolNames,folderSourceString,gridType)

%%%% display properties
fontSizeSmall=10; fontSizeMedium=12; fontSizeLarge=14;
backgroundColor = 'w';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots 
hEyePlots        = getPlotHandles(1,2,[0.05 0.60 0.4 0.3],0.05,0);
hEyeSigPlots     = getPlotHandles(1,2,[0.05 0.50 0.4 0.1],0.05,0);

hEyeMeanPlots    = getPlotHandles(2,1,[0.5 0.50 0.10 0.4],0.05,0);
hEyeMeanSigPlots = getPlotHandles(2,1,[0.6 0.50 0.10 0.4],0.05,0);
hEyePlotVsStimPos = getPlotHandles(2,1,[0.75 0.50 0.20 0.4],0.05,0);

hEyeMeanPlotsSingleDay = getPlotHandles(2,1,[0.75 0.05 0.2 0.4],0.05,0);

hMSPlot         = subplot('position',[0.275 0.05 0.175 0.4]);
%hMSSigPlot      = subplot('position',[0.05 0.05 0.175 0.1]);
hSpeedHistogram  = subplot('position',[0.05 0.05 0.175 0.4]);

hMSMeanPlot      = subplot('position',[0.5 0.05 0.10 0.4]);
hMSMeanSigPlot   = subplot('position',[0.6 0.05 0.10 0.4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls

hChoicePanel = uipanel('Title','Day','Unit','Normalized','fontSize', fontSizeLarge,'Position',[0 0.925 0.15 0.075]);

% Make protocolString
protocolString = getProtocolString(expDates,protocolNames);
uicontrol('Parent',hChoicePanel,'Unit','Normalized', 'Position',[0 0.5 0.5 0.5], ...
    'Style','text','String','Day','FontSize',fontSizeSmall);
hDay = uicontrol('Parent',hChoicePanel,'Unit','Normalized', 'Position',[0.5 0.5 0.5 0.5], ...
    'BackgroundColor', backgroundColor, ...
    'Style','popup','String',protocolString,'FontSize',fontSizeSmall);

% Speed cut off ranges
uicontrol('Parent',hChoicePanel,'Unit','Normalized', 'Position',[0 0 0.5 0.5], ...
    'Style','text','String','Cutoff(Deg/s):','FontSize',fontSizeSmall);
hCutoff = uicontrol('Parent',hChoicePanel,'Unit','Normalized', 'Position',[0.5 0 0.5 0.5], ...
    'BackgroundColor', backgroundColor, ...
    'Style','edit','String','25','FontSize',fontSizeSmall);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot, Rescale, Hold and cla
hPlottingPanel = uipanel('Title','Plot options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.15 0.925 0.15 0.075]);

uicontrol('Parent',hPlottingPanel,'Unit','Normalized', 'Position',[0 0.5 0.5 0.5], ...
    'Style','pushbutton','String','Plot','FontSize',fontSizeMedium, ...
    'Callback',{@plot_Callback});
% uicontrol('Parent',hPlottingPanel,'Unit','Normalized', 'Position',[0.5 0.5 0.5 0.5], ...
%     'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium, ...
%     'Callback',{@rescale_Callback});
% uicontrol('Parent',hPlottingPanel,'Unit','Normalized','Position',[0 0 0.5 0.5], ...
%     'Style','pushbutton','String','plot size','FontSize',fontSizeMedium, ...
%     'Callback',{@replotSizePlot_Callback});
uicontrol('Parent',hPlottingPanel,'Unit','Normalized', 'Position',[0.5 0.5 0.5 0.5], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

uicontrol('Parent',hPlottingPanel,'Unit','Normalized', 'Position',[0 0 0.5 0.5], ...
    'Style','text','String','Test Method','FontSize',fontSizeSmall);

[testMethodString,testMethods] = getTestMethods;
hTestMethod = uicontrol('Parent',hPlottingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position', [0.5 0 0.5 0.5], ...
    'Style','popup','String',testMethodString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the parameters. They should be the same for all the days
folderExtract = fullfile(folderSourceString,'data',monkeyName,gridType,expDates{1},protocolNames{1},'extractedData');
[~,cValsUnique,tValsUnique,eValsUnique,~,sValsUnique] = loadParameterCombinations(folderExtract);

% Get properties of the Stimulus
stimResults = loadStimResults(folderExtract);

parameterTextWidth = 0.25; parameterWidth = 0.25;
hParameterPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[0.3 0.925 0.4 0.075]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eotCode
EOTCodeString = getEOTCodeString(eValsUnique);
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 0.5 parameterTextWidth 0.5],...
    'Style','text','String','EOTCode','FontSize',fontSizeSmall);
hEOTCode = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position', [parameterTextWidth 0.5 parameterWidth 0.5], ...
    'Style','popup','String',EOTCodeString,'FontSize',fontSizeSmall);

% stimType
stimTypeString = getStimTypeString(sValsUnique);
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[0 0 parameterTextWidth 0.5], ...
    'Style','text','String','Stim Type','FontSize',fontSizeSmall);
hStimType = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position', [parameterTextWidth 0 parameterWidth 0.5], ...
    'Style','popup','String',stimTypeString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contrast
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[parameterTextWidth+parameterWidth 0.5 parameterTextWidth 0.5], ...
    'Style','text','String','Contrast (%)','FontSize',fontSizeSmall);

contrastString = getContrastString(cValsUnique,stimResults);
hContrast= uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position', [parameterTextWidth+parameterWidth+parameterTextWidth 0.5 parameterWidth 0.5], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeSmall);

% Temporal Frequency
uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'Position',[parameterTextWidth+parameterWidth 0 parameterTextWidth 0.5], ...
    'Style','text','String','TF (Hz)','FontSize',fontSizeSmall);

temporalFreqString = getTemporalFreqString(tValsUnique,stimResults);
hTemporalFreq = uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position', [parameterTextWidth+parameterWidth+parameterTextWidth 0 parameterWidth 0.5], ...
    'Style','popup','String',temporalFreqString,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Timing controls
hTimingPanel = uipanel('Title','Timing','Unit','Normalized','fontSize', fontSizeLarge, 'Position',[0.7 0.925 0.3 0.075]);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 0.5 0.25 0.5], ...
    'Style','text','String','Time Period','FontSize',fontSizeSmall);

hTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.25 0.5 0.25 0.5], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);

hTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [0.5 0.5 0.25 0.5], ...
    'Style','edit','String','0.4','FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    function cla_Callback(~,~)
        
        claGivenPlotHandle(hEyePlots); claGivenPlotHandle(hEyeSigPlots);
        claGivenPlotHandle(hEyePlotVsStimPos);
        claGivenPlotHandle(hEyeMeanPlots); claGivenPlotHandle(hEyeMeanSigPlots);
        claGivenPlotHandle(hEyeMeanPlotsSingleDay);
        
        claGivenPlotHandle(hMSMeanPlot); claGivenPlotHandle(hMSMeanSigPlot);
        claGivenPlotHandle(hMSPlot); claGivenPlotHandle(hSpeedHistogram);
    end
    function plot_Callback(~,~)
        
        claGivenPlotHandle(hEyePlots); claGivenPlotHandle(hEyeSigPlots);
        claGivenPlotHandle(hEyePlotVsStimPos);
        claGivenPlotHandle(hEyeMeanPlots); claGivenPlotHandle(hEyeMeanSigPlots);
        
        claGivenPlotHandle(hMSMeanPlot); claGivenPlotHandle(hMSMeanSigPlot);
        claGivenPlotHandle(hMSPlot); claGivenPlotHandle(hSpeedHistogram);
        %claGivenPlotHandle(hEyeMeanPlotsSingleDay);
        
        dayNum = get(hDay,'val');
        cutoff = str2double(get(hCutoff,'string'));
        useThisTestMethod = testMethods{get(hTestMethod,'val')};   
        useThisTimeRange=[str2double(get(hTMin,'string')) str2double(get(hTMax,'string'))];
        
        %Initialze
        clear expDate protocolName
        if dayNum > length(expDates)
            expDate = expDates;
            protocolName = protocolNames;
        else
            expDate = expDates{dayNum};
            protocolName = protocolNames{dayNum};
        end
        
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        e=get(hEOTCode,'val');
        s=get(hStimType,'val');
        
        clear eyeX eyeY
        colorNames = jet(2);
        [eyeX,eyeY,xs,MSData,allEyeSpeeds,numMSInRange] = getSortedEyeData(monkeyName,expDate,protocolName,folderSourceString,gridType,c,t,e,s,cutoff,useThisTimeRange);

        % Display Eye Data into categories
        compareAndDisplayData(eyeX,xs,hEyePlots(1),hEyeSigPlots(1),useThisTestMethod,colorNames);%,colorX,colorSig,displayPlot,showStdErr);
        compareAndDisplayData(eyeY,xs,hEyePlots(2),hEyeSigPlots(2),useThisTestMethod,colorNames);
        
        % Put labels
        ylabel(hEyePlots(1),'x eyedata (deg)'); ylabel(hEyePlots(2),'y eyedata (deg)');
        xlabel(hEyeSigPlots(1),'time (s)'); ylabel(hEyeSigPlots(1),'p-value');       
        xlabel(hEyeSigPlots(2),'time (s)'); ylabel(hEyeSigPlots(2),'p-value'); 
        
        %put legend
        legendStr1{1}= 'out'; legendStr1{2}= 'in';
        legendStr2 = cell(1,2);
        for legendPos=1:2
            legendStr2{legendPos} =  ['n=' num2str(size(eyeX{legendPos},1))];
        end
        legend(hEyePlots(1),legendStr1,'Location','NorthOutside');
        legend(hEyePlots(2),legendStr2,'Location','NorthOutside');
        
        % Significance analysis on Mean Data
        compareAndDisplayMeanData(eyeX,xs,hEyeMeanPlots(1),hEyeMeanSigPlots(1),useThisTestMethod,useThisTimeRange,colorNames);
        compareAndDisplayMeanData(eyeY,xs,hEyeMeanPlots(2),hEyeMeanSigPlots(2),useThisTestMethod,useThisTimeRange,colorNames);
        ylabel(hEyeMeanPlots(1),'x position (deg)'); ylabel(hEyeMeanPlots(2),'y position (deg)');
        
        % Plot eye position versus stimulus Pos
        %plotEyePosititionVsStimPos(folderSourceString,monkeyName,expDate,protocolName,hEyePlotVsStimPos(1),hEyePlotVsStimPos(2));

        % Show Microsaccades      
        % Histogram
        eyeSpeedCenters=0:100;
        histAllSpeeds = hist(allEyeSpeeds,eyeSpeedCenters);
        plot(hSpeedHistogram,eyeSpeedCenters,log10(histAllSpeeds),'k');
        hold(hSpeedHistogram,'on');
        MSPos = find(eyeSpeedCenters>cutoff);
        plot(hSpeedHistogram,eyeSpeedCenters(MSPos),log10(histAllSpeeds(MSPos)),'r');
        hold(hSpeedHistogram,'off');
        ylabel(hSpeedHistogram,'log(#)');
        xlabel(hSpeedHistogram,'speed (deg/s)');
        
 
        % # microsaccades /sec
        for ii=1:2
            [H,timeValsMS] = getPSTH(MSData{ii},20,[xs(1) xs(end)]);
            plot(hMSPlot,timeValsMS,H,'color',colorNames(ii,:));
            hold(hMSPlot,'on');
        end
        hold(hMSPlot,'off');
        axis(hMSPlot,'tight');
        xlabel(hMSPlot,'time (s)'); ylabel(hMSPlot,'MicroSaccades/s');
        
        % Significance test
        compareAndDisplayMeanData2(numMSInRange,hMSMeanPlot,hMSMeanSigPlot,useThisTestMethod,colorNames);
        
         % if expDates is a cell array (population data, report the p-values
        % of individual days as well.
        
        if iscell(expDate)
            K=length(expDate);
            pXAllDays=zeros(1,K);
            pYAllDays=zeros(1,K);
            for ii=1:K
                clear eyeX eyeY xs aziStr eleStr comparisonDataValues
                disp(['Performing significance analysis on ' num2str(ii) ': ' expDate{ii} protocolName{ii}]);
                [eyeX,eyeY,xs] = getSortedEyeData(monkeyName,expDate{ii},protocolName{ii},folderSourceString,gridType,c,t,e,s,cutoff,useThisTimeRange);
                pXAllDays(ii)=compareAndDisplayMeanData(eyeX,xs,[],[],useThisTestMethod,useThisTimeRange,[],0);
                pYAllDays(ii)=compareAndDisplayMeanData(eyeY,xs,[],[],useThisTestMethod,useThisTimeRange,[],0);
            end
           
            plotSignificanceData(pXAllDays,1:K,hEyeMeanPlotsSingleDay(1),'k',0.05);
            plotSignificanceData(pYAllDays,1:K,hEyeMeanPlotsSingleDay(2),'k',0.05);
            
            ylabel(hEyeMeanPlotsSingleDay(1),'p-value X data');
            ylabel(hEyeMeanPlotsSingleDay(2),'p-value Y data');
            xlabel(hEyeMeanPlotsSingleDay(2),'day number');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function protocolString = getProtocolString(expDates,protocolNames)
protocolString='';
for i=1:length(expDates)
    expDate = expDates{i};
    protocolName = protocolNames{i};

    protocolString = cat(2,protocolString,[num2str(i) ' ' expDate protocolName '|']);
end
protocolString = [protocolString 'all Days'];
end
function contrastString = getContrastString(cIndexValsUnique,stimResults)
if isfield(stimResults,'contrast0PC')
    [cVals0Unique,cVals1Unique] = getValsFromIndex(cIndexValsUnique,stimResults,'contrast');
    if length(cVals0Unique)==1
        contrastString = [num2str(cVals0Unique) ',' num2str(cVals1Unique)];
    else
        contrastString = '';
        for i=1:length(cVals0Unique)
            contrastString = cat(2,contrastString,[num2str(cVals0Unique(i)) ',' num2str(cVals1Unique(i)) '|']);
        end
        contrastString = [contrastString 'all'];
    end
    
else % Work with indices
    if length(cIndexValsUnique)==1
        if cIndexValsUnique ==0
            contrastString = '0';
        else
            contrastString = num2str(100/2^(7-cIndexValsUnique));
        end
        
    else
        contrastString = '';
        for i=1:length(cIndexValsUnique)
            if cIndexValsUnique(i) == 0
                contrastString = cat(2,contrastString,[ '0|']); %#ok<*NBRAK>
            else
                contrastString = cat(2,contrastString,[num2str(100/2^(7-cIndexValsUnique(i))) '|']);
            end
        end
        contrastString = [contrastString 'all'];
    end
end
end
function temporalFreqString = getTemporalFreqString(tIndexValsUnique,stimResults)

if isfield(stimResults,'temporalFreq0Hz')
    [tVals0Unique,tVals1Unique] = getValsFromIndex(tIndexValsUnique,stimResults,'temporalFreq');
    if length(tIndexValsUnique)==1
        temporalFreqString = [num2str(tVals0Unique) ',' num2str(tVals1Unique)];
    else
        temporalFreqString = '';
        for i=1:length(tIndexValsUnique)
            temporalFreqString = cat(2,temporalFreqString,[num2str(tVals0Unique(i)) ',' num2str(tVals1Unique(i)) '|']);
        end
        temporalFreqString = [temporalFreqString 'all'];
    end
else
    if length(tIndexValsUnique)==1
        if tIndexValsUnique ==0
            temporalFreqString = '0';
        else
            temporalFreqString = num2str(min(50,80/2^(7-tIndexValsUnique)));
        end
        
    else
        temporalFreqString = '';
        for i=1:length(tIndexValsUnique)
            if tIndexValsUnique(i) == 0
                temporalFreqString = cat(2,temporalFreqString,[ '0|']);
            else
                temporalFreqString = cat(2,temporalFreqString,[num2str(min(50,80/2^(7-tIndexValsUnique(i)))) '|']);
            end
        end
        temporalFreqString = [temporalFreqString 'all'];
    end
end
end
function EOTCodeString = getEOTCodeString(eValsUnique)

if length(eValsUnique)==1
    if eValsUnique == 0
        EOTCodeString = 'Correct';
    elseif eValsUnique == 1
        EOTCodeString = 'Wrong';
    elseif eValsUnique == 2
        EOTCodeString = 'Failed';
    elseif eValsUnique == 3
        EOTCodeString = 'Broke';
    elseif eValsUnique == 4
        EOTCodeString = 'Ignored';
    elseif eValsUnique == 5
        EOTCodeString = 'False Alarm';
    elseif eValsUnique == 6
        EOTCodeString = 'Distracted';
    elseif eValsUnique == 7
        EOTCodeString = 'Force Quit';
    else
        disp('Unknown EOT Code');
    end
else
    EOTCodeString = '';
    for i=1:length(eValsUnique)
        if eValsUnique(i) == 0
            EOTCodeString = [EOTCodeString 'Correct|']; %#ok<*AGROW>
        elseif eValsUnique(i) == 1
            EOTCodeString = [EOTCodeString 'Wrong|'];
        elseif eValsUnique(i) == 2
            EOTCodeString = [EOTCodeString 'Failed|'];
        elseif eValsUnique(i) == 3
            EOTCodeString = [EOTCodeString 'Broke|'];
        elseif eValsUnique(i) == 4
            EOTCodeString = [EOTCodeString 'Ignored|'];
        elseif eValsUnique(i) == 5
            EOTCodeString = [EOTCodeString 'False Alarm|'];
        elseif eValsUnique(i) == 6
            EOTCodeString = [EOTCodeString 'Distracted|'];
        elseif eValsUnique(i) == 7
            EOTCodeString = [EOTCodeString 'Force Quit|'];
        else
            disp('Unknown EOT Code');
        end
    end
    EOTCodeString = [EOTCodeString 'all'];
end
end

function stimTypeString = getStimTypeString(sValsUnique)

if length(sValsUnique)==1
    if sValsUnique == 0
        stimTypeString = 'Null';
    elseif sValsUnique == 1
        stimTypeString = 'Correct';
    elseif sValsUnique == 2
        stimTypeString = 'Target';
    elseif sValsUnique == 3
        stimTypeString = 'FrontPad';
    elseif sValsUnique == 4
        stimTypeString = 'BackPad';
    else
        disp('Unknown Stimulus Type');
    end
else
    stimTypeString = '';
    for i=1:length(sValsUnique)
        if sValsUnique(i) == 0
            stimTypeString = [stimTypeString 'Null|'];
        elseif sValsUnique(i) == 1
            stimTypeString = [stimTypeString 'Correct|'];
        elseif sValsUnique(i) == 2
            stimTypeString = [stimTypeString 'Target|'];
        elseif sValsUnique(i) == 3
            stimTypeString = [stimTypeString 'FrontPad|'];
        elseif sValsUnique(i) == 4
            stimTypeString = [stimTypeString 'BackPad|'];
        else
            disp('Unknown Stimulus Type');
        end
    end
    stimTypeString = [stimTypeString 'all'];
end
end
function claGivenPlotHandle(plotHandles)
[numRows,numCols] = size(plotHandles);
for ii=1:numRows
    for jj=1:numCols
        cla(plotHandles(ii,jj));
    end
end
end

function [parameterCombinations,cValsUnique,tValsUnique,eValsUnique,...
    aValsUnique,sValsUnique] = loadParameterCombinations(folderExtract) %#ok<*STOUT>

load(fullfile(folderExtract,'parameterCombinations.mat'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eyeX,eyeY,timeValsEyePos,MSData,allEyeSpeeds,numMSInRange] = getSortedEyeData(monkeyName,expDate,protocolName,folderSourceString,gridType,c,t,e,s,cutoff,timeRange)

FsEye=200;
numberOfCategories=2;

if ~iscell(expDate)
    expDates{1} = expDate;
    protocolNames{1} = protocolName;
else
    expDates=expDate;
    protocolNames=protocolName;
end

% Initialization
eyeX=cell(1,numberOfCategories);
eyeY=cell(1,numberOfCategories);
MSData=cell(1,numberOfCategories);
numMSInRange=cell(1,numberOfCategories);

allEyeSpeeds=[];

for i=1:length(expDates)
    
    folderExtract = fullfile(folderSourceString,'data',monkeyName,gridType,expDates{i},protocolNames{i},'extractedData');
    folderSegment = fullfile(folderSourceString,'data',monkeyName,gridType,expDates{i},protocolNames{i},'segmentedData');
    
    % get timevals for eye position. 
    if i==1
        load(fullfile(folderExtract,'EyeData.mat')); % returns the variable 'timeVals'
        if s==1 % valid
            timeValsEyePos   = (eyeRangeMS(1):1000/FsEye:eyeRangeMS(2)-1000/FsEye)/1000;
        elseif s==2
            timeValsEyePos = (0:1000/FsEye:eyeRangeMS(2)-1000/FsEye)/1000;
        end
    end

    % eyeData
    load(fullfile(folderSegment,'eyeData','eyeDataDeg'));
    
    % eyeSpeed
    load(fullfile(folderSegment,'eyeData','eyeSpeed'));
    
    % parameterCombinations
    parameterCombinations = loadParameterCombinations(folderExtract);

    numberOfCategories=2;
    
    for j=1:numberOfCategories
        eyeIndices = parameterCombinations{c,t,e,j,s};
        eyeX{j} = cat(1,eyeX{j},cell2Array(eyeDataDegX(eyeIndices),1)');
        eyeY{j} = cat(1,eyeY{j},cell2Array(eyeDataDegY(eyeIndices),1)');
        
        clear eyeSpeedMag speedXtmp speedYtmp
        
        speedXtmp = cell2Array(eyeSpeedX(eyeIndices),1)';
        speedYtmp = cell2Array(eyeSpeedY(eyeIndices),1)';
        eyeSpeedMag = FsEye*sqrt(speedXtmp.^2 + speedYtmp.^2);
        
        clear MStmp numMStmp
        [MStmp,numMStmp] = findMicroSaccades(eyeSpeedMag,cutoff,timeValsEyePos,timeRange);
        MSData{j} = [MSData{j} MStmp];
        numMSInRange{j}  = [numMSInRange{j} numMStmp];
        
        clear tmpEyeSpeeds
        tmpEyeSpeeds = eyeSpeedMag;
        allEyeSpeeds=cat(1,allEyeSpeeds,tmpEyeSpeeds(:));
    end
end
end
function Y = cell2Array(X,catRows)

if ~exist('catRows','var');          catRows=0;                         end

Y=[];
for i=1:length(X)
    
    if catRows
        Y=cat(2,Y,X{i}(:));
    else
        Y=cat(1,Y,X{i}(:));
    end
end
end
function [testMethodString,testMethods] = getTestMethods

testMethods{1} = 'default';
testMethods{2} = 'anova';
testMethods{3} = 'kruskalWallis';
testMethods{4} = 'ttest';

testMethodString='';
for i=1:length(testMethods)
    testMethodString = cat(2,testMethodString,[testMethods{i} '|']);
end
testMethodString = removeIfPresent(testMethodString,'|');
end
function stimResults = loadStimResults(folderExtract)
load (fullfile(folderExtract,'stimResults'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [valList0Unique,valList1Unique] = getValsFromIndex(indexListUnique,stimResults,fieldName)
if isfield(stimResults,[fieldName 'Index'])
    
    indexList = getfield(stimResults,[fieldName 'Index']); %#ok<*GFLD>
    if strcmpi(fieldName,'contrast')
        valList0 = getfield(stimResults,[fieldName '0PC']);
        valList1 = getfield(stimResults,[fieldName '1PC']);
    else
        valList0 = getfield(stimResults,[fieldName '0Hz']);
        valList1 = getfield(stimResults,[fieldName '1Hz']);
    end
    
    numList = length(indexListUnique);
    valList0Unique = zeros(1,numList);
    valList1Unique = zeros(1,numList);
    for i=1:numList
        valList0Unique(i) = unique(valList0(indexListUnique(i)==indexList));
        valList1Unique(i) = unique(valList1(indexListUnique(i)==indexList));
    end
else
    valList0Unique = indexListUnique;
    valList1Unique = indexListUnique;
end
end