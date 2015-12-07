% Display All Channels
function displayAllChannelsSRC(subjectName,expDate,protocolName,folderSourceString,gridType)

if ~exist('folderSourceString','var');   folderSourceString='E:';       end
if ~exist('gridType','var');             gridType='EEG';                end

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');

% load LFP Information
[analogChannelsStored,timeVals] = loadlfpInfo(folderLFP);
[neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes);

% Get Combinations
[parameterCombinations,cIndexValsUnique,tIndexValsUnique,eIndexValsUnique,...
    aIndexValsUnique,sIndexValsUnique] = loadParameterCombinations(folderExtract);

% Get properties of the Stimulus
stimResults = loadStimResults(folderExtract);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.25; panelStartHeight = 0.725;
staticPanelWidth = 0.25; staticStartPos = 0.025;
dynamicPanelWidth = 0.25; dynamicStartPos = 0.275;
timingPanelWidth = 0.25; timingStartPos = 0.525;
plotOptionsPanelWidth = 0.2; plotOptionsStartPos = 0.775;
backgroundColor = 'w';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Static Panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% staticTitle = [subjectName '_' expDate '_' protocolName];
% hStaticPanel = uipanel('Title','Information','fontSize', fontSizeLarge, ...
%    'Unit','Normalized','Position',[staticStartPos panelStartHeight staticPanelWidth panelHeight]);
% 
% if isfield(stimResults,'orientation')
%     oriStr = num2str(stimResults.orientation);
% else
%     oriStr = [num2str(stimResults.baseOrientation0) ',' num2str(stimResults.baseOrientation1)];
% end
% 
% staticText = [{ '   '};
%     {['Monkey Name: ' subjectName]}; ...
%     {['Date: ' expDate]}; ...
%     {['Protocol Name: ' protocolName]}; ...
%     {'   '}
%     {['Orientation  (Deg): ' oriStr]}; ...
%     {['Spatial Freq (CPD): ' num2str(stimResults.spatialFrequency)]}; ...
%     {['Azimuth      (Deg): ' num2str(stimResults.azimuth)]}; ...
%     {['Elevation    (Deg): ' num2str(stimResults.elevation)]}; ...
%     {['Sigma        (Deg): ' num2str(stimResults.sigma)]}; ...
%     {['Radius       (Deg): ' num2str(stimResults.radius)]}; ...
%     ];
% 
% disp(staticText);
% 
% uicontrol('Parent',hStaticPanel,'Unit','Normalized', ...
%     'Position',[0 0 1 1], 'Style','text','String',staticText,'FontSize',fontSizeSmall);

electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
hElectrodes = showElectrodeLocations(electrodeGridPos,[],[],[],1,0,gridType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.08; dynamicGap=0.02; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Contrast
contrastString = getContrastString(cIndexValsUnique,stimResults);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Contrast (%)','FontSize',fontSizeSmall);
hContrast = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeSmall);

% Temporal Frequency
temporalFreqString = getTemporalFreqString(tIndexValsUnique,stimResults);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeSmall);
hTemporalFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',temporalFreqString,'FontSize',fontSizeSmall);

% EOT Codes
EOTCodeString = getEOTCodeString(eIndexValsUnique);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','EOT Code','FontSize',fontSizeSmall);
hEOTCode = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-3*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',EOTCodeString,'FontSize',fontSizeSmall);

% Attend Loc
attendLocString = getAttendLocString(aIndexValsUnique);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Attended Location','FontSize',fontSizeSmall);
hAttendLoc = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-4*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',attendLocString,'FontSize',fontSizeSmall);

% Stimulus Type
stimTypeString = getStimTypeString(sIndexValsUnique);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Stimulus Type','FontSize',fontSizeSmall);
hStimType = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',stimTypeString,'FontSize',fontSizeSmall);

% Analysis Type
analysisTypeString = 'ERP|Firing Rate|FFT|delta FFT';
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-6*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

signalRange = [timeVals(1) timeVals(end)];
fftRange = [0 250];
baseline = [-0.2 0];
stimPeriod = [0 0.2];

% Signal Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Parameter','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Min','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeMedium);

% Stim Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);

% FFT Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','FFT Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Y Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-8*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Y Range','FontSize',fontSizeSmall);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','1','FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.1;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);

% Button for Plotting
[colorString, colorNames] = getColorString;
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Color','FontSize',fontSizeSmall);
hChooseColor = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',colorString,'FontSize',fontSizeSmall);

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 4*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale Y','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale X','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get main plot and message handles
if strcmpi(gridType,'ECoG')
    plotHandles = getPlotHandles(8,10,[0.05 0.05 0.9 0.625]);
elseif strcmpi(gridType,'EEG')
    plotHandles = getPlotHandles(9,11,[0.05 0.05 0.9 0.625]);
else
    plotHandles = getPlotHandles;
end

hMessage = uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],...
    'Style','text','String','DisplayAllChannelsSRC','FontSize',fontSizeLarge);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        e=get(hEOTCode,'val');
        a=get(hAttendLoc,'val');
        s=get(hStimType,'val');
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));

        blRange = [str2double(get(hBaselineMin,'String')) str2double(get(hBaselineMax,'String'))];
        stRange = [str2double(get(hStimPeriodMin,'String')) str2double(get(hStimPeriodMax,'String'))];

        goodPos = parameterCombinations{c,t,e,a,s};
        set(hMessage,'String',[num2str(length(goodPos)) ' stimuli found' ]);

        if ~isempty(goodPos)
            if analysisType == 2
                holdOnState = get(hHoldOn,'val');
                channelsStored = neuralChannelsStored;
                [baselineFiringRate,stimulusFiringRate] = plotSpikeData(plotHandles,channelsStored,goodPos,folderSpikes,...
                    timeVals,plotColor,SourceUnitID,holdOnState,blRange,stRange,gridType);
                
                responsiveElectrodes = (channelsStored(stimulusFiringRate>=5));
                rareElectrodes = setdiff((channelsStored(stimulusFiringRate>=2)),responsiveElectrodes);
                inhibitedElectrodes = (channelsStored(intersect(find(baselineFiringRate>stimulusFiringRate),find(baselineFiringRate>2))));
                disp(['rareUnits: ' num2str(rareElectrodes)]);
                disp(['responsive: ' num2str(responsiveElectrodes)]);
                disp(['inhibited : ' num2str(inhibitedElectrodes)]);
                
                showElectrodeLocations(electrodeGridPos,responsiveElectrodes,'b',hElectrodes,0,0,gridType,subjectName);
                showElectrodeLocations(electrodeGridPos,inhibitedElectrodes,'g',hElectrodes,1,0,gridType,subjectName);
            else
                channelsStored = analogChannelsStored;
                plotLFPData(plotHandles,channelsStored,goodPos,folderLFP,...
                    analysisType,timeVals,plotColor,blRange,stRange,gridType,subjectName);
            end
            
            if analysisType<=2 % ERP or spikes
                xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
            else
                xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
            end
            
            yRange = getYLims(plotHandles,channelsStored,gridType,subjectName);
            rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType == 2
            channelsStored = neuralChannelsStored;
        else
            channelsStored = analogChannelsStored;
        end
        
        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end

        yRange = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];
        rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType == 2
            channelsStored = neuralChannelsStored;
        else
            channelsStored = analogChannelsStored;
        end
        
        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end

        yRange = getYLims(plotHandles,channelsStored,gridType,subjectName);
        rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');

        [numRow,numCol] = size(plotHandles);

        if holdOnState
            for i=1:numRow
                for j=1:numCol
                    set(plotHandles(i,j),'Nextplot','add');
                end
            end
        else
            for i=1:numRow
                for j=1:numCol
                    set(plotHandles(i,j),'Nextplot','replace');
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        [numRow,numCol] = size(plotHandles);
        for i=1:numRow
            for j=1:numCol
                cla(plotHandles(i,j));
            end
        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that plots the data
function plotLFPData(plotHandles, channelsStored, goodPos, ...
    folderData, analysisType, timeVals, plotColor,blRange,stRange,gridType,subjectName)

if isempty(goodPos)
    disp('No entries for this combination..')
else
    if analysisType>2 % FFT
        Fs = round(1/(timeVals(2)-timeVals(1)));
        blPos = find(timeVals>=blRange(1),1)+ (1:diff(blRange)*Fs);
        stPos = find(timeVals>=stRange(1),1)+ (1:diff(stRange)*Fs);

        xsBL = 0:1/(diff(blRange)):Fs-1/(diff(blRange));
        xsST = 0:1/(diff(stRange)):Fs-1/(diff(stRange));
    end

    for i=1:length(channelsStored)
        disp(i)
        channelNum = channelsStored(i);

        % get position
        [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);

        if analysisType == 1        % compute ERP
            clear signal analogData
            load(fullfile(folderData ,['elec' num2str(channelNum)]));
            erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>

            %Plot
            plot(plotHandles(row,column),timeVals,erp,'color',plotColor);

        elseif analysisType == 2    % compute Firing rates
            disp('Use plotSpikeData instead of plotLFPData...');
        else
            clear signal analogData
            load(fullfile(folderData,['elec' num2str(channelNum)]));
            
            fftBL = abs(fft(analogData(goodPos,blPos),[],2)); 
            fftST = abs(fft(analogData(goodPos,stPos),[],2));

            if analysisType == 3
                plot(plotHandles(row,column),xsBL,log10(mean(fftBL)),'g'); 
                hold(plotHandles(row,column),'on');
                plot(plotHandles(row,column),xsST,log10(mean(fftST)),plotColor);
            end

            if analysisType == 4
                if xsBL == xsST %#ok<BDSCI>
                    plot(plotHandles(row,column),xsBL,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
                else
                    disp('Choose same baseline and stimulus periods..');
                end
            end
        end
    end
end
end
function [baselineFiringRate,stimulusFiringRate] = plotSpikeData(plotHandles,channelsStored,goodPos, ...
    folderData, timeVals, plotColor, SourceUnitID,holdOnState,blRange,stRange,gridType)

unitColors = ['r','m','y','c','g'];
binWidthMS = 10;

if isempty(goodPos)
    disp('No entries for this combination..')
else
    numChannels = length(channelsStored);
    baselineFiringRate = zeros(1,numChannels);
    stimulusFiringRate = zeros(1,numChannels);
    
    for i=1:numChannels
        %disp(i)
        channelNum = channelsStored(i);

        % get position
        [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);

        clear neuralInfo spikeData
        load([folderData 'elec' num2str(channelNum) '_SID' num2str(SourceUnitID(i))]);
        [psthVals,xs] = psth_SR(spikeData(goodPos),binWidthMS,timeVals(1),timeVals(end));
        
        % Compute the mean firing rates
        blPos = find(xs>=blRange(1),1)+ (1:(diff(blRange))/(binWidthMS/1000));
        stPos = find(xs>=stRange(1),1)+ (1:(diff(stRange))/(binWidthMS/1000));
        
        baselineFiringRate(i) = mean(psthVals(blPos));
        stimulusFiringRate(i) = mean(psthVals(stPos));
        
        if SourceUnitID(i)==0
            plot(plotHandles(row,column),xs,psthVals,'color',plotColor);
        elseif SourceUnitID(i)> 5
            disp('Only plotting upto 6 single units per electrode...')
        else
            plot(plotHandles(row,column),xs,smooth(psthVals),'color',unitColors(SourceUnitID(i)));
        end
        
        if (i<length(channelsStored))
            if channelsStored(i) == channelsStored(i+1)
                disp('hold on...')
                set(plotHandles(row,column),'Nextplot','add');
            else
                if ~holdOnState
                    set(plotHandles(row,column),'Nextplot','replace');
                end
            end
        end
    end
end
end   
function yRange = getYLims(plotHandles,channelsStored,gridType,subjectName)

% Initialize
yMin = inf;
yMax = -inf;

for i=1:length(channelsStored)
    channelNum = channelsStored(i);
    % get position
    [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);
    
    axis(plotHandles(row,column),'tight');
    tmpAxisVals = axis(plotHandles(row,column));
    if tmpAxisVals(3) < yMin
        yMin = tmpAxisVals(3);
    end
    if tmpAxisVals(4) > yMax
        yMax = tmpAxisVals(4);
    end
end
yRange = [yMin yMax];
end
function rescaleData(plotHandles,channelsStored,axisLims,gridType,subjectName)

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:length(channelsStored)
    channelNum = channelsStored(i);
    % get position
    [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);
    
    axis(plotHandles(row,column),axisLims);
    if (row==numRows && rem(column,2)==1)
        if column~=1
            set(plotHandles(row,column),'YTickLabel',[],'fontSize',labelSize);
        end
    elseif (rem(row,2)==0 && column==1)
        set(plotHandles(row,column),'XTickLabel',[],'fontSize',labelSize);
    else
        set(plotHandles(row,column),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
    end
end

% Remove Labels on the four corners
set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function EOTCodeString = getEOTCodeString(eIndexValsUnique)

if length(eIndexValsUnique)==1
    if eIndexValsUnique == 0
        EOTCodeString = 'Correct';
    elseif eIndexValsUnique == 1
        EOTCodeString = 'Wrong';
    elseif eIndexValsUnique == 2
        EOTCodeString = 'Failed';
    elseif eIndexValsUnique == 3
        EOTCodeString = 'Broke';
    elseif eIndexValsUnique == 4
        EOTCodeString = 'Ignored';
    elseif eIndexValsUnique == 5
        EOTCodeString = 'False Alarm';
    elseif eIndexValsUnique == 6
        EOTCodeString = 'Distracted';
    elseif eIndexValsUnique == 7
        EOTCodeString = 'Force Quit';
    else
        disp('Unknown EOT Code');
    end
else
    EOTCodeString = '';
    for i=1:length(eIndexValsUnique)
        if eIndexValsUnique(i) == 0
            EOTCodeString = cat(2,EOTCodeString,[ 'Correct|']);
        elseif eIndexValsUnique(i) == 1
            EOTCodeString = cat(2,EOTCodeString,[ 'Wrong|']);
        elseif eIndexValsUnique(i) == 2
            EOTCodeString = cat(2,EOTCodeString,[ 'Failed|']);
        elseif eIndexValsUnique(i) == 3
            EOTCodeString = cat(2,EOTCodeString,[ 'Broke|']);
        elseif eIndexValsUnique(i) == 4
            EOTCodeString = cat(2,EOTCodeString,[ 'Ignored|']);
        elseif eIndexValsUnique(i) == 5
            EOTCodeString = cat(2,EOTCodeString,[ 'False Alarm|']);
        elseif eIndexValsUnique(i) == 6
            EOTCodeString = cat(2,EOTCodeString,[ 'Distracted|']);
        elseif eIndexValsUnique(i) == 7
            EOTCodeString = cat(2,EOTCodeString,[ 'Force Quit|']);
        else
            disp('Unknown EOT Code');
        end
    end
    EOTCodeString = [EOTCodeString 'all'];
end
end
function attendLocString = getAttendLocString(aIndexValsUnique)

if length(aIndexValsUnique)==1
    if aIndexValsUnique == 0
        attendLocString = '0 (right)';
    elseif aIndexValsUnique == 1
        attendLocString = '1 (left)';
    else
        disp('Unknown attended location');
    end
else
    attendLocString = '';
    for i=1:length(aIndexValsUnique)
        if aIndexValsUnique(i) == 0
            attendLocString = cat(2,attendLocString,[ '0 (right)|']);
        elseif aIndexValsUnique(i) == 1
            attendLocString = cat(2,attendLocString,[ '1 (left)|']);
        else
            disp('Unknown attended location');
        end
    end
    attendLocString = [attendLocString 'Both'];
end
end
function stimTypeString = getStimTypeString(sIndexValsUnique)

if length(sIndexValsUnique)==1
    if sIndexValsUnique == 0
        stimTypeString = 'Null';
    elseif sIndexValsUnique == 1
        stimTypeString = 'Correct';
    elseif sIndexValsUnique == 2
        stimTypeString = 'Target';
    elseif sIndexValsUnique == 3
        stimTypeString = 'FrontPad';
    elseif sIndexValsUnique == 4
        stimTypeString = 'BackPad';
    else
        disp('Unknown Stimulus Type');
    end
else
    stimTypeString = '';
    for i=1:length(sIndexValsUnique)
        if sIndexValsUnique(i) == 0
            stimTypeString = cat(2,stimTypeString,['Null|']);
        elseif sIndexValsUnique(i) == 1
            stimTypeString = cat(2,stimTypeString,['Correct|']);
        elseif sIndexValsUnique(i) == 2
            stimTypeString = cat(2,stimTypeString,['Target|']);
        elseif sIndexValsUnique(i) == 3
            stimTypeString = cat(2,stimTypeString,['FrontPad|']);
        elseif sIndexValsUnique(i) == 4
            stimTypeString = cat(2,stimTypeString,['BackPad|']);
        else
            disp('Unknown Stimulus Type');
        end
    end
    stimTypeString = [stimTypeString 'all'];
end

end
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data
function [analogChannelsStored,timeVals,goodStimPos] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
end
function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
fileName = fullfile(folderSpikes,'spikeInfo.mat');
if exist(fileName,'file')
    load(fileName);
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end
function [parameterCombinations,cValsUnique,tValsUnique,eValsUnique,...
    aValsUnique,sValsUnique] = loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));
end
function stimResults = loadStimResults(folderExtract)
load (fullfile(folderExtract,'stimResults'));
end