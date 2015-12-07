function displayElectrodeImpedance(subjectName,expDates,gridType,folderSourceString,electrodeList)

if ~exist('electrodeList','var')
    if strcmp(gridType,'ECoG')
        electrodeList = [1:34 65:96];
    else
        electrodeList = 1:96;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
impedanceGrid = [0.05 0.05 0.6 0.9];
impedanceVsDaysGrid = [0.7 0.8 0.25 0.15];
tableGrid = [0.7 0.05 0.25 0.5]; numTableColumns = 3;
electrodeGrid = [0.7 0.56 0.25 0.175]; 

numPlots = length(expDates);
dY = impedanceGrid(4)/numPlots;

impedanceHighCutOff = 2500;
impedanceModerateCutOff = 1500;
impedanceYLim=2500;

goodExpDates='';

meanImpedances =zeros(1,numPlots);
semImpedances  =zeros(1,numPlots);

for i=1:numPlots
    
    subplot('Position',[impedanceGrid(1) impedanceGrid(2)+(numPlots-i)*dY impedanceGrid(3) dY]); 
    expDate = expDates{i};
    
    clear impedanceValuesAll impedanceValues
    fileName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,'impedanceValues.mat');
    impedanceValuesAll = getImpedanceValues(fileName);
    impedanceValues = impedanceValuesAll(electrodeList);
    
    stem(electrodeList,impedanceValues); ylim([0 impedanceYLim]);
    hold on;
    plot(electrodeList,impedanceModerateCutOff,'k--');
    
    if i~=numPlots
        set(gca,'XTickLabel',[]);
         goodExpDates=cat(2,goodExpDates,[expDate '|']);
    else
        xlabel('Electrode Number');
        goodExpDates=cat(2,goodExpDates,expDate);
    end
    ylabel('K\Omega'); 
    
    % Write stuff
    clear moderatelyBadChannelPos badChannelPos
    moderatelyBadChannelPos = find(impedanceValues>impedanceModerateCutOff);
    stem(electrodeList(moderatelyBadChannelPos),impedanceValues(moderatelyBadChannelPos),'m');
    badChannelPos = find(impedanceValues>impedanceHighCutOff);
    stem(electrodeList(badChannelPos),impedanceValues(badChannelPos),'r');
    text(0.025,0.9,expDate,'Unit','Normalized');
    
    goodImpedances = impedanceValues(impedanceValues<impedanceHighCutOff);
    meanImpedances(i) = mean(goodImpedances);
    semImpedances(i) = std(goodImpedances)/sqrt(length(goodImpedances));
end

hImpedance=subplot('Position',impedanceVsDaysGrid);
if ~exist('dayIndices','var')
    dayIndices=1:numPlots;
end
errorbar(hImpedance,dayIndices,meanImpedances,semImpedances);

% Make expDate chooser
hExpDate = uicontrol('Unit','Normalized', 'Position',[0.7 0.9 0.25 0.1], ...
    'Style','popup','String',goodExpDates,'FontSize',10,'Callback',{@plot_Callback});

    function plot_Callback(~,~)
        num = get(hExpDate,'val');
        clear expDate impedanceValues
        expDate = expDates{num};
        
        clear impedanceValuesAll impedanceValues
        fileName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,'impedanceValues.mat');
        impedanceValuesAll = getImpedanceValues(fileName);
        impedanceValues = impedanceValuesAll(electrodeList);
        
        moderatelyBadChannels = electrodeList(impedanceValues>impedanceModerateCutOff);
        badChannels = electrodeList(impedanceValues>impedanceHighCutOff);
        
        showElectrodeLocations(electrodeGrid,moderatelyBadChannels,'m',[],0,0,gridType,subjectName);
        showElectrodeLocations(electrodeGrid,badChannels,'r',[],1,0,gridType,subjectName);
        
        numRows = ceil(length(electrodeList)/numTableColumns);
        entries = cell(numRows,numTableColumns);
        for ii=1:numRows
            for jj=1:numTableColumns
                entries{ii,2*jj-1} = num2str(electrodeList(ii+numRows*(jj-1)));
                entries{ii,2*jj} = num2str(impedanceValues(ii+numRows*(jj-1)));
            end
        end
        makeTable(tableGrid,entries,1); hold off
    end
end
function impedanceValues = getImpedanceValues(fileName) %#ok<STOUT>
load(fileName);
end