% This function is used to extract the digital data from raw .nev data files


function [hFile,digitalTimeStamps,digitalEvents]=extractDigitalDataBlackrock(subjectName,expDate,protocolName,folderSourceString,gridType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = [subjectName expDate protocolName '.nev'];

dataLog{1,2} = subjectName;
dataLog{2,2} = gridType;
dataLog{3,2} = expDate;
dataLog{4,2} = protocolName;

[~,folderName]=getFolderDetails(dataLog);

% folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
makeDirectory(folderName);
folderIn = fullfile(folderSourceString,'data','rawData',[subjectName expDate]);
folderExtract = fullfile(folderName,'extractedData');
makeDirectory(folderExtract);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the NEV file

% Load the appropriate DLL
dllName = 'C:\Users\LabComputer6\Documents\MATLAB\Programs\SRAYLab Programs\SoftwareMAP (Only required ones)\NeuroShare\nsNEVLibrary64.dll';

[nsresult] = ns_SetLibrary(dllName);
if (nsresult ~= 0);      error('DLL was not found!');                   end

% Load data file and display some info about the file open data file
[nsresult, hFile] = ns_OpenFile(fullfile(folderIn,fileName));
if (nsresult ~= 0);      error('Data file did not open!');              end

% Get file information
[nsresult, fileInfo] = ns_GetFileInfo(hFile);
% Gives you entityCount, timeStampResolution and timeSpan
if (nsresult ~= 0);      error('Data file information did not load!');  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Digital Codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, entityInfo] = ns_GetEntityInfo(hFile,1:fileInfo.EntityCount);
eventList = find([entityInfo.EntityType] == 1);

digitalEventName = 'digin';
for i =1:length(eventList)
    if strcmp(digitalEventName,[entityInfo(eventList(i)).EntityLabel]);
        digID=i;
        break;
    end
end

if ~exist('digID','var')    
    error(['No event named ' digitalEventName]);
else
    numDigitalEvents=entityInfo(eventList(digID)).ItemCount;
end
disp(['Number of digital events: ' num2str(numDigitalEvents)]);

% Get the digital events
[~,digitalTimeStamps,digitalEvents] = ns_GetEventData(hFile,eventList(digID),1:numDigitalEvents);

% Blackrock collects data at 30,000 samples per second. Sometimes a digital
% code transition takes longer than this, and is counted twice. First we
% find out the double counts.

deltaLimit = 1.5/30000; 
dt = diff(digitalTimeStamps);
badDTPos = find(dt<=deltaLimit);

if ~isempty(badDTPos)
    disp([num2str(length(badDTPos)) ' of ' num2str(length(digitalTimeStamps)) ' (' num2str(100*length(badDTPos)/length(digitalTimeStamps),2) '%) are repeats and will be discarded']);
    digitalTimeStamps(badDTPos)=[];
    digitalEvents(badDTPos)=[];
end

save(fullfile(folderExtract,'\NEVFileInfo.mat'), 'fileInfo', 'entityInfo');
end