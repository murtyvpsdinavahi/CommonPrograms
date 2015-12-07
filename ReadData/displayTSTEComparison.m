% Compares TrialStart (TS) and TrialEnd (TE) obtained from different data
% streams and uses LL data to constrcut stimResults

function maxDiffTrialStartTimePercent=displayTSTEComparison(folderExtract)

load(fullfile(folderExtract,'LL.mat'));
load(fullfile(folderExtract,'digitalEvents.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trialStart
clear xD xL
[~,xD] = getDigitalData(digitalCodeInfo,'TS');
xL = LL.startTime;

maxDiffTrialStartTimePercent=compareTimes(xD,xL,'Start');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TrialEnd: Match EOTCodes
if (sum(getDigitalData(digitalCodeInfo,'TE'))>0) % When the digital codes contain the values
    xD = getDigitalData(digitalCodeInfo,'TE');
    xL=LL.eotCode;
    compareCodes(xD,xL);
else
    disp('EOT Codes not available in the digital data, using trialEnd times instead...');
    [~,xD] = getDigitalData(digitalCodeInfo,'TE');
    xL = LL.endTime;
    compareTimes(xD,xL,'End');
end
end
function [value,time] = getDigitalData(digitalCodeInfo,codeName)

for i=1:length(digitalCodeInfo)
    if strcmp(digitalCodeInfo(i).codeName,codeName)
        endPos=i;
        break;
    end
end

if ~exist('endPos','var')
    error(['No trialEvent named ' codeName]);
else
    value=digitalCodeInfo(endPos).value;
    time=digitalCodeInfo(endPos).time;
end

end
function maxDiffTimePercent = compareTimes(xD,xL,str)

xD = diff(xD); lxD = length(xD); xD=xD(:);
xL = diff(xL); lxL = length(xL); xL=xL(:);

if strcmpi(str,'Start')
    h1 = subplot(221); h2 = subplot(223);
else
    h1 = subplot(222); h2 = subplot(224);
end

if lxD == lxL
    disp(['Number of ' str ' Trials: ' num2str(lxD+1)]);
    axes(h1);
    plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
    ylabel(['Difference in ' str ' Times (s)']);
    legend('Dig','LL','Location','SouthEast');
    
    axes(h2);
    plot(1000*(xD-xL),'b');
    ylabel('Digital-Lablib times (ms)');
    xlabel('Trial Number');
    
    maxDiffTimePercent = 100*max(abs(xD-xL) ./ xD);
else
    error(['Num of ' str 'Trials: digital: ' num2str(lxD+1) ' , LL: ' num2str(lxL+1)]);
end
end
function compareCodes(xD,xL)

lxD = length(xD); xD=xD(:);
lxL = length(xL); xL=xL(:);

if lxD == lxL
    disp(['Number of eotCodes: ' num2str(lxD)]);
    subplot(222)
    plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
    ylabel('eotCode number');
    %legend('Dig','LL','Location','SouthEast');
    
    subplot(224)
    plot(xD-xL,'b');
    ylabel('\delta eotCode number');
    xlabel('Trial Number');
    
    if max(abs(xD-xL))>0
        error('EOT codes do not match');
    end
else
    error(['Number of stimOnset: digital: ' num2str(lxD) ' , LL: ' num2str(lxL)]);
end
end