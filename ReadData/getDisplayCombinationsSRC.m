% This function returns a N dimentional cell array containing the stim
% numbers of a particular type of stimulus condition for the SRC protocol. 
% No changes from the previous version used for microelectrode data

function parameterCombinations = getDisplayCombinationsSRC(folderOut,goodStimNums)

load(fullfile(folderOut,'stimResults.mat'));

% Five parameters are chosen:
% 1. Contrast
% 2. Temporal Frequency
% 3. eotCodes
% 4. Attention Location - 0 or 1
% 5. Stimulus type - front padding or valid

% Parameters index
parameters{1} = 'contrastIndex';
parameters{2} = 'temporalFreqIndex';
parameters{3} = 'eotCodes';
parameters{4} = 'attendLoc';
parameters{5} = 'type0'; %#ok<NASGU> % only 1 and 3, these values are same for type1 for goodStimNums 

% get Contrast
cValsAll  = stimResults.contrastIndex;
tValsAll  = stimResults.temporalFreqIndex;
eValsAll  = stimResults.eotCodes;
aValsAll  = stimResults.attendLoc;
sValsAll  = stimResults.type0;

cValsGood = cValsAll(goodStimNums);
tValsGood = tValsAll(goodStimNums);
eValsGood = eValsAll(goodStimNums);
aValsGood = aValsAll(goodStimNums);
sValsGood = sValsAll(goodStimNums);

cValsUnique = unique(cValsGood); cLen = length(cValsUnique);
tValsUnique = unique(tValsGood); tLen = length(tValsUnique);
eValsUnique = unique(eValsGood); eLen = length(eValsUnique);
aValsUnique = unique(aValsGood); aLen = length(aValsUnique);
sValsUnique = unique(sValsGood); sLen = length(sValsUnique);

% If more than one value, make another entry with all values
if (cLen > 1);           cLen=cLen+1;                    end
if (tLen > 1);           tLen=tLen+1;                    end
if (eLen > 1);           eLen=eLen+1;                    end
if (aLen > 1);           aLen=aLen+1;                    end
if (sLen > 1);           sLen=sLen+1;                    end

allPos = 1:length(goodStimNums);

% display
disp(['Number of unique contrasts: ' num2str(cLen)]);
disp(['Number of unique TFs: ' num2str(tLen)]);
disp(['Number of unique eotCodes: ' num2str(eLen)]);
disp(['Number of unique attendLoc: ' num2str(aLen)]);
disp(['Number of unique stimtypes: ' num2str(sLen)]);
disp(['total combinations: ' num2str((cLen)*(tLen)*(eLen)*(aLen)*(sLen))]);

for c=1:cLen
    if c==cLen
        cPos = allPos;
    else
        cPos = find(cValsGood == cValsUnique(c));
    end
    
    for t=1:tLen
        if t==tLen
            tPos = allPos;
        else
            tPos = find(tValsGood == tValsUnique(t));
        end
        
        for e=1:eLen
            if e==eLen
                ePos = allPos;
            else
                ePos = find(eValsGood == eValsUnique(e));
            end
            
            for a=1:aLen
                if a==aLen
                    aPos = allPos;
                else
                    aPos = find(aValsGood == aValsUnique(a));
                end
                
                for s=1:sLen
                    if s==sLen
                        sPos = allPos;
                    else
                        sPos = find(sValsGood == sValsUnique(s));
                    end
                    
                    ctPos = intersect(cPos,tPos);
                    ctePos = intersect(ctPos,ePos);
                    cteaPos = intersect(ctePos,aPos);
                    parameterCombinations{c,t,e,a,s} = intersect(cteaPos,sPos); %#ok<*AGROW>
                end
            end
        end
    end
end

% save 
save(fullfile(folderOut,'parameterCombinations.mat'),'parameters','parameterCombinations', ... 
    'cValsUnique','tValsUnique','eValsUnique','aValsUnique','sValsUnique');