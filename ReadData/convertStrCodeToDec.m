function decList = convertStrCodeToDec(codeInStrList,useSingelITC18Flag,useSimpleCodeFlag)

if ~exist('useSingelITC18Flag','var');      useSingelITC18Flag=1;       end
if ~exist('useSimpleCodeFlag','var');       useSimpleCodeFlag=0;        end

numCodes = size(codeInStrList,1);
decList = zeros(1,numCodes);

for i=1:numCodes
    if useSimpleCodeFlag
        if strcmpi(codeInStrList(i,:),'TS')
            decList(i)= 2;
        elseif strcmpi(codeInStrList(i,:),'FI')
            decList(i)= 4;
        elseif strcmpi(codeInStrList(i,:),'ON')
            decList(i)= 8;
        elseif strcmpi(codeInStrList(i,:),'OF')
            decList(i)= 16;
        elseif strcmpi(codeInStrList(i,:),'TO')
            decList(i)= 32;
        elseif strcmpi(codeInStrList(i,:),'SA')
            decList(i)= 64;
        elseif strcmpi(codeInStrList(i,:),'TE')
            decList(i)= 128;
        end 
    else
        if useSingelITC18Flag       %8-14 represent the first letter, 1-7 represent the second.
            decList(i)=256*double(codeInStrList(i,1))+2*double(codeInStrList(i,2)); %#ok<*AGROW>
        else                        %8 bits are used for a hex letter
            decList(i)=256*double(codeInStrList(i,1))+double(codeInStrList(i,2));
        end
    end
end
end