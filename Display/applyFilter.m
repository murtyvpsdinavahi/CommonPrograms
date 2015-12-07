function [y,filterStr] = applyFilter(x,Fs,filterKind,filterType,filterOrder,cutoff)

if ~exist('Fs','var');                  Fs = 2000;                      end
if ~exist('filterKind','var');          filterKind = 'butter';          end
if ~exist('filterType','var');          filterType = 'high';            end
if ~exist('filterOrder','var');         filterOrder = 1;                end
if ~exist('cutoff','var');              cutoff = 3;                     end

if strcmp(filterKind,'butter')
    [B,A] = butter(filterOrder,cutoff/(Fs/2),filterType);  % 1st order high pass
    filterStr = ['_' filterKind '_' filterType '_order' num2str(filterOrder) '_cutoff' num2str(cutoff)]; 
end

if ~isempty(x)
    y = filtfilt(B,A,x')';
end
end