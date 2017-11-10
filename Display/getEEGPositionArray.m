
function electrodeArray = getEEGPositionArray(capName,refType)

% Get chanlocs
if strcmpi(refType,'unipolar')
    chL = load(fullfile(pwd,'Montages','Layouts',capName,[capName '.mat']));
    chanlocs = chL.chanlocs;
elseif strcmpi(refType,'bipolar')
    chL = load(fullfile(pwd,'Montages','Layouts',capName,['bipolarChanlocs' upper(capName(1)) capName(2:end) '.mat']));
    chanlocs = chL.eloc;
end

% Get x and y co-ordinates
[~,~,Th,Rd] = readlocs(chanlocs);
Th = pi/180*Th; 
[x,y]     = pol2cart(Th,Rd);

% Shift origin to 0,0 of figure
x = x + 0.5;
y = y + 0.5;

% Find row and column numbers for each electrode by rescaling x and y and
% rounding them off to mearest integer.
colNum = round(y*100);
rowNum = repmat(100,1) - round(x*100);

% Create electrodeArray
electrodeArray = zeros(100,100);

for iElec = 1:length(x)
    electrodeArray(rowNum(iElec),colNum(iElec)) = iElec;
end

% Remove extra zeros
electrodeArray(:,all(electrodeArray==0,1)) = [];
electrodeArray(all(electrodeArray==0,2),:) = [];

end