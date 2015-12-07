function [x,xVal] = plotPsychometricFunction(X,Y,colorName,startPt)

uniqueX = unique(X);
lX = length(uniqueX);
xVal = zeros(1,lX);
yVal = zeros(1,lX);
numVal = zeros(1,lX);
for i=1:lX
    xVal(i) = uniqueX(i);
    yVal(i) = mean(Y(X==uniqueX(i)));
    numVal(i) = length(find(X==uniqueX(i)));
end

plot(xVal,yVal,'Marker','o','color',colorName,'lineStyle','none'); hold on;

% Fit data
opts = optimset('TolX',1e-6,'TolFun',1e-6,'MaxIter',500,...
    'Display','off','LargeScale','off','MaxFunEvals',500);

x = fminsearch(@(x) fitWeibullcdf(x,xVal,yVal,numVal),startPt,opts);

xs = 0:(xVal(2)-xVal(1))/100:xVal(end);
ys = x(3)*wblcdf(xs,x(1),x(2));
plot(xs,ys,'color',colorName);

end

function f = fitWeibullcdf(x,xVal,yVal,numVal)

a = x(1); b = x(2); c = x(3);
f = sum(numVal.* ((yVal-c*wblcdf(xVal,a,b)).^2));
end