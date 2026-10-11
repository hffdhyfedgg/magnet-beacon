
function [x, y] = Find(f, targetValue, targetAccuracy, start, n,  step, maxCountSteps)
arguments
    f(1,1) function_handle     
    targetValue (1,1) double
    targetAccuracy (1,1) double
    start (1,1) double = 5
    n (1,1) {mustBeInteger} = 10    
    step (1,1) double = 1e-2
    maxCountSteps (1,1) {mustBeInteger} = 1e3
end


startMin = start;
startMax = start;
values = f(startMin);
%{
figure;
yy = f(linspace(1e-3, 1, 100));
plot(linspace(1e-3, 1, 100), abs(yy{:}), 'black')
%}

i =0;
while abs(values{:}) < targetValue && i < (start/step) && startMin > 0
    startMax = startMin;
    startMin = startMin - step;
    values = f(startMin);
    i = i + 1;    
end

xdata = linspace(startMin, startMax, n);
values = f(xdata);
[accuracy, ind] = min(abs(abs(values{:})-targetValue));
i = 0;
while round(accuracy,5) > targetAccuracy && i < maxCountSteps
    
    if abs(values{:}(ind)) > targetValue && ind < n
        xdata = linspace(xdata(ind), xdata(ind+1), n);
        values = f(xdata);
        [accuracy, ind] = min(abs(abs(values{:}) - targetValue));
    elseif ind > n
        xdata = linspace(xdata(ind-1), xdata(ind), n);
        values = f(xdata);
        [accuracy, ind] = min(abs(abs(values{:}) - targetValue));
    else
        x = xdata(ind);
        y = values{:}(ind);
        return
    end
    i = i+1;
end
x = xdata(ind);
y = values{:}(ind);
end