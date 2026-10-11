clear; close all; clc
load('..\GenV.mat')


nt = 25;
phase = linspace(0, 2*pi,  nt);

n = 2.5e3;
bias = 0;
x = linspace(-4, 4, n);
y = linspace(-4, 4, n);

d = 0.1;
l = 0.25;
bm = 1;

GenVS = CreateABParameterGenerator2({'l', [d + l/2, 0, 0], 'cos', l, [1, 0, 0];...
                                    'l', [-d - l/2, 0, 0], 'cos', l, [1, 0, 0];...
                                    'l', [0, d + l/2, 0], 'sin', l, [0, 1, 0];...
                                    'l', [0, -d - l/2, 0], 'sin', l, [0, 1, 0]});
%%
%GenV(bm*l, 0, 0, 0, 5, 5, 5)
%GenVS(bm/2, 0, 0, 0, 5, 5, 5)
var = MgFieldL(0, x, y, bm/2, l, GenVS);
var2 = MgFieldL(0, x, y, bm*l, 1, GenV);
%%
fg = figure;

ax = axes(fg);
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
ax.NextPlot = 'add';
ax.Box = 'on';

ang = pi/12;
of_ang = pi/64;
n_ang = 10;
stx = 1.09*(d+l)*cos([linspace(-ang-of_ang, -of_ang, n_ang), linspace(of_ang, ang+of_ang, n_ang)] + pi/2);
sty = 1.09*(d+l)*sin([linspace(-ang-of_ang, -of_ang, n_ang), linspace(of_ang, ang+of_ang, n_ang)] + pi/2);

stx2 = zeros(1, length(stx));
sty2 = zeros(1, length(stx));

st = stream2(x, y, var(:,:,1), var(:,:,2), stx, sty, [0.1, 1e5]);


plot(ax, [d, d+l], [0, 0], 'black')
plot(ax, [0, 0], [d, d+l], 'black')
plot(ax, [-d, -d-l], [0, 0], 'black')
plot(ax, [0, 0], [-d, -d-l], 'black')


for i = 1:numel(st)
    [m, ind] = min(sqrt( (st{i}(:, 1) - stx(i)).^2 + (st{i}(:, 2) + sty(i)).^2 ));
    [m2, ind2] = min(abs(st{i}(1:ind, 2)));
    stx2(i) = st{i}(ind2, 1);
    
    plot(ax, st{i}(1:ind, 1), st{i}(1:ind, 2), 'b')
end



%fg1 = figure;

%ax1 = axes(fg1);
%ax1.XTick = [];
%ax1.YTick = [];
%ax1.ZTick = [];
%ax1.NextPlot = 'add';
%ax1.Box = 'on';

%var = MgFieldL(0, x, y, bm*l, 1, GenV);

%ang = pi/12;
%of_ang = pi/6;
%n_ang = 10;
%stx = (d+l)*cos([linspace(-ang-of_ang, -of_ang, n_ang), linspace(of_ang, ang+of_ang, n_ang)] + pi/2);
%sty = (d+l)*sin([linspace(-ang-of_ang, -of_ang, n_ang), linspace(of_ang, ang+of_ang, n_ang)] + pi/2);
%stx = 0.1.*stx;
%sty = 0.1.*sty;


st = stream2(x, y, var2(:,:,1), var2(:,:,2), stx2, sty2,[0.1, 1e5]);

for i = 1:numel(st)
    %[m, ind] = min(sqrt( (st{i}(:, 1) - stx(i)).^2 + (st{i}(:, 2) + sty(i)).^2 ));
    [m, ind] = min(sqrt( (st{i}(:, 1)).^2 + (st{i}(:, 2)).^2 ));
    plot(ax, st{i}(1:ind, 1), st{i}(1:ind, 2), '--black')
    plot(ax, st{i}(1:ind, 1), -st{i}(1:ind, 2), '--black')
end
ax.XLim = [-3.6 3.6];
ax.YLim = [-3.6 3.6];
%axis([-3.6, 3.6, -3.6, 3.6])
%% Функции

function var = MgFieldL(t, x, y, bm, l, genVs)

result = zeros(numel(x), numel(y), 3);
for i = 1:numel(x) 
    for j = 1:numel(y)
        ab = genVs(bm/l, 0, 0, 0, y(j), x(i), 1e-3);
        vec = ab*[cos(t); sin(t)];
        result(i, j, 1) = vec(1);
        result(i, j, 2) = vec(2);
        result(i, j, 3) = vec(3);
    end
end

var = result;

end