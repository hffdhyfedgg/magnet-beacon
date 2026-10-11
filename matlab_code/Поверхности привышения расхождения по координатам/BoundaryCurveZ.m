load('GenV.mat')
bm = 1e-6;

n = 75;



target = 0.1;

l = (linspace((0.1), (3), 30));
yy = zeros(numel(l), n);
zz = zeros(numel(l), n);
parfor i = 1:numel(l)
z = zeros(1, n);
y = exp(linspace(log(l(i)/2), log(30), n));
yy(i,:) = y;
GenVS = CreateABParameterGenerator({'l', 'cos', l(i), [1, 0, 0]; 'l', 'sin', l(i), [0, 1, 0]});
    for j = 1:n

        f = @(t) CalcErr(bm, l(i), repmat(1e-3, size(t)), repmat(y(j), size(t)), t, zeros(size(t)), zeros(size(t)), zeros(size(t)), GenV, GenVS, [16, 16]);%{1 rtrue, 2 xtrue, 3 ytrue, 4 ztrue, 5 rsol, 6 xsol, 7 ysol, 8 zsol, 9 rer, 10 xer, 11 yer, 12 zer, 13 rreler, 14 xreler, 15 yreler, 16 zreler}

        [r, v] = Find(f, target, 2e-2, 4*l(i));

        z(j) = r;

    end
zz(i, :) = z;
disp(i)
end
%%
fg = figure;
ax = axes(fg, 'NextPlot', 'add');
hold on
ymax = 10;

ry = ymax/45;
rx = 10/45;
circy = ry*sin(linspace(-pi, pi, 40));
for i = 1:30   
    
    circx = yy(i, end-3)+rx*cos(linspace(-pi, pi, 40));
    
    plot(ax, yy(i,:), zz(i, :), 'k')
    %fill(ax, circx, circy+zz(i, end-3), "white");
    %plot(ax, circx, circy+zz(i, round(n/1.85 ,0)), 'k')
    %text(ax, yy(i, end-3), zz(i, end-3),  num2str(round(l(i), 2)), "HorizontalAlignment", "center" , 'FontSize',7, "VerticalAlignment", "middle");
end
ylim([0, ymax])
grid on

figure;
plot(l, zz(:, end), 'k')