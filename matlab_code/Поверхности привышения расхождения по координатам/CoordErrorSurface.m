

load('GenV.mat')
bm = 1e-6;
ll = linspace(0.5/12, 3, 12);
%dd = linspace(0, 0.5, 12);
%ll = 0.25/2;
dd = 0;

curve_x = zeros(1, numel(ll));
curve_y = zeros(1, numel(ll));
curve_z = zeros(1, numel(ll));
for k = 1:numel(ll)
    
    for m = 1:numel(dd)
    l = ll(k);
    GenVS = CreateABParameterGenerator({'l', [dd(m) + l/2, 0, 0], 'cos', l, [1, 0, 0];...
                                        'l', [-dd(m) - l/2, 0, 0], 'cos', l, [1, 0, 0];...
                                        'l', [0, dd(m) + l/2, 0], 'sin', l, [0, 1, 0];...
                                        'l', [0, -dd(m) - l/2, 0], 'sin', l, [0, 1, 0]});
    n = 50;
    alf = linspace(1e-3, pi/1.99, n);
    bet = linspace(1e-3, pi/1.99, n);
    erbcord_x = zeros(n, n);
    erbcord_y = zeros(n, n);
    erbcord_z = zeros(n, n);
    target = 0.1;

    sigs = zeros(n, n);

    max_x = zeros(1, n);
    max_y = zeros(1, n);
    max_z = zeros(1, n);
    %alf_m = 0;
    %bet_m = 0;
    parfor i = 1:n
        xx = zeros(1, n);
        yy = zeros(1, n);
        zz = zeros(1, n);
        ssigs = zeros(1, n);
        maxx_x = 0;
        maxx_y = 0;
        maxx_z = 0;
        for j = 1:n

            f = @(t) CalcErr(bm, l, t*cos(alf(i))*cos(bet(j)), t*sin(alf(i))*cos(bet(j)), t*sin(bet(j)), zeros(size(t)), zeros(size(t)), zeros(size(t)), GenV, GenVS, [13, 13]);
            %{1 rtrue, 2 xtrue, 3 ytrue, 4 ztrue, 5 rsol, 6 xsol, 7 ysol, 8 zsol, 9 rer, 10 xer, 11 yer, 12 zer, 13 rreler, 14 xreler, 15 yreler, 16 zreler}
            %tic

            [r, v] = Find(f, target, 1e-2, 15*l+7*dd(m)-0.0156);
            %toc
            %disp("r = "+num2str(r) + " | v = "+num2str(v))

            xx(j) = r.*cos(alf(i))*cos(bet(j));
            yy(j) = r.*sin(alf(i))*cos(bet(j));
            zz(j) = r.*sin(bet(j));
            %{
            if i > 23 && i < 27 && j == 2
                Graphic(@CalcErr, bm, l, xx(j), yy(j), zz(j), GenV, GenVS);
            end
            %}
            if r < 0
                xx(j) = nan;
                yy(j) = nan;
                zz(j) = nan;
            end

            if r > sqrt(maxx_x^2 + maxx_y^2 + maxx_z^2)
                maxx_x = xx(j);
                maxx_y = yy(j);
                maxx_z = zz(j);
                %alf_m = alf(i);
                %bet_m = bet(j);
            end


            %erV(i, j) = sqrt(sum((Gen_s(bm/l, 0, 0, 0, yy(j), xx(j), zz(j)) - GenV(bm, 0, 0, 0, yy(j), xx(j), zz(j))).^2, 'all'));
            %{
            r = linspace(0.01, 10, 100);
            figure;
            semilogy(r, abs(CalcErr(bm, l, r.*cos(alf(i)).*cos(bet(j)), r.*sin(alf(i)).*cos(bet(j)), r.*sin(bet(j)), zeros(size(r)), zeros(size(r)), zeros(size(r)), GenV, GenVS)), 'red' );
            %}
        end
        erbcord_x(i, :) = xx;
        erbcord_y(i, :) = yy;
        erbcord_z(i, :) = zz;
        sigs(i, :) = ssigs;
        disp(i)

        max_x(i) = maxx_x;
        max_y(i) = maxx_y;
        max_z(i) = maxx_z;
    end
    [a, ind] = max(max_x.^2 + max_y.^2 + max_z.^2);
    curve_x(k, m) = max_x(ind);
    curve_y(k, m) = max_y(ind);
    curve_z(k, m) = max_z(ind);
    %disp("alf = " + alf_m + " | bet = " + bet_m)
    %{
    tMap = repmat((sqrt(erbcord_x.^2 + erbcord_y.^2 + erbcord_z.^2) - 0.9*min(sqrt(erbcord_x.^2 + erbcord_y.^2 + erbcord_z.^2), [], 'all')), 1, 1, 3);
    tMap = -0.1+tMap/max(tMap, [], 'all');
    fg1 = figure;
    surface(erbcord_x, erbcord_y, erbcord_z, tMap, "EdgeColor", 'none', "FaceColor", 'texturemap')
    hold on
    %plot3(curve_x(k), curve_y(k), curve_z(k), "x black")
    xlim([0, max(erbcord_x, [], 'all')]); ylim([0, max(erbcord_y, [], 'all')]); zlim([0, max(erbcord_z, [], 'all')]);
    title("Поверхность превышения погрешности в " + num2str(100*target) + " % R, L = " + num2str(round(ll(k),2)))
    xlabel('x')
    hold off
    %}
    disp(k);
    end
end
%%

figure;
plot3(curve_x, curve_y, curve_z, "x");
title("Граничный радиус. Точки максимума.")

figure;
plot(ll, sqrt(curve_x.^2 + curve_y.^2 + curve_z.^2), 'k');
title("Граничный радиус.")

%{
tMap = repmat((sqrt(curve_x.^2 + curve_y.^2 + curve_z.^2) - 0.9*min(sqrt(curve_x.^2 + curve_y.^2 + curve_z.^2), [], 'all')), 1, 1, 3);
tMap = -0.1+tMap/max(tMap, [], 'all');

figure;
surf(ll, dd, sqrt(curve_x.^2 + curve_y.^2 + curve_z.^2), tMap, "EdgeColor", 'none', "FaceColor", 'texturemap')
xlabel('Длинна соленоидов [м].'); ylabel('Расстояние между соленоидами [м].'); title("Граничный радиус.");
%}
%Graphic(@CalcErr, bm, l, 0.2654, 0.0652, 0.7102, GenV, GenVS);
%Graphic(@CalcErr, bm, l, 0.1969, 0.1955, 0.01355, GenV, GenVS)
%Graphic(@CalcErr, bm, l, 0.2886, 0.9553, 0.06438, GenV, GenVS);
%%

function [r] = fun(f, t)
a = f(t);
r = a{:};
end

function [alpha, beta] = Graphic(f, bm, l, x, y, z, GenV, GenVS)

bet = asin(z/sqrt(x^2 + y^2 + z^2));
alf = acos(x/sqrt(x^2 + y^2 + z^2)*cos(bet));

r = linspace(0.01, 0.5, 1e3);
%figure;
var = f(bm, l, r.*cos(alf).*cos(bet), r.*sin(alf).*cos(bet), r.*sin(bet), zeros(size(r)), zeros(size(r)), zeros(size(r)), GenV, GenVS, [13, 13]);
%{1 rtrue, 2 xtrue, 3 ytrue, 4 ztrue, 5 rsol, 6 xsol, 7 ysol, 8 zsol, 9 rer, 10 xer, 11 yer, 12 zer, 13 rreler, 14 xreler, 15 yreler, 16 zreler}

figure;
plot(r, var{:})
hold on
plot([min(r), max(r)], [0.1 0.1])
plot(sqrt(sum([x,y,z].^2)), 0.1, "x")
title("\alpha = " + num2str(alf) + " \beta = " + num2str(bet))
hold off

disp(alf);
disp(bet);
end