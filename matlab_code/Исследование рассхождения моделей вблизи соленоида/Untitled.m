m = 1;
N = 100;
l = 0.25;
x = linspace(-4, 4, N);
y = linspace(-4, 4, N);
z = zeros(size(x));
dip = MagneticMoment(0, 0, 0, m, 0, 0);
Line = CreateABParameterGenerator({'l', [0.5, 0, 0], 'cos', l, [1, 0, 0];...
                                    'l', [-0.5, 0, 0], 'cos', l, [1, 0, 0];...
                                    'l', [0, 0.5, 0], 'sin', l, [0, 1, 0];...
                                    'l', [0, -0.5, 0], 'sin', l, [0, 1, 0]});
alf = 0;
bet = 0;
gam = 0;
bx = zeros(N,N);
by = zeros(N,N);
bz = zeros(N,N);



bx3 = zeros(N, N);
by3 = zeros(N, N);
bz3 = zeros(N, N);

for i = 1:numel(x)
    for j = 1:numel(y)
        [bx(i, j), by(i, j), bz(i, j)] = dip.Field(x(i), y(j), 0);

        [bx3(i, j), by3(i, j), bz3(i, j)] = CalcAmpls(x(i), y(j), z(1), alf, bet, gam, m, l, Line);
    end
end

%%
figure;
[X, Y] = meshgrid(x, y);
surf(X, Y, sqrt(bx3.^2 + by3.^2 + bz3.^2));

%%
figure;
nrm = sqrt(bx.^2 + by.^2);
[X, Y] = meshgrid(x, y);
%quiver(X(1:n:end), Y(1:n:end), bx(1:n:end)./nrm(1:n:end), -by(1:n:end)./nrm(1:n:end))

%strmsl1 = streamslice(X, Y, -bx, by);
hold on


nrm = sqrt(bx3.^2 + by3.^2);
strmsl3 = streamslice(X, Y, bx3, by3);
for i = 1:numel(strmsl3)
    strmsl3(i).Color = [1,0,0];
end


%%
%{
figure;
beacon = MagneticMoment(0, 0, 0, 1, 0, 0);
f = @(x, y, z) abs(beacon.Field(x,y,z) - fn_Bx(1, x, y, z))-1e-10;
fimplicit3(f, [-10 10 -10 10 -10 10])
title('Bx 1e-10 L = 1')
hold on
xlabel('x')
%}
alf = 0;
bet = 0;
gam = 0;
bm = 1;
L = 0.25;
Gen_s = CreateABParameterGenerator({'l', 'cos', L, [1, 0, 0]; 'l', 'sin', L, [0, 1, 0]});

figure;
f = @(x, y, z) CalcRelativeError(x, y, z, alf, bet, gam, bm, L, Gen_s)-50;
fimplicit3(f, [0 10 0 10 0 10])
title('Поверхность привышения расхождения по x в 10 процентов от r')
hold on
xlabel('x')


%%

debug_fg = figure('Name', 'debug');
hold on
x0 = 1;
y0 = 1;
z0 = 1;
alf = deg2rad(1);
bet = deg2rad(0);
gam = deg2rad(0);
R =...
   [cos(gam ),   -sin(gam ),  0;...
    sin(gam ),   cos(gam ),   0;...
    0,          0,          1]*...
   [cos(bet ),   0,          sin(bet );...
    0,          1,          0;...
    -sin(bet ),  0,          cos(bet );]*...
   [1,          0,          0;...
    0,          cos(alf ),   -sin(alf );...
    0,          sin(alf ),   cos(alf );];
px = 1.45;
py = 0.43;
pz = -2;
plot3(px, py, pz, "xblue")
p = [px, py, pz];
plot3(0,0,0,"xblack")
plot3([0 1].*px, [0 0].*py, [0 0].*pz, "black")
plot3([0 0].*px, [0 1].*py, [0 0].*pz, "black")
plot3([0 0].*px, [0 0].*py, [0 1].*pz, "black")
plot3([1 1].*px, [0 0].*py, [0 1].*pz, "black")
plot3([1 1].*px, [0 1].*py, [0 0].*pz, "black")
plot3([0 1].*px, [1 1].*py, [0 0].*pz, "black")
plot3([0 0].*px, [1 1].*py, [0 1].*pz, "black")
plot3([0 1].*px, [0 0].*py, [1 1].*pz, "black")
plot3([0 0].*px, [0 1].*py, [1 1].*pz, "black")
plot3([1 1].*px, [1 0].*py, [1 1].*pz, "black")
plot3([1 1].*px, [1 1].*py, [1 0].*pz, "black")
plot3([1 0].*px, [1 1].*py, [1 1].*pz, "black")
 
alt_p = R*(p-[x0,y0,z0])';
plot3(x0,y0,z0,"xred")
xx = R(1,:).*alt_p(1);
yy = R(2,:).*alt_p(2);
zz = R(3,:).*alt_p(3);
plot3([x0, x0+xx(1)],[y0, y0+xx(2)],[z0, z0+xx(3)], "red")
plot3([x0, x0+yy(1)],[y0, y0+yy(2)],[z0, z0+yy(3)], "red")
plot3([x0, x0+zz(1)],[y0, y0+zz(2)],[z0, z0+zz(3)], "red")

%% Функции
function err = CalcRelativeError(x, y, z, alf, bet, gam, bm, l, Gen_s)
arguments
    x (1, :) {mustBeNumeric, mustBeReal}
    y (1, :) {mustBeNumeric, mustBeReal}
    z (1, :) {mustBeNumeric, mustBeReal}
    alf (1, 1) {mustBeNumeric, mustBeReal}
    bet (1, 1) {mustBeNumeric, mustBeReal}
    gam (1, 1) {mustBeNumeric, mustBeReal}
    bm (1, 1) {mustBeNumeric, mustBeReal}
    l (1, 1) {mustBeNumeric, mustBeReal}
    Gen_s (1, 1) function_handle
end
err = zeros(size(x));

for i = 1:numel(x)
    sol_v = Gen_s(bm/l, alf, bet, gam, x(i), y(i), z(i));
    err(i) = 100*abs((CalcX(bm, sol_v(1,1), sol_v(2,1), sol_v(3,1), sol_v(1,2), sol_v(2,2), sol_v(3,2)) - x)/sqrt(x^2 + y^2 + z^2));
end
end

function [Ax, Ay, Az] = CalcAmpls(x, y, z, alf, bet, gam, bm, l, Gen_s)
arguments
    x (1, :) {mustBeNumeric, mustBeReal}
    y (1, :) {mustBeNumeric, mustBeReal}
    z (1, :) {mustBeNumeric, mustBeReal}
    alf (1, 1) {mustBeNumeric, mustBeReal}
    bet (1, 1) {mustBeNumeric, mustBeReal}
    gam (1, 1) {mustBeNumeric, mustBeReal}
    bm (1, 1) {mustBeNumeric, mustBeReal}
    l (1, 1) {mustBeNumeric, mustBeReal}
    Gen_s (1, 1) function_handle
end

Ax = zeros(size(x));
Ay = zeros(size(x));
Az = zeros(size(x));

for i = 1:numel(x)
    sol_v = Gen_s(bm/l, alf, bet, gam, y(i), x(i), z(i));
    %err(i) = 100*abs((CalcX(bm, sol_v(1,1), sol_v(2,1), sol_v(3,1), sol_v(1,2), sol_v(2,2), sol_v(3,2)) - x)/sqrt(x^2 + y^2 + z^2));
    Ax(i) = sqrt(sol_v(1,1)^2 + sol_v(1,2)^2);
    Ay(i) = sqrt(sol_v(2,1)^2 + sol_v(2,2)^2);
    Ay(i) = sqrt(sol_v(3,1)^2 + sol_v(3,2)^2);
end
end