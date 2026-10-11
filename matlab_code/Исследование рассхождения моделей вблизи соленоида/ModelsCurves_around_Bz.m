bm = 1e-7;
l = 0.25;
mynorm = @(r) (sqrt(sum(r.^2)));
step = 0.4;
erer = 0.1;

angles_n = 35;
alf = linspace(1e-5, pi/2 - 1e-5, angles_n);
bet = linspace(1e-5, pi/2 - 1e-5, angles_n);
erbcord_x = zeros(angles_n, angles_n);
erbcord_y = zeros(angles_n, angles_n);
erbcord_z = zeros(angles_n, angles_n);
%{
syms m0 l real positive
syms x y z real positive
r = [x,y,z]; % Радиус-вектор до точки
sym_ro = sym('ro', 'real'); %Плотность магнитного момента
rm = [l,0,0]; %Радиус-вектор до элменета катушки
assumptions
vec_ro = [sym_ro,0,0];
syms bx by bz b real
dB = (m0/(4*pi)) * (3*(r-rm)*dot(vec_ro,(r-rm))/mynorm(r-rm)^5-vec_ro/mynorm(r-rm)^3);
%Часть интегрирование Bx
fn = 1/((l-x)^2+(y)^2+(z)^2)^sym(3/2);
half_Bx1(l) = int(fn,l);
fn = (l-x)*(l-x)/(mynorm([l-x,y,z]))^5;
half_Bx2(l) = 3*int(fn,l);



disp("Магнитный момент равен: " + num2str(m)+ "[А*м^2] | Длинна соленоида равна: " + num2str(var_L));
angles_n = 100;
alf = deg2rad( linspace(0.01, 89.99, angles_n));
bet = deg2rad( linspace(0.01, 89.99, angles_n));
R0 = 0.2;

erer = 1e-1;
L = var_L;
L0 = -L/2;
ro = m/L;


%disp(num2str(i) + " - ro = " + num2str(ro) + "[А*м]" + " | L = " + num2str(L) + "[м]")

%Заканчиваем интегрировать Bx
Bx = (m0*sym_ro/(4*pi))*((half_Bx2(L0+L)-half_Bx2(L0))-(half_Bx1(L0+L)-half_Bx1(L0)));
fn_Bx = matlabFunction(subs(simplify(Bx), m0, 4e-7*pi));

%Интегрирование By
By = simplify(subs(int(dB(2), l), l, L0+L) - subs(int(dB(2), l), l, L0));
fn_By = matlabFunction(subs(simplify(By), m0, 4e-7*pi));

%Интегрирование Bz
Bz = simplify(subs(int(dB(3),l),l,L0+L) - subs(int(dB(3),l),l,L0));
fn_Bz = matlabFunction(subs(simplify(Bz), m0, 4e-7*pi));
beacon = MagneticMoment(0, 0, 0, ro*L, 0, 0);

%}
Gen_s = CreateABParameterGenerator({'l', 'cos', l, [1, 0, 0]; 'l', 'sin', l, [0, 1, 0]});
%erV = zeros(angles_n, angles_n);
for i = 1:angles_n
    xx = zeros(1, angles_n);
    yy = zeros(1, angles_n);
    zz = zeros(1, angles_n);
    for j = 1:angles_n

        f = @(t) CalcAmplErr(t*cos(alf(i))*cos(bet(j)), t*cos(alf(i))*sin(bet(j)), t*sin(bet(j)), 0, 0, 0, bm, l, Gen_s, GenV);
        [r, v] = FindValue(f, 10, erer, step, 4, 100);
        xx(j) = r*cos(alf(i))*cos(bet(j));
        yy(j) = r*sin(alf(i))*cos(bet(j));
        zz(j) = r*sin(bet(j));
        %erV(i, j) = sqrt(sum((Gen_s(bm/l, 0, 0, 0, yy(j), xx(j), zz(j)) - GenV(bm, 0, 0, 0, yy(j), xx(j), zz(j))).^2, 'all'));

    end
    erbcord_x(i, :) = xx;
    erbcord_y(i, :) = yy;
    erbcord_z(i, :) = zz;
    disp(i);
end

%%
fg1 = figure;
surface(erbcord_x, erbcord_y, erbcord_z, "EdgeColor", 'none')
xlim([0, max(erbcord_x, [], 'all')]); ylim([0, max(erbcord_y, [], 'all')]); zlim([0, max(erbcord_z, [], 'all')]);
title("Поверхность превышения погрешности в 10 % по норме амплитуд")
xlabel('x')
%{
fg2 = figure;
for j = 1:angles_n
    plot3(erbycord(1,:, j), erbycord(2,:, j), erbycord(3,:, j))
    if j == 1
        hold on
    end
end
%}

%% fimplicit
%{
Gen_s = CreateABParameterGenerator({'l', 'cos', l, [1, 0, 0]; 'l', 'sin', l, [0, 1, 0]});
f = @(x, y, z) CalcABerr(x, y, z, zeros(size(x)), zeros(size(x)), zeros(size(x)), bm, l, Gen_s, GenV, "ax");
figure;
fimplicit3(f, [-10 10 -10 10 -10 10])
%}



%% Функции
function [t, v]= FindValue(f, targetValue, er, step_size, start_t, max_count_of_steps)
arguments
    f (1, 1) function_handle
    targetValue (1,1) double
    er (1,1) double
    step_size (1,1) double
    start_t (1,1) double
    max_count_of_steps (1,1) double = 1000
end
if targetValue <= 0
    error("targetValue должен быть больше 0.")
end

t_min = start_t;
t_max =start_t;

v1 = f(start_t);
v2 = 11;
step = 0;
while ~(v2 <= targetValue && v1 >= targetValue)
    if step > 2e4
        %error("Превышенно максимальное число шагов - "+num2str(max_count_of_steps)+".")
        xdata = 0:(step_size/10):start_t;
        ydata = zeros(1, numel(xdata));
        for i = 1:numel(xdata)
            ydata(i) = f(xdata(i));
        end
        plot(xdata, ydata)
        t = nan;
        v = nan;
        return
    end
    t_max = t_min;
    t_min = t_min-step_size;
    v2 = v1;
    v1 = f(t_min);
    step = step + 1;
    %disp("Шаг: " +num2str(step)+ " | t - [" + num2str(t_min) + ", " + num2str(t_max) + "] | v - [" + num2str(v1) + ", " + num2str(v2) + "]")%Отладка
end
step = 0;
tt = (t_max - t_min)/2;
vv = f(tt);

while abs(vv - targetValue) > er
    if step > max_count_of_steps
        %error("Превышенно максимальное число шагов - "+num2str(max_count_of_steps)+".")

        break
    end
    if vv > targetValue
        t_min = tt;
        tt = t_min+(t_max - t_min)/2;
    else
        t_max = tt;
        tt = t_min+(t_max - t_min)/2;
    end
    vv = f(tt);
    step = step + 1;
    %disp("Шаг: " +num2str(step)+ " | t - " + num2str(tt) + " | v - " + num2str(vv))%Отладка
end
t = tt;
v = vv;
end

function err = CalcAmplErr(x, y, z, alf, bet, gam, bm, l ,Gen_s, Gen_d)
arguments
    x (:,:) double
    y (:,:) double
    z (:,:) double
    alf (:,:) double
    bet (:,:) double
    gam (:,:) double
    bm (1,1) double
    l (1,1) double
    Gen_s (1, 1) function_handle
    Gen_d (1, 1) function_handle
end
if ~(numel(x) == numel(y) && numel(x) == numel(z) && numel(x) == numel(alf) && numel(x) == numel(bet) && numel(x) == numel(gam))
    error("Не подходящая размерность аргументов.")
end
norm_d = zeros(size(x));
for i = 1:numel(x)
    vd = Gen_d(bm, x(i), y(i), z(i), alf(i), bet(i), gam(i));
    vs = Gen_s(bm, x(i), y(i), z(i), alf(i), bet(i), gam(i));
    norm_d(i) = abs((sqrt(sum(vd.^2, 'all')) - sqrt(sum(vs.^2, 'all')))/sqrt(sum(vd.^2, 'all')))*100;
end
    err = norm_d;
end

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

sol_v = Gen_s(bm/l, alf, bet, gam, y, x, z);
err = 100*abs((CalcZ(bm, sol_v(1,1), sol_v(2,1), sol_v(3,1), sol_v(1,2), sol_v(2,2), sol_v(3,2)) - x)/sqrt(x.^2 + y.^2 + z.^2));

end
function abEr = CalcABerr(x, y, z, alf, bet, gam, bm, l, Gen_s, Gen_d, parameter)
arguments
    x (:,:) double
    y (:,:) double
    z (:,:) double
    alf (:,:) double
    bet (:,:) double
    gam (:,:) double
    bm (1,1) double
    l (1,1) double
    Gen_s (1, 1) function_handle
    Gen_d (1, 1) function_handle
    parameter {mustBeMember(parameter, ["ax", "ay", "az", "bx", "by", "bz"])}
end
if ~(numel(x) == numel(y) && numel(x) == numel(z) && numel(x) == numel(alf) && numel(x) == numel(bet) && numel(x) == numel(gam))
    error("Не подходящая размерность аргументов.")
end
erV = zeros(3, 2, numel(x));
for i = 1:numel(x)
    erV(:, :, i) = abs(Gen_d(bm, x(i), y(i), z(i), alf(i), bet(i), gam(i)) - Gen_s(bm, x(i), y(i), z(i), alf(i), bet(i), gam(i)));
end
switch parameter
    case "ax"
        abEr = erV(1, 1, :);
    case "ay"
        abEr = erV(2, 1, :);
    case "az"
        abEr = erV(3, 1, :);
    case "bx"
        abEr = erV(1, 2, :);
    case "by"
        abEr = erV(2, 2, :);
    case "bz"
        abEr = erV(3, 2, :);
end
end
