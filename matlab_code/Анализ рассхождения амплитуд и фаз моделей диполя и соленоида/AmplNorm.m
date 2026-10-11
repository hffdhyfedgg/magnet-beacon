bm = 1e-6;
l = 0.25;
Gen_s = CreateABParameterGenerator({'l', 'cos', l, [1, 0, 0]; 'l', 'sin', l, [0, 1, 0]});

step = 0.001;
erer = 0.1;

angles_n = 100;
alf = linspace(1e-5, pi/2 - 1e-5, angles_n);
bet = linspace(1e-5, pi/2 - 1e-5, angles_n);
n = 50;
y = linspace(0, 1, n);
erbcord_x = zeros(angles_n, angles_n);
erbcord_y = zeros(angles_n, angles_n);
erbcord_z = zeros(angles_n, angles_n);


parfor i = 1:angles_n
    xx = zeros(1, angles_n);
    yy = zeros(1, angles_n);
    zz = zeros(1, angles_n);
    for j = 1:angles_n
        %{
        f = @(t) CalcAmplEr(t*cos(alf(i))*cos(bet(j)), t*, t*sin(bet(i)), 1e-3, 1e-3, 1e-3, bm, l, GenV, Gen_s );
        [r, v] = FindValue(f, 10, erer, step, 25, 1e7);
        xx(j) = r*cos(bet(i));
        yy(j) = y(j);
        zz(j) = r*sin(bet(i));
        %erV(i, j) = sqrt(sum((Gen_s(bm/l, 0, 0, 0, yy(j), xx(j), zz(j)) - GenV(bm, 0, 0, 0, yy(j), xx(j), zz(j))).^2, 'all'));
         %}
        
        t = linspace(1e-3, 1, 1e5);
        result_error = CalcAmplEr(t*cos(alf(i))*cos(bet(j)), t*sin(alf(i))*cos(bet(j)), t*sin(bet(j)), repmat(1e-4, size(t)), repmat(1e-4, size(t)), repmat(1e-4, size(t)), bm, l, GenV, Gen_s);
        [val, index] = max(result_error);
        [val, index] = min(abs(result_error(index:end)-10));
        %{
        semilogy(t, result_error)
        hold on
        semilogy(t(index), result_error(index), 'x')
        hold off
        %}
        xx(j) = t(index)*cos(alf(i))*cos(bet(j));
        yy(j) = t(index)*sin(alf(i))*cos(bet(j));
        zz(j) = t(index)*sin(bet(j));
    end
    erbcord_x(i, :) = xx;
    erbcord_y(i, :) = yy;
    erbcord_z(i, :) = zz;
    disp(i);
end

%%
fg1 = figure;
surface(erbcord_x, erbcord_y, erbcord_z, "EdgeColor", 'none')
%xlim([0, max(erbcord_x, [], 'all')]); ylim([0, max(erbcord_y, [], 'all')]); zlim([0, max(erbcord_z, [], 'all')]);
title("Поверхность превышения погрешности в 10 % по норме амплитуд")
xlabel('x')




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

t_min = 0;
t_max = 0;        

v1 = targetValue*0.9;
v2 = f(start_t);
step = 0;
while ~(v2 <= targetValue && v1 >= targetValue)
    %{
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
    %}
    t_min = t_max;
    t_max = t_max + step_size;
    v1 = v2;
    v2 = f(t_max);
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

function err = CalcAmplEr(x, y, z, alf, bet, gam, bm, l, gen_d, gen_s)
arguments
    x  double
    y  double
    z  double
    alf  double
    bet  double
    gam  double
    bm (1,1) double
    l (1,1) double
    gen_d (1, 1) function_handle
    gen_s (1, 1) function_handle
end

if ~(numel(x) == numel(y) && numel(x) == numel(z) && numel(x) == numel(alf) && numel(x) == numel(bet) && numel(x) == numel(gam))
    error("Не подходящая размерность аргументов.")
end
dip = zeros(3,2);
sol = zeros(3,2);
ampl_dip = zeros(3,2);
ampl_sol = zeros(3,2);
abs_er = zeros(size(x));
rel_er = zeros(size(x));

for i = 1:numel(x)
    dip = gen_d(bm, alf(i), bet(i), gam(i), y(i),  x(i), z(i));
    sol = gen_s(bm/l, alf(i), bet(i), gam(i),y(i),  x(i), z(i));
    ampl_dip = AmplAngle(dip(1), dip(2), dip(3), dip(4), dip(5), dip(6));
    ampl_sol = AmplAngle(sol(1), sol(2), sol(3), sol(4), sol(5), sol(6));
    abs_er(i) = sqrt(sum([ampl_dip{:,1}].^2)) - sqrt(sum([ampl_sol{:,1}].^2));
    rel_er(i) = 100*abs(abs_er(i)./sqrt(sum([ampl_dip{:,1}].^2)));
end
err = rel_er;
end

function amplPhase = AmplAngle(ax, ay, az, bx, by, bz)
aamplx = zeros(size(ax));
pphx = zeros(size(ax));
aamply = zeros(size(ax));
pphy = zeros(size(ax));
aamplz = zeros(size(ax));
pphz = zeros(size(ax));

for i = 1:numel(ax)
    aamplx(i) = sqrt(ax(i)^2 + bx(i)^2);
    pphx(i) = atan(bx(i)/ax(i));
    aamply(i) = sqrt(ay(i)^2 + by(i)^2);
    pphy(i) = atan(by(i)/ay(i));
    aamplz(i) = sqrt(az(i)^2 + bz(i)^2);
    pphz(i) = atan(bz(i)/az(i));
end
amplPhase = {aamplx, pphx; aamply, pphy; aamplz, pphz};
end