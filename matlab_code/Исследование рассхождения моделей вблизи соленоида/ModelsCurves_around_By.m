m = 1;
mylog = @(a, b)log(b)./log(a);
p = 1;
var_L = 0.25;
mynorm = @(r) (sqrt(sum(r.^2)));
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
%}

disp("Магнитный момент равен: " + num2str(m)+ "[А*м^2] | Длинна соленоида равна: " + num2str(var_L));
angles_n = 100;
alf = deg2rad( linspace(0.01, 89.99, angles_n));
bet = deg2rad( linspace(0.01, 89.99, angles_n));

erbcord = zeros(3, angles_n, angles_n);
erer = 0.01;
L = var_L;
L0 = -L/2;
ro = m/L;

%{
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
%}
beacon = MagneticMoment(0, 0, 0, m, 0, 0);
lin = MagneticLine(0, 0, var_L, [0,0,0], m/L);
step = 0.001;
nrm = lin.Field(1, 1, 1);
parfor i = 1:angles_n
    for j = 1:angles_n
        
        %Метод половинного деления
        
        rm = 10;
        rM =10;        
        X = rm*cos(bet(j))*cos(alf(i));
        Y = rm*cos(bet(j))*sin(alf(i));
        Z = rm*sin(bet(j));
        [bf_x, bf_y, bf_z] = beacon.Field(X, Y, Z);
        [b1, b_aprox, b2] = lin.Field(X, Y, Z);   %b_aprox = fn_By(ro, X, Y, Z);
        er1 = 100*abs( (b_aprox - bf_y)./nrm );
        er2 = 11;

        while ~(er2 <= 10 && er1 >= 10)
            
            rM = rm;
            rm = rm-step;
            er2 = er1;
            X = rm*cos(bet(j))*cos(alf(i));
            Y = rm*cos(bet(j))*sin(alf(i));
            Z = rm*sin(bet(j));
            [bf_x, bf_y, bf_z] = beacon.Field(X, Y, Z);
            [b1, b_aprox, b2] = lin.Field(X, Y, Z);   %b_aprox = fn_By(ro, X, Y, Z);
            er1 = 100*abs( (b_aprox - bf_y)./nrm );
        end

        r = rm + (rM-rm)/2;
        
        X = r*cos(bet(j))*cos(alf(i));
        Y = r*cos(bet(j))*sin(alf(i));
        Z = r*sin(bet(j));        
        [bf_x, bf_y, bf_z] = beacon.Field(X, Y, Z);
        [b1, b_aprox, b2] = lin.Field(X, Y, Z);   %b_aprox = fn_By(ro, X, Y, Z);
        erb = 100*abs( (b_aprox - bf_y)./nrm ); 

        while abs(erb-10) > erer
            if erb > 10
                rm = r;
                r = rm+(rM - rm)/2;                
            else
                rM = r;
                r = rm+(rM - rm)/2;
            end
            
            X = r*cos(bet(j))*cos(alf(i));
            Y = r*cos(bet(j))*sin(alf(i));
            Z = r*sin(bet(j));
            [b1, b_aprox, b2] = lin.Field(X, Y, Z);   %b_aprox = fn_By(ro, X, Y, Z);
            [bf_x, bf_y, bf_z] = beacon.Field(X, Y, Z);
            erb = 100*abs( (b_aprox - bf_y)./nrm );
            %disp("erb = " + num2str(erb) + "rm = " + num2str(rm) + "rM = " + num2str(rM));
        end

        erbcord(:, i, j) = [X, Y, Z];
        
    end
    disp(i);
end

%%
fg1 = figure;
Xdata = zeros(angles_n, angles_n);
Ydata = zeros(angles_n, angles_n);
Zdata = zeros(angles_n, angles_n);
absErB = zeros(angles_n, angles_n);
for i = 1:angles_n
    for j = 1:angles_n  
        b_aprox = fn_By(ro, erbcord(1, i, j), erbcord(2, i, j), erbcord(3, i, j));
        [bf_x, bf_y, bf_z] = beacon.Field(erbcord(1, i, j), erbcord(2, i, j), erbcord(3, i, j));
        absErB(i,j) = 100*abs( (b_aprox - bf_y) );
        
        %plot3(erbcord(1,:, j), erbcord(2,:, j), erbcord(3,:, j))
        Xdata(i, j) = erbcord(1, i, j);
        Ydata(i, j) = erbcord(2, i, j);
        Zdata(i, j) = erbcord(3, i, j);
    end
end
%surface(Xdata, Ydata, Zdata, "EdgeColor", 'none')
title("B_y L = 0.25m") 
disp("Среднее ошибок расхождения моделей по By = " + num2str(mean(absErB,  'all')))
disp("Стандартное отклонение моделей по By = " + num2str(std(absErB, [],  'all')))
%{
fg2 = figure; fimplicit3
for j = 1:angles_n    
    plot3(erbycord(1,:, j), erbycord(2,:, j), erbycord(3,:, j))
    if j == 1
        hold on
    end
end
%}