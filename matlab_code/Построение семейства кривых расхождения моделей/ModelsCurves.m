m = 1;
nL = 8;
mylog = @(a, b)log(b)./log(a);
p = 1;
var_L = ( linspace( (10e-2)^p, (1.5)^p, nL)).^(1/p);
mynorm = @(r) (sqrt(sum(r.^2)));

syms m0 l real positive
syms x y z real positive
r = [x,y,z]; % Радиус-вектор до точки
sym_ro = sym('ro', 'real'); %Плотность магнитного момента
rm = [l,0,0]; %Радиус-вектор до элементарного элменета катушки
assumptions
vec_ro = [sym_ro,0,0];
syms bx by bz b real
dB = (m0/(4*pi)) * (3*(r-rm)*dot(vec_ro,(r-rm))/mynorm(r-rm)^5-vec_ro/mynorm(r-rm)^3);
%Часть интегрирование Bx
fn = 1/((l-x)^2+(y)^2+(z)^2)^sym(3/2);
half_Bx1(l) = int(fn,l);
fn = (l-x)*(l-x)/(mynorm([l-x,y,z]))^5;
half_Bx2(l) = 3*int(fn,l);

NN = 100;

accros_fgX = figure('Name', 'accros_fgX');
acX_ax = axes(accros_fgX, 'YScale', 'log', 'XGrid', 'on', 'YGrid', 'on', 'NextPlot', 'add');
xlabel("y [м]"); ylabel("|\delta|B_x [%]");
accros_fgY = figure('Name', 'accros_fgY');
acY_ax = axes(accros_fgY, 'YScale', 'log', 'XGrid', 'on', 'YGrid', 'on', 'NextPlot', 'add');
xlabel("y [м]"); ylabel("|\delta|B_y [%]");
accros_fgZ = figure('Name', 'accros_fgZ');
acZ_ax = axes(accros_fgZ, 'YScale', 'log', 'XGrid', 'on', 'YGrid', 'on', 'NextPlot', 'add');
xlabel("y [м]"); ylabel("|\delta|B_z [%]");

along_fgX = figure('Name', 'along_fgX');
alX_ax = axes(along_fgX, 'YScale', 'log', 'XGrid', 'on', 'YGrid', 'on', 'NextPlot', 'add');
xlabel("x [м]"); ylabel("|\delta|B_x [%]");
along_fgY = figure('Name', 'along_fgY');
alY_ax = axes(along_fgY, 'YScale', 'log', 'XGrid', 'on', 'YGrid', 'on', 'NextPlot', 'add');
xlabel("x [м]"); ylabel("|\delta|B_y [%]");
along_fgZ = figure('Name', 'along_fgZ');
alZ_ax = axes(along_fgZ, 'YScale', 'log', 'XGrid', 'on', 'YGrid', 'on', 'NextPlot', 'add');
xlabel("x [м]"); ylabel("|\delta|B_z [%]");

disp("Магнитный момент равен: " + num2str(m)+ "[А*м^2]");
for i = 1:nL
    
    L = var_L(i);
    L0 = -L/2;
    ro = m/L;
    N = NN;
    
    disp(num2str(i) + " - ro = " + num2str(ro) + "[А*м]" + " | L = " + num2str(L) + "[м]") 
    
    %Заканчиваем интегрировать Bx
    Bx = (m0*sym_ro/(4*pi))*((half_Bx2(L0+L)-half_Bx2(L0))-(half_Bx1(L0+L)-half_Bx1(L0)));
    fn_Bx = matlabFunction(subs(simplify(Bx),m0,4e-7*pi));
    
    %Интегрирование By 
    By = simplify(subs(int(dB(2), l), l, L0+L) - subs(int(dB(2), l), l, L0));
    fn_By = matlabFunction(subs(simplify(By),m0,4e-7*pi));
    
    %Интегрирование Bz
    Bz = simplify(subs(int(dB(3),l),l,L0+L) - subs(int(dB(3),l),l,L0));
    fn_Bz = matlabFunction(subs(simplify(Bz),m0,4e-7*pi));
    
    
    %Получаем данные кривых поперёк
    N = 100;
    Y = linspace(10e-2, 111e-2, N);
    X = repmat(0.0001, size(Y));
    Z = repmat(0.0001, size(X));
    beacon = MagneticMoment(0, 0, 0, ro*L, 0, 0);
    [bf_x, bf_y, bf_z] = beacon.Field(X, Y, Z);
    bx_aprox = fn_Bx(repelem(ro,N), X, Y, Z);
    by_aprox = fn_By(repelem(ro,N), X, Y, Z);
    bz_aprox = fn_Bz(repelem(ro,N), X, Y, Z);
    
    %accros_xdata(i, :) = Y;
    accros_erbx = 100*abs( (bx_aprox - bf_x)./bx_aprox );
    accros_erby = 100*abs( (by_aprox - bf_y)./by_aprox );
    accros_erbz = 100*abs( (bz_aprox - bf_z)./bz_aprox );
    
    %acX_centers(:,i) = [Y(round(N/2)), accros_erbx(round(N/2))]' ;
    ind = 50;%N/2
    plot(acX_ax, Y, accros_erbx, "black");
    %[tmp, ind] = min(abs(accros_erbx-2));
    text(acX_ax, Y(ind), accros_erbx(ind), num2str(i), 'BackgroundColor', 'white', 'FontSize', 5, 'Margin', 1)
    %acY_centers(:,i) = [Y(round(N/2)), accros_erby(round(N/2))]' ;
    plot(acY_ax, Y, accros_erby, "black");
    %[tmp, ind] = min(abs(accros_erby-2));
    text(acY_ax, Y(ind), accros_erby(ind), num2str(i), 'BackgroundColor', 'white', 'FontSize',5, 'Margin', 1)
    %acZ_centers(:,i) = [Y(round(N/2)), accros_erbz(round(N/2))]' ;
    plot(acZ_ax, Y, accros_erbz, "black");
    %[tmp, ind] = min(abs(accros_erbz-2));
    text(acZ_ax, Y(ind), accros_erbz(ind), num2str(i), 'BackgroundColor', 'white', 'FontSize',5, 'Margin', 1)
    
    
    %Получаем данные кривых вдоль
    X = (L0 + 1.45*L):1e-2:10*L;
    Y = repmat(0.0001, size(X));
    Z = repmat(0.0001, size(X));
    [bf_x, bf_y, bf_z] = beacon.Field(X, Y, Z);
    N = size(X,2);
    bx_aprox = fn_Bx(repelem(ro,N), X, Y, Z);
    by_aprox = fn_By(repelem(ro,N), X, Y, Z);
    bz_aprox = fn_Bz(repelem(ro,N), X, Y, Z);
    
    %along_xdata = X;
    along_erbx = 100*abs( (bx_aprox - bf_x)./bx_aprox );
    along_erby = 100*abs( (by_aprox - bf_y)./by_aprox );
    along_erbz = 100*abs( (bz_aprox - bf_z)./bz_aprox );
    
    %alX_centers(:,i) = [X(round(N/2)), along_erbx(round(N/2))]' ;
    plot(alX_ax, X, along_erbx, "black");
    [tmp, ind] = min(abs(along_erbx-2));
    text(alX_ax, X(ind), along_erbx(ind), num2str(i), 'BackgroundColor', 'white', 'FontSize', 5, 'Margin', 1)
    %alY_centers(:,i) = [X(round(N/2)), along_erby(round(N/2))]' ;
    plot(alY_ax, X, along_erby, "black");
    [tmp, ind] = min(abs(along_erby-2));
    text(alY_ax, X(ind), along_erby(ind), num2str(i), 'BackgroundColor', 'white', 'FontSize', 5, 'Margin', 1)
    %alZ_centers(:,i) = [X(round(N/2)), along_erbz(round(N/2))]' ;
    plot(alZ_ax, X, along_erbz, "black");
     [tmp, ind] = min(abs(along_erbz-2));
    text(alZ_ax, X(ind), along_erbz(ind), num2str(i), 'BackgroundColor', 'white', 'FontSize', 5, 'Margin', 1)

end

%{
accros_fgX.Children.XLabel.Text.String = {'y [м]'};
accros_fgX.Children.YLabel.Text.String = {'|\delta|B_x [%]'};
accros_fgY.Children.XLabel.Text.String = {'y [м]'};
accros_fgY.Children.YLabel.Text.String = {'|\delta|B_y [%]'};
accros_fgZ.Children.XLabel.Text.String = {'y [м]'};
accros_fgZ.Children.YLabel.Text.String = {'|\delta|B_z [%]'};


along_fgX.Children.XLabel.Text.String = {'x [м]'};
along_fgX.Children.YLabel.Text.String = {'|\delta|B_x [%]'};
along_fgY.Children.XLabel.Text.String = {'x [м]'};
along_fgY.Children.YLabel.Text.String = {'|\delta|B_y [%]'};
along_fgZ.Children.XLabel.Text.String = {'x [м]'};
along_fgZ.Children.YLabel.Text.String = {'|\delta|B_z[%]'};
%}

%% Не работает
for i = 1:numel(nL)
    for it = {{acX_ax, acX_centers}, {acY_ax, acY_centers}, {acZ_ax, acZ_centers}, {alX_ax, alX_centers}, {alY_ax, alY_centers}, {alZ_ax, alZ_centers}}
        d = min(abs([it{1}{1}.XLim(2)-it{1}{1}.XLim(1), it{1}{1}.YLim(2)-it{1}{1}.YLim(1)]))/200;
        n = 10;
        doub = it{1}{2};        
        x = (doub(1,i)+d/2).*cos(linspace(0, 2*pi, n));
        y = (doub(2,i)+d/2).*sin(linspace(0, 2*pi, n));
        
        fill(it{1}{1}, x, y, "white");
        plot(it{1}{1}, x, y, "black");
        text(it{1}{1}, doub(1,i), doub(2,i), num2str(i), "HorizontalAlignment", "center" , 'FontSize',7, "VerticalAlignment", "middle");
    end
end