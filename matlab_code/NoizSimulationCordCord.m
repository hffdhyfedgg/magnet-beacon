%function [out] = NoizSimulation(sampl_size)

%Зададим углы и координаты
%Ток на катушке - 0.5 А, Сопротивление на 1кГц - 250Ом
%clear
%other_r = @(b, ax, ay, az, bx, by, bz) 2^(1/6)*(1/(-sqrt((ax*bx+ay*by+az*bz)^2*4+(ax^2+ay^2+az^2-bx^2-by^2-bz^2)^2)+ax^2+ay^2+az^2+bx^2+by^2+bz^2))^(1/6)*(b^2)^(1/6);

load D:\Bezmen\PavelSheptuhin\GenVS.mat
N =100;
sampl_size = 1e5;
%f = 1e3; %частота
%fs = 1e6;
%T = 1/f;
x = linspace( 0.1, 10, N);
y = linspace( 0.1, 10, N);
z = linspace( 0.1, 10, N);
gam = repelem( deg2rad(0), N);
bet = repelem( deg2rad(0), N);
alf =   repelem( deg2rad(0), N);
Bm = 2*1e-7*700*pi*(23e-3/2)^2*0.5;%700 витков, диаметр 23 мм, ток 0.5А;

rr1 = zeros(size(linspace(log(0.5*1e-9/22.42), log(300*1e-9/22.42), 10)));
rr2 = zeros(size(linspace(log(0.5*1e-9/22.42), log(300*1e-9/22.42), 10)));
sig = 0;%exp(linspace(log(0.5*1e-9/22.42), log(300*1e-9/22.42), 10));
%sig = 1e-9/22.42;


for j = 1:numel(sig)
tag = "otherSigm";


true_V = zeros(3, 2, N);

for i = 1:N
    
    true_V(:,:,i) =  GenVS(Bm/0.25, alf(i), bet(i), gam(i), y(i), x(i), z(i));
end

aa = zeros(3, N, sampl_size);
bb = zeros(3, N, sampl_size);

sigm = sig(j);


parfor k = 1:sampl_size
    
    %V = true_V(1, 1, :);
    for i = 1:N        
        %aa(1,i,k) = true_V(1,1,i);   aa(2,i,k) = true_V(2,1,i) ;   aa(3,i,k) = true_V(3,1,i);
        %bb(1,i,k) = true_V(1,2,i);   bb(2,i,k) = true_V(2,2,i);   bb(3,i,k) = true_V(3,2,i);
        Ax  = normrnd(sqrt(true_V(1, 1, i)^2 + true_V(1, 2, i)^2), sigm); 
        phx = normrnd(atan(true_V(1, 2, i)/true_V(1, 1, i)), sigm);
        
        Ay  = normrnd(sqrt(true_V(2, 1, i)^2 + true_V(2, 2, i)^2), sigm); 
        phy = normrnd(atan(true_V(2, 2, i)/true_V(2, 1, i)), sigm); 
        
        Az  = normrnd(sqrt(true_V(3, 1, i)^2 + true_V(3, 2, i)^2), sigm); 
        phz = normrnd(atan(true_V(3, 2, i)/true_V(3, 1, i)), sigm);
        
        aa(:, i, k) = [Ax*cos(phx), Ay*cos(phy), Az*cos(phz)];
        bb(:, i, k) = [Ax*sin(phx), Ay*sin(phy), Az*sin(phz)];      
 %       f(i,k) = parfeval(@AddNoize, 2, true_V(:,:,i), sigm);
    end
   % aa(:,:,k) = a;
    %bb(:,:,k) = b;
    %disp("Паралельное вычисление ab параметров: выполнен шаг №"+num2str(k)+". Всего шагов: "+num2str(sampl_size)+" Время завершения шага: "+string(datetime))
end




%toc


%{
disp("Последовательное вычисление ab параметров")
tic
a = zeros(3,N);
b = zeros(3,N);
for k = 1:sampl_size
    
    for i = 1:N
        
        %aa(1,i,k) = true_V(1,1,i);   aa(2,i,k) = true_V(2,1,i) ;   aa(3,i,k) = true_V(3,1,i);
        %bb(1,i,k) = true_V(1,2,i);   bb(2,i,k) = true_V(2,2,i);   bb(3,i,k) = true_V(3,2,i);
        Ax  = normrnd(sqrt(true_V(1,1,i)^2 + true_V(1,2,i)^2), sigm); 
        phx = normrnd(atan(true_V(1,2,i)/true_V(1,1,i)), sigm);
        
        Ay  = normrnd(sqrt(true_V(2,1,i)^2 + true_V(2,2,i)^2), sigm); 
        phy = normrnd(atan(true_V(2,2,i)/true_V(2,1,i)), sigm); 
        
        Az  = normrnd(sqrt(true_V(3,1,i)^2 + true_V(3,2,i)^2), sigm); 
        phz = normrnd(atan(true_V(3,2,i)/true_V(3,1,i)), sigm);
        
        aa(:, i, k) = [Ax*cos(phx), Ay*cos(phy), Az*cos(phz)];
        bb(:, i, k) = [Ax*sin(phx), Ay*sin(phy), Az*sin(phz)];   
    end

end
toc
%}
%{
err.r = CalcR(repelem(Bm,N), mean(aa(1,:,1),3), mean(aa(2,:,1),3), mean(aa(3,:,1),3), mean(bb(1,:,1),3), mean(bb(2,:,1),3), mean(bb(3,:,1),3));
err.x = CalcX(repelem(Bm,N), mean(aa(1,:,1),3), mean(aa(2,:,1),3), mean(aa(3,:,1),3), mean(bb(1,:,1),3), mean(bb(2,:,1),3), mean(bb(3,:,1),3));
err.y = CalcY(repelem(Bm,N), mean(aa(1,:,1),3), mean(aa(2,:,1),3), mean(aa(3,:,1),3), mean(bb(1,:,1),3), mean(bb(2,:,1),3), mean(bb(3,:,1),3));
err.z = CalcZ(repelem(Bm,N), mean(aa(1,:,1),3), mean(aa(2,:,1),3), mean(aa(3,:,1),3), mean(bb(1,:,1),3), mean(bb(2,:,1),3), mean(bb(3,:,1),3));
err.z1 = sqrt(err.r.^2-err.x.^2-err.y.^2);
err.z2 = err.r.^2-err.x.^2-err.y.^2;
%}

%% ab Для мнимых значений
%{
for i = 1:numel(err.z)
    if imag(err.z(i)) ~= 0
        imag(err.z(i));
        %ima = [aa(1,i,1), aa(2,i,1), aa(3,i,1)];
        %imb = [bb(1,i,1), bb(2,i,1), bb(3,i,1)];
    end
end
%}
%% Координаты
errR = zeros(1,N);   abserrR = zeros(1,N);
errX = zeros(1,N);   abserrX = zeros(1,N);
errY = zeros(1,N);   abserrY = zeros(1,N);
errZ = zeros(1,N);   abserrZ = zeros(1,N);
errAlf = zeros(1,N);   abserrAlf = zeros(1,N);
errBet = zeros(1,N);   abserrBet = zeros(1,N);
errGam = zeros(1,N);   abserrGam = zeros(1,N);
erX = zeros(1, N);
erY = zeros(1, N);
erZ = zeros(1, N);
erR = zeros(1, N);
var = zeros(1, N);
fim1 = zeros(1, N); %Вероятность, что мнимая часть привысит реальную.
fim2 = zeros(1, N); %Вероятность появления мнимой чати.
%tic

flag1 = true;
flag2 = true;

for i = 1:N

    a = zeros(3,sampl_size);
    b = zeros(3,sampl_size);
    var_a1 = aa(1,i,:);
    var_a2 = aa(2,i,:);
    var_a3 = aa(3,i,:);
    var_b1 = bb(1,i,:);
    var_b2 = bb(2,i,:);
    var_b3 = bb(3,i,:);
    trueR = sqrt(x(i)^2+y(i)^2+z(i)^2);
    trueX = x(i);
    trueY = y(i);
    trueZ = z(i);
    trueAlf = repelem(alf(i),sampl_size);
    trueBet = repelem(bet(i),sampl_size);
    trueGam = repelem(gam(i),sampl_size);
    
    R_abserr = zeros(1,sampl_size);
    R_err = zeros(1,sampl_size);
    R_er = zeros(1,sampl_size);
    
    X_abserr = zeros(1,sampl_size);
    X_err = zeros(1,sampl_size);
    X_er = zeros(1,sampl_size);
    
    Y_abserr = zeros(1,sampl_size);
    Y_err = zeros(1,sampl_size);
    Y_er = zeros(1,sampl_size);
    
    Z_abserr = zeros(1,sampl_size);
    Z_err = zeros(1,sampl_size);
    Z_er = zeros(1,sampl_size);
    
    er_alf = zeros(1,sampl_size);
    er_bet = zeros(1,sampl_size);
    er_gam = zeros(1,sampl_size);
    
    ffim2 = 0;
    ffim1 = 0;
    
    parfor k = 1:sampl_size
        %a(:,k) = [var_a1(k), var_a2(k), var_a3(k)];
        %b(:,k) = [var_b1(k), var_b2(k), var_b3(k)];
        R_er(k) = CalcR(Bm, var_a1(k), var_a2(k), var_a3(k), var_b1(k), var_b2(k), var_b3(k));
        R_abserr(k) = abs(real(R_er(k)) - trueR);
        R_err(k) = 100*abs((R_er(k)-trueR)./trueR);
        
        X_er(k) = CalcX(Bm, var_a1(k), var_a2(k), var_a3(k), var_b1(k), var_b2(k), var_b3(k));
        X_abserr(k) = abs(real(X_er(k))-trueX);
        X_err(k) = 100*abs((real(X_er(k))-trueX)./trueX);
        
        Y_er(k) = CalcY(Bm, var_a1(k), var_a2(k), var_a3(k), var_b1(k), var_b2(k), var_b3(k));
        Y_abserr(k) = abs(real(Y_er(k))-trueY);
        Y_err(k) = 100*abs((real(Y_er(k))-trueY)./trueY);
        
        Z_er(k) = CalcZ(Bm, var_a1(k), var_a2(k), var_a3(k), var_b1(k), var_b2(k), var_b3(k));
        Z_abserr(k) = abs(real(Z_er(k)) - trueZ);
        Z_err(k) = 100*abs((real(Z_er(k)) - trueZ)./trueZ);        
        
        if imag(Z_er(k))   ~= 0
            ffim2 = ffim2+1;
        end
        
        if real(Z_er(k)) > imag(Z_er(k))
            ffim1 = ffim1+1;
        end
        %rotM = CalcRotMatrix(Bm, var_a1(k), var_a2(k), var_a3(k), var_b1(k), var_b2(k), var_b3(k));
        %[er_alf(k), er_bet(k), er_gam(k)] = CalcAngles(real(rotM(:,:)));
    end
    %var(i) = sqrt(true_V(1,1,i).^2+true_V(2,1,i).^2+true_V(3,1,i).^2).*sqrt(true_V(1,2,i).^2+true_V(2,2,i).^2+true_V(3,2,i).^2)-(true_V(1,1,i).*true_V(1,2,i) + true_V(2,1,i).*true_V(2,2,i) + true_V(3,1,i).*true_V(3,2,i))
    %var(i) = CalcX(Bm, true_V(1,1,i), true_V(2,1,i), true_V(3,1,i), true_V(1,2,i), true_V(2,2,i), true_V(3,2,i))./CalcR(Bm, true_V(1,1,i), true_V(2,1,i), true_V(3,1,i), true_V(1,2,i), true_V(2,2,i), true_V(3,2,i)) - sqrt( (CalcR(Bm, true_V(1,1,i), true_V(2,1,i), true_V(3,1,i), true_V(1,2,i), true_V(2,2,i), true_V(3,2,i))^6 * (true_V(1,1,i).^2+true_V(2,1,i).^2+true_V(3,1,i).^2))/(3*Bm^2) - 1/3);
    %{
    var1 = CalcR(repelem(Bm,N),a(1,:),a(2,:),a(3,:),b(1,:),b(2,:),b(3,:));
    var2 = sqrt(x.^2+y.^2+z.^2);
    er = mean(100*abs((var1-var2)./var2));
    %}
    
    
    fim1(i) = ffim1/sampl_size;
    fim2(i) = ffim2/sampl_size;
    
    abserrR(i) =  mean(R_abserr);
    errR(i) = mean(R_err);
    erR(i) = mean(R_er);

    abserrX(i) =  mean(X_abserr);
    errX(i) = mean(X_err);
    erX(i) = mean(X_er);

    abserrY(i) =  mean(Y_abserr);
    errY(i) = mean(Y_err);
    erY(i) = mean(Y_er);

    abserrZ(i) =  mean(Z_abserr);
    errZ(i) = mean(Z_err);
    erZ(i) = mean(Z_er);
   %Частота появления мнимой части
    
    if i ~= 1
        if flag1
            if ( real(erZ(i-1) ) >= imag( erZ(i-1)) ) && ( real(erZ(i) ) < imag( erZ(i)) )
                flag1 = false;
                rr1(j) = trueR;

            end
        end

        if flag2
            if imag( erZ(i) ) ~= 0
                flag2 = false;
                rr2(j) = trueR;

            end
        end
    end
    abserrAlf(i) = mean(abs(er_alf-trueAlf));
    errAlf(i) = mean(100*abs((er_alf-trueAlf)./trueAlf));

    abserrBet(i) = mean(abs(er_bet-trueBet));
    errBet(i) = mean(100*abs((er_bet-trueBet)./trueBet));

    abserrGam(i) = mean(abs(er_gam - trueGam));
    errGam(i) = mean(100*abs((er_gam-trueGam)./trueGam));
    
    disp("Паралельное вычисление координат и погрешностей параметров. Выполнено: "+num2str(100*i/N)+"%"+" Время завершения шага: "+string(datetime))
end
%{
Графики в-ей появления или превышения мнимой части
figure;
title( "Dispertion = "+num2str(sigm) + ", вероятности")
hold on
plot(sqrt(x.^2+y.^2+z.^2), fim1, "blue", sqrt(x.^2+y.^2+z.^2), fim2, "red");
plot([rr1(j) rr1(j)], [0 1], "--black");
plot([rr2(j) rr2(j)], [0 1], ":black");
legend("Реальная часть больше мнимой", "Появление мнимой части.", "Превышение реальной части", "Появление мнимой части");
hold off
%}

if j == 1
% Графики измерения координат


saveFlag = false;
xdata = sqrt(x.^2+y.^2+z.^2);
inftag = "solenoid";
xlabeltxt = "r [м]";
fg = figure( "Name", "Dispertion = "+num2str(sigm) );
hold on
%title("sigm = " + num2str(sigm))
title("Соленоид")

plot(xdata, real(erZ), "blue", xdata, imag(erZ), "red", xdata, z, "black" );
plot([rr1(j) rr1(j)], [min([real(erZ), imag(erZ)]), max([real(erZ), imag(erZ)])], "--black");
plot([rr2(j) rr2(j)], [min([real(erZ), imag(erZ)]), max([real(erZ), imag(erZ)])], ":black");
xlabel("r [м]"); ylabel("z [м]");
legend("Реальная часть расч. Z", "Мнимая часть расч. Z", "Истинная Z")
hold off
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+".png")
end

fg = figure;
title("Относительная погрешность R")
zdata = errR;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\delta|R [%]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Относительная погрешность R.png")
end

fg = figure;
title("Абсолютная погрешность R")
zdata = abserrR;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\Delta|R [м]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Абсолютная погрешность R.png")
end

fg = figure;
title("Относительная погрешность X")
plot(xdata, errX,...
     [rr1(j) rr1(j)], [min(errX), max(errX)], "--black",...
     [rr2(j) rr2(j)], [min(errX), max(errX)], ":black");
xlabel(xlabeltxt); ylabel("|\delta|X [%]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Относительная погрешность X.png")
end

fg = figure;
title("Абсолютная погрешность X")
zdata = abserrX;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\Delta|X [м]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Абсолютная погрешность X.png")
end

fg = figure;
title("Относительная погрешность Y")
plot(xdata, errY,...
     [rr1(j) rr1(j)], [min(errY), max(errY)], "--black",...
     [rr2(j) rr2(j)], [min(errY), max(errY)], ":black");
xlabel(xlabeltxt); ylabel("|\delta|Y [%]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Относительная погрешность Y.png")
end

fg = figure;
title("Абсолютная погрешность Y")
zdata = abserrY;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\Delta|Y [м]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Абсолютная погрешность Y.png")
end

fg = figure;
title("Относительная погрешность Z")
plot(xdata, errZ,...
     [rr1(j) rr1(j)], [min(errZ), max(errZ)], "--black",...
     [rr2(j) rr2(j)], [min(errZ), max(errZ)], ":black");
xlabel(xlabeltxt); ylabel("|\delta|Z [%]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Относительная погрешность Z.png")
end

fg = figure;
title("Абсолютная погрешность Z")
zdata = abserrZ;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\Delta|Z [м]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Абсолютная погрешность Z.png")
end

end

clear R_abserr R_err R_er X_abserr X_err X_er Y_abserr Y_err Y_er Z_abserr Z_err Z_er var_a1 var_a2 var_a3 var_b1 var_b2 var_b3




end

fg = figure;
title("Граничный радиус")
plot(Bm./sig, rr1, "blue", Bm./sig, rr2, "red");
xlabel("Bm/sig"); ylabel("r [м]");
legend("Превышение реальной части z", "Появление мнимой части z");

saveas(fg, "гран.png")
%{
%% Графики измерения координат


saveFlag = true;
xdata = alf;
inftag = "x4y4z4";
xlabeltxt = "\alpha [rad]";
fg = figure( "Name", "Dispertion = "+num2str(sigm) );
hold on
title("sigm = " + num2str(sigm))

plot(xdata, real(erZ), "blue", xdata, imag(erZ), "red", xdata, z, "black" );
plot([rr1(j) rr1(j)], [min([real(erZ), imag(erZ)]), max([real(erZ), imag(erZ)])], "--black");
plot([rr2(j) rr2(j)], [min([real(erZ), imag(erZ)]), max([real(erZ), imag(erZ)])], ":black");
xlabel("r [м]"); ylabel("z [м]");
legend("Реальная часть расч. Z", "Мнимая часть расч. Z", "Истинная Z")
hold off
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+".png")
end

fg = figure;
title("Относительная погрешность R")
zdata = errR;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\delta|R [%]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Относительная погрешность R.png")
end

fg = figure;
title("Абсолютная погрешность R")
zdata = abserrR;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\Delta|R [м]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Абсолютная погрешность R.png")
end

fg = figure;
title("Относительная погрешность X")
plot(xdata, errX,...
     [rr1(j) rr1(j)], [min(errX), max(errX)], "--black",...
     [rr2(j) rr2(j)], [min(errX), max(errX)], ":black");
xlabel(xlabeltxt); ylabel("|\delta|X [%]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Относительная погрешность X.png")
end

fg = figure;
title("Абсолютная погрешность X")
zdata = abserrX;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\Delta|X [м]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Абсолютная погрешность X.png")
end

fg = figure;
title("Относительная погрешность Y")
plot(xdata, errY,...
     [rr1(j) rr1(j)], [min(errY), max(errY)], "--black",...
     [rr2(j) rr2(j)], [min(errY), max(errY)], ":black");
xlabel(xlabeltxt); ylabel("|\delta|Y [%]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Относительная погрешность Y.png")
end

fg = figure;
title("Абсолютная погрешность Y")
zdata = abserrY;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\Delta|Y [м]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Абсолютная погрешность Y.png")
end

fg = figure;
title("Относительная погрешность Z")
plot(xdata, errZ,...
     [rr1(j) rr1(j)], [min(errZ), max(errZ)], "--black",...
     [rr2(j) rr2(j)], [min(errZ), max(errZ)], ":black");
xlabel(xlabeltxt); ylabel("|\delta|Z [%]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Относительная погрешность Z.png")
end

fg = figure;
title("Абсолютная погрешность Z")
zdata = abserrZ;
plot(xdata, zdata,...
     [rr1(j) rr1(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], "--black",...
     [rr2(j) rr2(j)], [min([real(zdata), imag(zdata)]), max([real(zdata), imag(zdata)])], ":black");
xlabel(xlabeltxt); ylabel("|\Delta|Z [м]");
if saveFlag 
    saveas(fg, inftag+" dispertion_"+num2str(sigm)+"Абсолютная погрешность Z.png")
end




%% Графики углов

figure(fg(9))
plot(axes(fg(9)), xdata, errAlf);
xlabel(xlabeltxt); ylabel("|\delta| \alpha [%]");

figure(fg(10))
plot(axes(fg(10)), xdata, rad2deg(abserrAlf));
xlabel(xlabeltxt); ylabel("|\Delta| alpha [rad]");

figure(fg(11))
plot(axes(fg(11)), xdata, errBet);
xlabel(xlabeltxt); ylabel("|\delta| \beta [%]");
figure(fg(12))
plot(axes(fg(12)), xdata, rad2deg(abserrBet));
xlabel(xlabeltxt); ylabel("|\Delta| \beta [rad]");
figure(fg(13))
plot(axes(fg(13)), xdata, errGam);
xlabel(xlabeltxt); ylabel("|\delta|\gamma [%]");
figure(fg(14))
plot(axes(fg(14)), xdata, rad2deg(abserrGam));
xlabel(xlabeltxt); ylabel("|\Delta|\gamma [grad]");

%% сохранение
if false
    k = 1;
    for i = fg        
        saveas(i, "figure_" + fg(k).Name + ".png");
        k = k+1;        
    end
end 
%plot(sqrt(x.^2+y.^2+z.^2), sqrt(x.^2+y.^2+z.^2), "blue", sqrt(x.^2+y.^2+z.^2), CalcR(repelem(Bm,N), a(1,:), a(2,:), a(3,:), b(1,:), b(2,:), b(3,:)));
%=normrnd()
out = fg;
%end
%%


%{
function [aa, bb] = AddNoize( true_V, sigm )
        Ax  = normrnd(sqrt(true_V(1,1)^2 + true_V(1,2)^2), sigm); 
        phx = normrnd(atan(true_V(1,2)/true_V(1,1)), sigm);
        
        Ay  = normrnd(sqrt(true_V(2,1)^2 + true_V(2,2)^2), sigm); 
        phy = normrnd(atan(true_V(2,2)/true_V(2,1)), sigm); 
        
        Az  = normrnd(sqrt(true_V(3,1)^2 + true_V(3,2)^2), sigm); 
        phz = normrnd(atan(true_V(3,2)/true_V(3,1)), sigm);
        
        aa(1) = Ax*cos(phx);
        bb(1) = Ax*sin(phx);
        aa(2) = Ay*cos(phy);
        bb(2) = Ay*sin(phy);
        aa(3) = Az*cos(phz);
        bb(3) = Az*sin(phz);
end
%}

%}