N = 1000;%[2:10,10:10:300];
samples = 1e4;
f = 1;
T = 1/f;
ph = deg2rad(25);
ampl = 1;
sigm = linspace(ampl, 100*ampl, 10);
%asinc1 = parfeval(@myStd, 2, 5e4);
%[ans1, ans2] = fetchOutputs(asinc1)
%[ans1, ans2] = myStd()
%% Вычисления погрешности a и b параметров
res_A = zeros(samples,size(N, 2));
res_B = zeros(samples,size(N, 2));
m_A = zeros(1,size(N, 2));
m_B = zeros(1,size(N, 2));
resAmpl = zeros(samples,size(N, 2));
resPhi = zeros(samples,size(N, 2));
stdAmpl = zeros(1,N);
stdPhi = zeros(1,N);

std_A = zeros(size(sigm));
std_B = zeros(size(sigm));

%{
fg = figure();
ax1 = axes(fg);
ax1.Title.String = "Дисперсия a и b параметров.";
pl = plot(ax1, [0 0], [0 0], [0 0], [0 0]);
legend("std A", "std B"); xlabel("n"); ylabel("a,b");
pl(1).XData = []; pl(1).YData = [];
pl(2).XData = []; pl(2).YData = [];

fg2 = figure();
ax21 = subplot(2,1,1); ax21.Title.String = "Дисперсия амплитуды."; ax21.XLabel.String = "n"; ax21.YLabel.String = "A [Тл]";
ax22 = subplot(2,1,2); ax22.Title.String = "Дисперсия фазы."; ax22.XLabel.String = "n"; ax22.YLabel.String = "Фаза [град.]";

pl_ampl = plot(ax21, [0 0], [0 0]);
pl_ampl.XData = []; pl_ampl.YData = [];
pl_ph = plot(ax22, [0 0], [0 0]);
pl_ph.XData = []; pl_ph.YData = [];
%}


for j = 1:size(sigm,2)
for i = 1:size(N,2)
    if i ~= 0
        disp(i);
    end
    sig = sigm(j);
    parfor k = 1:samples
        %clear A B data fn minim pr 
        data = normrnd(ampl*cos(2*pi*f.*linspace(0, T, N(i))+ph), sig);
        
        A = optimvar('A', 1);
        B = optimvar('B', 1);
        fn = A*cos(2*pi*f*linspace(0, T, N(i)))+B*sin(2*pi*f*linspace(0, T, N(i)));
        minim = sum((fn-data).^2);
        pr = optimproblem("Objective", minim);
        clc;
        disp(k);
        %show(pr)
        sl = solve(pr);
        res_A(k,i) = sl.A;
        res_B(k,i) = sl.B;
        resAmpl(k, i) = sqrt( res_A(k, i)^2 + res_B(k, i)^2 );
        resPhi(k, i) = atan( res_B(k, i)/res_A(k, i) );
        
        %Графики дисперсии амплитуды и фазы
        %pl_ampl.XData = [pl_ampl.XData,k]; pl_ampl.YData = [pl_ampl.YData, sig/std(resAmpl(1:k,i))];
        %pl_ph.XData = [pl_ph.XData,k]; pl_ph.YData = [pl_ph.YData, sig/std(resPhi(1:k,i))];
        %Графики дисперсии a и b параметров
        %pl(1).XData = [pl(1).XData,k]; pl(1).YData = [pl(1).YData, std(res_A(1:k,i))];
        %pl(2).XData = [pl(2).XData,k]; pl(2).YData = [pl(2).YData, std(res_B(1:k,i))];
        %plot(linspace(0, T, N(i)),data,"blue",linspace(0, T, N(i)),sl.A*cos(2*pi*f.*linspace(0, T, N(i))) + sl.B*sin(2*pi*f.*linspace(0, T, N(i))),"red")
        
        %stdAmpl(k) = sig/std(resAmpl(1:k,i));
        %stdPhi(k) = sig/std(resPhi(1:k,i));
    end
    m_A(i) = mean(res_A(:,i));
    m_B(i) = mean(res_B(:,i));
    %plot(linspace(0, T, N(i)),data,"blue",linspace(0, T, 100),m_A(i)*cos(2*pi*f.*linspace(0, T, 100)) + m_B(i)*sin(2*pi*f.*linspace(0, T, 100)),"red")

end
std_A(j) = std(res_A);
std_B(j) = std(res_B);

end



%mean(pl(1).YData)
%mean(pl(2).YData)
%% Графики
fg1 = figure;
semilogy(sigm,std_A);
title("СКО А")


fg11 = figure;
semilogy(sigm,sigm./std_A)
title("Отношение изначальной дисперсии к СКО А")

fg2 = figure;
semilogy(sigm,std_B);
title("СКО В")


fg21 = figure;
semilogy(sigm,sigm./std_B)
title("Отношение изначальной дисперсии к СКО B")
%%
fg3 = figure;
hold on
title("Все графики")
for  i = 3:size(N,2)
plot3(repelem(N(i),N(i)), linspace(0, T, N(i)), normrnd(ampl*cos(2*pi*f.*linspace(0, T, N(i))+ph), sig), "blue");

plot3(repelem(N(i),100), linspace(0, T, 100), m_A(i).*cos(2*pi*f.*linspace(0, T, 100))+m_B(i).*sin(2*pi*f.*linspace(0, T,100)), "red", "LineWidth", 2);
end



function [stdA, stdPh] = myStd(samples)
N = 1000;%[2:10,10:10:300];
%samples = 100;
f = 1;
T = 1/f;
ph = deg2rad(25);
ampl = 1;
sig = 1e-1;

res_A = zeros(samples,size(N, 2));
res_B = zeros(samples,size(N, 2));
m_A = zeros(1,size(N, 2));
m_B = zeros(1,size(N, 2));
resAmpl = zeros(samples,size(N, 2));
resPhi = zeros(samples,size(N, 2));
stdAmpl = zeros(1,samples);
stdPhi = zeros(1,samples);
for i = 1:size(N,2)
    if i ~= 0
        disp(i);
    end
    
    for k = 1:samples
        clear A B data fn minim pr 
        data = normrnd(ampl*cos(2*pi*f.*linspace(0, T, N(i))+ph), sig);
        
        A = optimvar('A', 1);
        B = optimvar('B', 1);
        fn = A*cos(2*pi*f*linspace(0, T, N(i)))+B*sin(2*pi*f*linspace(0, T, N(i)));
        minim = sum((fn-data).^2);
        pr = optimproblem("Objective", minim);
        clc;
        disp(k);
        %show(pr)
        sl = solve(pr);
        res_A(k,i) = sl.A;
        res_B(k,i) = sl.B;
        resAmpl(k, i) = sqrt( res_A(k, i)^2 + res_B(k, i)^2 );
        resPhi(k, i) = atan( res_B(k, i)/res_A(k, i) );
        
        %Графики дисперсии амплитуды и фазы
        %pl_ampl.XData = [pl_ampl.XData,k]; pl_ampl.YData = [pl_ampl.YData, sig/std(resAmpl(1:k,i))];
        %pl_ph.XData = [pl_ph.XData,k]; pl_ph.YData = [pl_ph.YData, sig/std(resPhi(1:k,i))];
        %Графики дисперсии a и b параметров
        %pl(1).XData = [pl(1).XData,k]; pl(1).YData = [pl(1).YData, std(res_A(1:k,i))];
        %pl(2).XData = [pl(2).XData,k]; pl(2).YData = [pl(2).YData, std(res_B(1:k,i))];
        %plot(linspace(0, T, N(i)),data,"blue",linspace(0, T, N(i)),sl.A*cos(2*pi*f.*linspace(0, T, N(i))) + sl.B*sin(2*pi*f.*linspace(0, T, N(i))),"red")
        
        stdAmpl(k) = sig/std(resAmpl(1:k, i));
        stdPhi(k) = sig/std(resPhi(1:k, i));
    end
    m_A(i) = mean(res_A(:,i));
    m_B(i) = mean(res_B(:,i));
    %plot(linspace(0, T, N(i)),data,"blue",linspace(0, T, 100),m_A(i)*cos(2*pi*f.*linspace(0, T, 100)) + m_B(i)*sin(2*pi*f.*linspace(0, T, 100)),"red")

end
stdA = mean(stdAmpl(round(samples/5):end));
stdPh = mean(stdPhi(round(samples/5):end));
end
