load rr1.mat
load rr2.mat
boot_xdata = 2*1e-7*700*pi*(23e-3/2)^2*0.5./exp(linspace(log(0.5*1e-9/22.42), log(300*1e-9/22.42), 10));

x = boot_xdata;

a = optimvar('a', 2);
n = 3;
fn = a(2)*(x.^(1/n))+a(1)*(x.^1);

minim = sum((fn-rr1).^2); 
pr = optimproblem("Objective", minim);
sl = solve(pr);

minim = sum((fn-rr2).^2); 
pr = optimproblem("Objective", minim);
sl2 = solve(pr);

disp("---------------------------------------------------")
disp("Константы параметра превышения мнимой a1 = " + num2str( sl.a(1) ) + " a2 = " + num2str( sl.a(2) ) );
disp("Константы параметра появления мнимой a1 = " + num2str( sl2.a(1) ) + " a2 = " + num2str( sl2.a(2) ) );

%%
fg_rr1 = figure('Name', 'fg_rr1');
xlabel("Bm/sig"); ylabel("r [м]");
%legend("Рез-ы числ. моделирования", "Рез-ы числ. моделирования","Апроксимаци","Апроксимация");
hold on

plot(boot_xdata, rr1, "xred");
plot(boot_xdata, rr2, "xblue");
xx = linspace( min(boot_xdata), max(boot_xdata), 200);
plot(xx, sl.a(2)*(xx.^(1/n))+sl.a(1)*(xx.^1), "red");
plot(xx, sl2.a(2)*(xx.^(1/n))+sl2.a(1)*(xx.^1), "blue");

delt = sl.a(2)*(boot_xdata.^(1/n))+sl.a(1)*(boot_xdata.^1)-rr1;
disp("Превышение мнимой над реальной")
disp("Разность суммы квадратов: "+num2str(sum(delt.^2)));
disp("Абсолютная максимальная ошибка апроксимации: " + num2str( max( abs(delt) ) ));
disp("Относительная максимальная ошибка апроксимации: " + num2str( max( 100*abs(delt/rr1) ) ) + " %");

delt = sl2.a(2)*(boot_xdata.^(1/n))+sl2.a(1)*(boot_xdata.^1)-rr2;
disp("Появление мнимой части")
disp("Разность суммы квадратов: "+num2str(sum(delt.^2)));
disp("Абсолютная максимальная ошибка апроксимации: " + num2str( max( abs(delt) ) ));
disp("Относительная максимальная ошибка апроксимации: " + num2str( max( 100*abs(delt/rr2) ) ) + " %");

