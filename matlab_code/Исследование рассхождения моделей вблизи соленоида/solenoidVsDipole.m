N = 100;
x = linspace( 0.1, 4, N);
y = linspace( 0.1, 4, N);
z = linspace( 0.1, 4, N);
r = sqrt(x.^2 + y.^2 + z.^2);
gam = repelem( deg2rad(0), N);
bet = repelem( deg2rad(0), N);
alf = repelem( deg2rad(0), N);

Bm = 1e6*2*1e-7*700*pi*(23e-3/2)^2*0.5;%700 витков, диаметр 23 мм, ток 0.5А;
V = zeros(3, 2);
true_V = zeros(3, 2);
aa = zeros(3, N);
bb = zeros(3, N);
erV = zeros(1, N);
ervv = zeros(3, 2, N);

errR = zeros(1,N);   abserrR = zeros(1,N);
errX = zeros(1,N);   abserrX = zeros(1,N);
errY = zeros(1,N);   abserrY = zeros(1,N);
errZ = zeros(1,N);   abserrZ = zeros(1,N);
erX = zeros(1, N);
erY = zeros(1, N);
erZ = zeros(1, N);
erR = zeros(1, N);

for i = 1:N
    true_V(:,:) = GenV(Bm, alf(i), bet(i), gam(i), y(i), x(i), z(i));
    V(:,:) =  GenVS(Bm/0.25, alf(i), bet(i), gam(i), y(i), x(i), z(i));
    aa(1,i) = V(1,1);   aa(2, i) = V(2,1) ;   aa(3,i) = V(3,1);
    bb(1,i) = V(1,2);   bb(2,i) = V(2,2);   bb(3,i) = V(3,2);
    erV(i) = 100*abs(sqrt(true_V(1)^2+true_V(2)^2+true_V(3)^2+true_V(4)^2+true_V(5)^2+true_V(6)^2) - sqrt(V(1)^2+V(2)^2+V(3)^2+V(4)^2+V(5)^2+V(6)^2))/Bm;
    ervv(:, :, i) = true_V - V;
    
    erR(i) = CalcR(Bm, V(1,1), V(2,1), V(3,1), V(1,2), V(2,2), V(3,2));
    erX(i) = CalcX(Bm, V(1,1), V(2,1), V(3,1), V(1,2), V(2,2), V(3,2));
    erY(i) = CalcY(Bm, V(1,1), V(2,1), V(3,1), V(1,2), V(2,2), V(3,2));
    erZ(i) = CalcZ(Bm, V(1,1), V(2,1), V(3,1), V(1,2), V(2,2), V(3,2));

    abserrR(i) = abs( r(i) - erR(i) );
    abserrX(i) = abs( x(i) - erX(i) );
    abserrY(i) = abs( y(i) - erY(i) );
    abserrZ(i) = abs( z(i) - erZ(i) );

    errR(i) = 100*abserrR(i)./abs(r(i));
    errX(i) = 100*abserrX(i)./abs(x(i));
    errY(i) = 100*abserrY(i)./abs(y(i));
    errZ(i) = 100*abserrZ(i)./abs(z(i));
    
end


%%
%{
xstr = "r [м]";
ystr = "|\delta| ? [%]";
strtitle = "Bm = "+num2str(Bm);
figure
semilogy(r,  errX)
xlabel(xstr); ylabel(strrep(ystr, '?', 'x')); title(strtitle);

figure
semilogy(r,  errY)
xlabel(xstr); ylabel(strrep(ystr, '?', 'y')); title(strtitle);

figure
semilogy(r,  errZ)
xlabel(xstr); ylabel(strrep(ystr, '?', 'z')); title(strtitle);

figure
semilogy(r,  errR)
xlabel(xstr); ylabel(strrep(ystr, '?', 'r')); title(strtitle);

figure
semilogy(r, erV)
xlabel(xstr); ylabel("|(V_s - V)|/B_m"); title(strtitle);
%}
titles = {'ax', 'bx'; 'ay', 'by'; 'az', 'bz'};

for i = 1:3
    for k = 1:2
        ydata = zeros(1,N);
        for j = 1:N
            ydata(j) = ervv(i, k, j);
        end
        figure;
        plot(r, ydata)
        title(titles{i, k})
    end
end