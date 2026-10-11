%Симуляция шумов с МНК
N =100;
f = 1e3; %частота
fs = 1e6; %частота дискретизации
T = 1/f;
sampl_size = round(fs/f);

x = linspace(0,10,N);
y = linspace(0,10,N);
z = linspace(0,10,N);
gam = repelem( deg2rad(10), N);
bet = repelem( deg2rad(10), N);
alf = repelem( deg2rad(10), N);
Bm = 2*1e-7*700*pi*(23e-3/2)^2*0.5;%700 витков, диаметр 23 мм, ток 0.5А;

sig = 1e-9;

true_V = zeros(3,2,N);
for i = 4:N
    true_V(:,:,i) = GenV(Bm, alf(i), bet(i), gam(i), x(i), y(i), z(i));
    vx = true_V(1,1,i).*cos(2.*pi.*f.*linspace(0,T,sampl_size)) + true_V(1,2,i).*sin(2.*pi.*f.*linspace(0,T,sampl_size));
    plot(normrnd(vx));
end