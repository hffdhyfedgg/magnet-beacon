clear
load D:\Bezmen\PavelSheptuhin\GenV.mat
Nr =20;
Nang  =100;
sampl_size = 5e3;
sig = 1e-9/22.42/10;
l = linspace(1, 10, Nr);
direction_angle1 = pi/1000;%От -пи до +пи
direction_angle2 = pi/1000;%От -пи/2 до +пи/2
n = [cos(direction_angle2)*cos(direction_angle1), cos(direction_angle2)*sin(direction_angle1), sin(direction_angle2)];
clear direction_angle1 direction_angle2
r = l.*n';
alf =linspace(-pi, pi, Nang);
bet = zeros(1, Nang);
gam =zeros(1, Nang);
bm = 2*1e-7*700*pi*(23e-3/2)^2*0.5;%700 витков, диаметр 23 мм, ток 0.5А
graphicData.alf = zeros(Nr, Nang);
graphicData.bet = zeros(Nr, Nang);
graphicData.gam = zeros(Nr, Nang);

parfor i = 1:Nr
    erAlf = zeros(1, sampl_size);
    erBet = zeros(1, sampl_size);
    erGam = zeros(1, sampl_size);
    for j = 1:Nang
        for k = 1:sampl_size
            a = zeros(1, 3);
            b = zeros(1, 3);
            v = GenV(bm, alf(j), bet(j), gam(j), r(2, i), r(1, i),  r(3, i));

            Ax  = normrnd(sqrt(v(1, 1)^2 + v(1, 2)^2), sig); 
            phx = normrnd(atan(v(1, 2)/v(1, 1)), sig);

            Ay  = normrnd(sqrt(v(2, 1)^2 + v(2, 2)^2), sig); 
            phy = normrnd(atan(v(2, 2)/v(2, 1)), sig); 

            Az  = normrnd(sqrt(v(3, 1)^2 + v(3, 2)^2), sig); 
            phz = normrnd(atan(v(3, 2)/v(3, 1)), sig);

            a(:) = [Ax*cos(phx), Ay*cos(phy), Az*cos(phz)];
            b(:) = [Ax*sin(phx), Ay*sin(phy), Az*sin(phz)];
            
            rtm = real(CalcRotMatrix(bm, a(1), a(2), a(3), b(1), b(2), b(3)));
            [erAlf(k), erBet(k), erGam(k)] = CalcAngles(rtm);

        end
        graphicData(i, j).alf = mean(100*abs((erAlf-alf(j))/alf(j)));
        graphicData(i, j).bet = mean(100*abs((erBet-bet(j))/gam(j)));
        graphicData(i, j).gam = mean(100*abs((erGam-bet(j))/gam(j)));

    end
    disp(i);
end
%% Графики
txtTag = "alfa change (pi)";
zAlf = zeros(Nang, Nr);
zBet = zeros( Nang, Nr);
zGam = zeros( Nang, Nr);

for i = 1:Nang
    for j = 1:Nr
        zAlf(i,j) = graphicData(j, i).alf;
        zBet(i,j) = graphicData(j, i).bet;
        zGam(i,j) = graphicData(j, i).gam;
    end
end
fg = figure("Name", txtTag);
%tab1 = uitab('Title', 'Погрешность радиуса');
parent = fg;
ax1 = axes(parent, 'Position',[0.1 0.1 0.8 0.8]);%, 'title', txtTag,'xlabel', "r [м]",'ylabel', "angle [rad]",'zlabel', "|\delta|r [%]"
ax2 = axes(parent, 'Position',[0.75 0.75 0.20 0.20]);
title("Направление")
[x,y] = meshgrid(l, alf);
plot3(ax1, x, y, zAlf);
plot3(ax2, [0, n(1)], [0, n(2)], [0, n(3)]);
title("Направление")
ax1.XLabel.String =  "r [м]";
ax1.YLabel.String =   "angle [rad]";
ax1.ZLabel.String =  "|\delta|\alpha [%]";
ax1.ZLim = [0 100];
ax1.YTick = [-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi];
ax1.YTickLabel = {'-\pi', '-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2', '\pi'};
saveas(fg, "r_"+txtTag+".png");

fg = figure("Name", txtTag);
%tab2 = uitab('Title', 'Погрешность X');
parent = fg;
ax1 = axes(parent, 'Position',[0.1 0.1 0.8 0.8]);
title(txtTag)
xlabel("r [м]"); ylabel("angle [rad]"); zlabel("|\delta|x [%]");
ax2 = axes(parent, 'Position',[0.75 0.75 0.20 0.20]);
title("Направление")
plot3(ax1, x, y, zBet);
plot3(ax2, [0, n(1)], [0, n(2)], [0, n(3)]);
title("Направление")
ax1.XLabel.String =  "r [м]";
ax1.YLabel.String =   "angle [rad]";
ax1.ZLabel.String =  "|\delta|\beta [%]";
ax1.ZLim = [0 200];
ax1.YTick = [-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi];
ax1.YTickLabel = {'-\pi', '-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2', '\pi'};
saveas(fg, "x_"+txtTag+".png");

fg = figure("Name", txtTag);
%tab3 = uitab('Title', 'Погрешность Y');
parent = fg;
ax1 = axes(parent, 'Position',[0.1 0.1 0.8 0.8]);
title(txtTag)
ax2 = axes(parent, 'Position',[0.75 0.75 0.20 0.20]);
plot3(ax1, x, y, zGam);
plot3(ax2, [0, n(1)], [0, n(2)], [0, n(3)]);
title("Направление")
ax1.XLabel.String =  "r [м]";
ax1.YLabel.String =   "angle [rad]";
ax1.ZLabel.String =  "|\delta|\gamma [%]";
ax1.ZLim = [0 200];
ax1.YTick = [-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi];
ax1.YTickLabel = {'-\pi', '-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2', '\pi'};
saveas(fg, "y_"+txtTag+".png");