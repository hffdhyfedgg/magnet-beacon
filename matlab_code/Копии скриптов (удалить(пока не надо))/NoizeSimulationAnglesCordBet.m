clear
load D:\Bezmen\PavelSheptuhin\GenV.mat
Nr =20;
Nang  =100;
sampl_size = 1e4;
sig = 1e-9/22.42;
l = linspace(0.1, 10, Nr);
direction_angle1 = pi/4;%От -пи до +пи
direction_angle2 = pi/4;%От -пи/2 до +пи/2
n = [1/sqrt(3),1/sqrt(3),1/sqrt(3),];%n = [cos(direction_angle2)*cos(direction_angle1), cos(direction_angle2)*sin(direction_angle1), sin(direction_angle2)];
clear direction_angle1 direction_angle2
r = l.*n';
alf = zeros(1, Nang);
bet = linspace(-pi/2, pi/2, Nang);
gam = zeros(1, Nang);
bm = 2*1e-7*700*pi*(23e-3/2)^2*0.5;%700 витков, диаметр 23 мм, ток 0.5А
graphicData.r = zeros(Nr, Nang);
graphicData.x = zeros(Nr, Nang);
graphicData.y = zeros(Nr, Nang);
graphicData.z = zeros(Nr, Nang);


parfor i = 1:Nr
    
    erR = zeros(1, sampl_size);
    erX = zeros(1, sampl_size);
    erY = zeros(1, sampl_size);
    erZ = zeros(1, sampl_size);
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

            erR(k) = CalcR(bm, a(1), a(2), a(3), b(1), b(2), b(3));
            erX(k) = CalcX(bm, a(1), a(2), a(3), b(1), b(2), b(3));
            erY(k) = CalcY(bm, a(1), a(2), a(3), b(1), b(2), b(3));
            erZ(k) = CalcZ(bm, a(1), a(2), a(3), b(1), b(2), b(3));
        end
        graphicData(i, j).r = mean(100*abs((erR-l(i))/l(i)));
        graphicData(i, j).x = mean(100*abs((erX-r(1, i))/r(1, i)));
        graphicData(i, j).y = mean(100*abs((erY-r(2, i))/r(2, i)));
        graphicData(i, j).z = mean(100*abs((erZ-r(3, i))/r(3, i)));
    end
    disp(i);
end
%% Графики
txtTag = "AngCordBeta_mil";
r = zeros(Nang, Nr);
zdx = zeros( Nang, Nr);
zdy = zeros( Nang, Nr);
zdz = zeros( Nang, Nr);
for i = 1:Nang
    for j = 1:Nr
        r(i,j) = graphicData(j, i).r;
        zdx(i,j) = graphicData(j, i).x;
        zdy(i,j) = graphicData(j, i).y;
        zdz(i,j) = graphicData(j, i).z;
    end
end
fg = figure("Name", txtTag);
%tab1 = uitab('Title', 'Погрешность радиуса');
parent = fg;
ax1 = axes(parent, 'Position',[0.1 0.1 0.8 0.8]);%, 'title', txtTag,'xlabel', "r [м]",'ylabel', "angle [rad]",'zlabel', "|\delta|r [%]"
ax2 = axes(parent, 'Position',[0.75 0.75 0.20 0.20]);
title("Направление")
[x,y] = meshgrid(l, bet);
plot3(ax1, x, y, r);
plot3(ax2, [0, n(1)], [0, n(2)], [0, n(3)]);
title("Направление")
ax1.XLabel.String =  "r [м]";
ax1.YLabel.String =   "angle [rad]";
ax1.ZLabel.String =  "|\delta|r [%]";
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
plot3(ax1, x, y, zdx);
plot3(ax2, [0, n(1)], [0, n(2)], [0, n(3)]);
title("Направление")
ax1.XLabel.String =  "r [м]";
ax1.YLabel.String =   "angle [rad]";
ax1.ZLabel.String =  "|\delta|x [%]";
ax1.ZLim = [0 200];
ax1.YTick = [-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi];
ax1.YTickLabel = {'-\pi', '-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2', '\pi'};
saveas(fg, "x_"+txtTag+".png");

fg = figure("Name", txtTag);
%tab3 = uitab('Title', 'Погрешность Y');
parent = fg;
ax1 = axes(parent, 'Position',[0.1 0.1 0.8 0.8]);
title(txtTag)
xlabel("r [м]"); ylabel("angle [rad]"); zlabel("|\delta|y [%]");
ax2 = axes(parent, 'Position',[0.75 0.75 0.20 0.20]);
plot3(ax1, x, y, zdy);
plot3(ax2, [0, n(1)], [0, n(2)], [0, n(3)]);
title("Направление")
ax1.XLabel.String =  "r [м]";
ax1.YLabel.String =   "angle [rad]";
ax1.ZLabel.String =  "|\delta|y [%]";
ax1.ZLim = [0 200];
ax1.YTick = [-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi];
ax1.YTickLabel = {'-\pi', '-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2', '\pi'};
saveas(fg, "y_"+txtTag+".png");

fg = figure("Name", txtTag);
%tab4 = uitab('Title', 'Погрешность Z');
parent = fg;
ax1 = axes(parent, 'Position',[0.1 0.1 0.8 0.8]);
title(txtTag)
xlabel("r [м]"); ylabel("angle [rad]"); zlabel("|\delta|z [%]");
ax2 = axes(parent, 'Position',[0.75 0.75 0.20 0.20]);
plot3(ax1, x, y, zdz);
plot3(ax2, [0, n(1)], [0, n(2)], [0, n(3)]);
title("Направление")
ax1.XLabel.String =  "r [м]";
ax1.YLabel.String =   "angle [rad]";
ax1.ZLabel.String =  "|\delta|z [%]";
ax1.ZLim = [0 200];
ax1.YTick = [-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi];
ax1.YTickLabel = {'-\pi', '-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2', '\pi'};
saveas(fg, "z_"+txtTag+".png");