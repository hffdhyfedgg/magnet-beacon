
n1 = 25;
n2 = 1e6;

sig = 1e-2;



ssig = [1e-3, 1e-2, 1e-1, 10, 100, 1000];


a_alf = struct(...
    'median', cell(1, length(ssig)),...
    'std', cell(1, length(ssig)),...
    'mean', cell(1, length(ssig)),...
    'max', cell(1, length(ssig)),...
    'min', cell(1, length(ssig))...
    );

b_bet = struct(...
    'median', cell(1, length(ssig)),...
    'std', cell(1, length(ssig)),...
    'mean', cell(1, length(ssig)),...
    'max', cell(1, length(ssig)),...
    'min', cell(1, length(ssig))...
    );

g_gam = struct(...
    'median', cell(1, length(ssig)),...
    'std', cell(1, length(ssig)),...
    'mean', cell(1, length(ssig)),...
    'max', cell(1, length(ssig)),...
    'min', cell(1, length(ssig))...
    );


parfor ii = 1:length(ssig)
    sig = ssig(ii);
x = ones(1, n1);
y = ones(1, n1);
z = 2*ones(1, n1);

alf = zeros(1, n1);
bet = deg2rad(linspace(-180, 180, n1));
gam = zeros(1, n1);

v = repmat(GenV_o(1, alf, bet, gam, x, y, z),1,1,1,n2);
v = normrnd(v, sig);
sz = size(v);

a = CalcAll(1, v);
%%
angles = rotm2eul(a.rotmat(:,:,:), 'ZYX');
% gam_r = rad2deg(mean(  abs( reshape(angles(:,1), n1, n2) -repmat(gam', 1, n2))  , 2));%repmat(gam', 1, n2)
% bet_r = rad2deg(mean(  abs( reshape(angles(:,2), n1, n2) - repmat(bet(1:end)', 1, n2))  , 2));
% alf_r = rad2deg(mean(  abs( reshape(angles(:,3), n1, n2) - repmat(alf(1:end)', 1, n2))  , 2));

gam_r = rad2deg(reshape(angles(:,1), n1, n2));
bet_r = rad2deg(reshape(angles(:,2), n1, n2));
alf_r = rad2deg(reshape(angles(:,3), n1, n2));


%%
a_alf(ii) = FillFields(alf_r);
b_bet(ii) = FillFields(bet_r);
g_gam(ii) = FillFields(gam_r);

% boxplot(ax(1), rad2deg(gam), bet_r)
% 
% 
% boxplot(ax(2), rad2deg(gam), alf_r)
% 
% 
% boxplot(ax(3), rad2deg(gam), gam_r)


% plot(ax(1), rad2deg(gam), bet_r, colors(ii))
% 
% 
% plot(ax(2), rad2deg(gam), alf_r, colors(ii))
% 
% 
% plot(ax(3), rad2deg(gam), gam_r, colors(ii))

end
%%
tit = 'std';

figure;
ax(1) = axes();
hold on
ylabel('\alpha')
xlabel('\beta')
title(tit)

figure;
ax(2) = axes();
hold on
ylabel('\beta')
xlabel('\beta')
title(tit)

figure;
ax(3) = axes();
hold on
ylabel('\gamma')
xlabel('\beta')
title(tit)

colors = ['b', 'r', 'g', 'm', 'c', 'k'];
for i = 1:length(ssig)
   plot(ax(1), linspace(-180, 180, n1), a_alf(i).std, colors(i))
   plot(ax(2), linspace(-180, 180, n1), b_bet(i).std, colors(i))
   plot(ax(3), linspace(-180, 180, n1), g_gam(i).std, colors(i))
end


function res_struct = FillFields(var)
    res_struct.median = median(var, 2);
    res_struct.std = std(var, [], 2);
    res_struct.mean = mean(var, 2);
    res_struct.max = max(var, [], 2);
    res_struct.min = min(var, [], 2);
end