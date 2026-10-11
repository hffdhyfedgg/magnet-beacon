function v = GenV_o(bm, alf, bet, gam, x, y, z)
arguments
    bm  (:,1) {mustBeNumeric, mustBeReal};
    alf (:,1) {mustBeNumeric, mustBeReal};
    bet (:,1) {mustBeNumeric, mustBeReal};
    gam (:,1) {mustBeNumeric, mustBeReal};
    x (:,1) {mustBeNumeric, mustBeReal};
    y (:,1) {mustBeNumeric, mustBeReal};
    z (:,1) {mustBeNumeric, mustBeReal};
end
v = genV(bm(:), alf(:), bet(:), gam(:), x(:), y(:), z(:), size(alf));
end


function v = genV(bm, alf, bet, gam, x, y, z, sz)

if all(size(bm) ~= [1,1])
    bm = repmat(bm, 1, 6);
end

var1 = 1.0 ./ (abs(y).^2 + x.^2 + z.^2).^(2.5).*3.0;
var2 = ((y.^2.*3.0) ./ (abs(y).^2 + x.^2 + z.^2) - 1.0).*1.0 ./ (abs(y).^2 + x.^2 + z.^2).^(1.5);
r =  (abs(y).^2 + x.^2 + z.^2);
x2 = x.^2 .* 3.0;
xy = x .* y;
xz = x .* z;
yz = y .* z;


c_a = cos(alf);
c_b = cos(bet);
c_g = cos(gam);
s_a = sin(alf);
s_b = sin(bet);
s_g = sin(gam);

var3 = s_a .* s_g + c_a .* c_g .* s_b;
var4 = c_a .* s_g - c_g .* s_a .* s_b;
var5 = c_a .* c_g + s_a .* s_b .* s_g;
var6 = c_g .* s_a - c_a .* s_b .* s_g;



% row_v = bm*[...
%   - var4 .* var2 + yz .* var3 .* var1 + xy .* c_b .* c_g .* var1,...
%     var5 .* var2 - yz .* var6 .* var1 + xy .* c_b .* s_g .* var1,...
%     xy .* s_b .* 1.0 ./ r.^(2.5) .*  - 3.0 + c_b .* s_a .* var2 + yz .* c_a .* c_b .* var1,...
%     c_b .* c_g .* (x2 ./ r - 1.0) .* 1.0 ./ r.^(1.5) - xy .* var4 .* var1 + xz .* var3 .* var1,...
%     c_b .* s_g .* (x2 ./ r - 1.0) .* 1.0 ./ r.^(1.5) + xy .* var5 .* var1 - xz .* var6 .* var1,...
%   - s_b .* (x2 ./ r - 1.0) .* 1.0 ./ r.^(1.5) + xz .* c_a .* c_b .* var1 + xy .* c_b .* s_a .* var1];

% disp(row_v);
% v = row_v;
v = reshape( bm*[...
  - var4 .* var2 + yz .* var3 .* var1 + xy .* c_b .* c_g .* var1,...
    var5 .* var2 - yz .* var6 .* var1 + xy .* c_b .* s_g .* var1,...
    xy .* s_b .* 1.0 ./ r.^(2.5) .*  - 3.0 + c_b .* s_a .* var2 + yz .* c_a .* c_b .* var1,...
    c_b .* c_g .* (x2 ./ r - 1.0) .* 1.0 ./ r.^(1.5) - xy .* var4 .* var1 + xz .* var3 .* var1,...
    c_b .* s_g .* (x2 ./ r - 1.0) .* 1.0 ./ r.^(1.5) + xy .* var5 .* var1 - xz .* var6 .* var1,...
  - s_b .* (x2 ./ r - 1.0) .* 1.0 ./ r.^(1.5) + xz .* c_a .* c_b .* var1 + xy .* c_b .* s_a .* var1]',...
[3,2,sz]);

% @(Bm,alf,bet,gam,x,y,z)reshape(
% [-Bm.*(cos(alf).*sin(gam)-cos(gam).*sin(alf).*sin(bet)).*((y.^2.*3.0)./(abs(y).^2+x.^2+z.^2)-1.0).*1.0./(abs(y).^2+x.^2+z.^2).^(3.0./2.0)+Bm.*y.*z.*(sin(alf).*sin(gam)+cos(alf).*cos(gam).*sin(bet)).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0+Bm.*x.*y.*cos(bet).*cos(gam).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0,
%     Bm.*(cos(alf).*cos(gam)+sin(alf).*sin(bet).*sin(gam)).*((y.^2.*3.0)./(abs(y).^2+x.^2+z.^2)-1.0).*1.0./(abs(y).^2+x.^2+z.^2).^(3.0./2.0)-Bm.*y.*z.*(cos(gam).*sin(alf)-cos(alf).*sin(bet).*sin(gam)).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0+Bm.*x.*y.*cos(bet).*sin(gam).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0,
%     Bm.*x.*y.*sin(bet).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*-3.0+Bm.*cos(bet).*sin(alf).*((y.^2.*3.0)./(abs(y).^2+x.^2+z.^2)-1.0).*1.0./(abs(y).^2+x.^2+z.^2).^(3.0./2.0)+Bm.*y.*z.*cos(alf).*cos(bet).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0,
%     Bm.*cos(bet).*cos(gam).*((x.^2.*3.0)./(abs(y).^2+x.^2+z.^2)-1.0).*1.0./(abs(y).^2+x.^2+z.^2).^(3.0./2.0)-Bm.*x.*y.*(cos(alf).*sin(gam)-cos(gam).*sin(alf).*sin(bet)).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0+Bm.*x.*z.*(sin(alf).*sin(gam)+cos(alf).*cos(gam).*sin(bet)).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0,
%     Bm.*cos(bet).*sin(gam).*((x.^2.*3.0)./(abs(y).^2+x.^2+z.^2)-1.0).*1.0./(abs(y).^2+x.^2+z.^2).^(3.0./2.0)+Bm.*x.*y.*(cos(alf).*cos(gam)+sin(alf).*sin(bet).*sin(gam)).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0-Bm.*x.*z.*(cos(gam).*sin(alf)-cos(alf).*sin(bet).*sin(gam)).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0,
%     -Bm.*sin(bet).*((x.^2.*3.0)./(abs(y).^2+x.^2+z.^2)-1.0).*1.0./(abs(y).^2+x.^2+z.^2).^(3.0./2.0)+Bm.*x.*z.*cos(alf).*cos(bet).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0+Bm.*x.*y.*cos(bet).*sin(alf).*1.0./(abs(y).^2+x.^2+z.^2).^(5.0./2.0).*3.0],
% [3,2])

end
