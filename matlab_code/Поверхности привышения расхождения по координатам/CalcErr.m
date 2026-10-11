function results = CalcErr(bm, l, x, y, z, a, b, g, genV, genVS, ind)%a - alpha, b - beta, g - gamma
xtrue = zeros(1, numel(x));
ytrue = zeros(1, numel(x));
ztrue = zeros(1, numel(x));
xsol = zeros(1, numel(x));
ysol = zeros(1, numel(x));
zsol = zeros(1, numel(x));

xer = zeros(1, numel(x));
yer = zeros(1, numel(x));
zer = zeros(1, numel(x));

rtrue = zeros(1, numel(x));
rsol = zeros(1, numel(x));
rer = zeros(1, numel(x));
rreler = zeros(1, numel(x));

for i = 1:numel(x)
    v = genV(bm, a(i), b(i), g(i), y(i), x(i), z(i));
    vs = genVS(0.5*bm/l, a(i), b(i), g(i), y(i), x(i), z(i));
    
    xtrue(i) = CalcX(bm, v(1), v(2), v(3), v(4), v(5), v(6));
    ytrue(i) = CalcY(bm, v(1), v(2), v(3), v(4), v(5), v(6));
    %ztrue(i) = CalcZ(bm, v(1), v(2), v(3), v(4), v(5), v(6));
    rtrue(i) = CalcR(bm, v(1), v(2), v(3), v(4), v(5), v(6));
    
    ztrue(i) = sqrt(1 - (xtrue(i)/rtrue(i))^2 - (ytrue(i)/rtrue(i))^2);
    
    xsol(i) = CalcX(bm, vs(1), vs(2), vs(3), vs(4), vs(5), vs(6));
    ysol(i) = CalcY(bm, vs(1), vs(2), vs(3), vs(4), vs(5), vs(6));
    %zsol(i) = CalcZ(bm, vs(1), vs(2), vs(3), vs(4), vs(5), vs(6));
    rsol(i) = CalcR(bm, vs(1), vs(2), vs(3), vs(4), vs(5), vs(6));
    
    
    zsol(i) = sqrt(1 - (xsol(i)/rsol(i))^2 - (ysol(i)/rsol(i))^2);
    
    xer(i) = xsol(i) - xtrue(i);
    yer(i) = ysol(i) - ytrue(i);
    zer(i) = zsol(i) - ztrue(i);
    rer(i) = rsol(i) - rtrue(i);
    
    xreler(i) = xer(i)/abs(xtrue(i));
    yreler(i) = yer(i)/abs(ytrue(i));
    zreler(i) = zer(i)/abs(ztrue(i));
    rreler(i) = rer(i)/abs(rtrue(i)); 
end
%results = {[xtrue; ytrue; ztrue], [xsol; ysol; zsol], [xer; yer; zer], [xreler; yreler; zreler]};
al = {rtrue, xtrue, ytrue, ztrue, rsol, xsol, ysol, zsol, rer, xer, yer, zer, rreler, xreler, yreler, real(zreler)};
results = al(ind(1):ind(2));
end


