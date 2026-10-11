load('GenV.mat')
bm = 1e-6;

n = 75;
target = 0.1;
l = (linspace(0.1, 3, n));
rr = zeros(1, n);
parfor i = 1:numel(l)

    GenVS = CreateABParameterGenerator({'l', 'cos', l(i), [1, 0, 0]; 'l', 'sin', l(i), [0, 1, 0]});

    f = @(t) CalcErr(bm, l(i), t, repmat(1e-3, size(t)), repmat(1e-3, size(t)), zeros(size(t)), zeros(size(t)), zeros(size(t)), GenV, GenVS, [14, 14]);
    %{1 rtrue, 2 xtrue, 3 ytrue, 4 ztrue, 5 rsol, 6 xsol, 7 ysol, 8 zsol, 9 rer, 10 xer, 11 yer, 12 zer, 13 rreler, 14 xreler, 15 yreler, 16 zreler}

    [r, v] = Find(f, target, 2e-2, l(i)*4);

    rr(i) = r;

end

%%
figure;
plot(l, rr)

%%
GenVS = CreateABParameterGenerator({'l', 'cos', l(i), [1, 0, 0]; 'l', 'sin', l(i), [0, 1, 0]});