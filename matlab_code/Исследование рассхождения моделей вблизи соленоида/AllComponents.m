m = 1;
L = 0.25;
ro = m/L;

beacon = MagneticMoment(0,0,0,m,0,0);

N = 100;
x = linspace(0.01, 5, N);
y = linspace(0.01, 5, N);
z = linspace(0.01, 5, N);

B_coil = zeros(3, numel(x), numel(y), numel(z));%Катушка
B_m = zeros(3, numel(x), numel(y), numel(z));%Магнитный момент

parfor i = 1:numel(x)
    for j = 1:numel(y)
        for k = 1:numel(z)
            B_coil(:, i, j, k) = beacon.Field(x(i), y(i), z(i));
            
        end
    end
end