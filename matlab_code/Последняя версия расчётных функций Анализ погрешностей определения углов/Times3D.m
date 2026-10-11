function m = Times3D(m1,m2)

if size(m1,1) ~= size(m1,2)
    error('Не квадратная матрица.')
end
if ~all(size(m1) == size(m2))
    error('Размеры матриц не совпадают.')
end

if length(size(m1)) ~= 3 || length(size(m2)) ~= 3
    error('Матрицы должны быть трёхмерными.')
end

n = size(m1, 1);
k = size(m1, 3);


aa = repelem(m1, n,1,1);
bb = repmat(permute(m2, [2,1,3]), n,1,1);
pre_cc = aa.*bb;
pre_cc = sum(pre_cc,2);
m = permute(reshape(pre_cc,3,3,k), [2,1,3]);

end

