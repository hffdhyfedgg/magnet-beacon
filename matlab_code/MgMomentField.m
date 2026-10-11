function [bx,by,bz] = MgMomentField(x,y,z,M,n,dim)

arguments
    x (:,:) double;
    y (:,:) double;
    z (:,:) double;
    n (1,3) double;
    M (1,1) double;
    dim (1,1) int8 = 0;
end

B = @(x,y,z) (1e-7) * (3*([x,y,z])*dot(M*n,[x,y,z])/norm([x,y,z])^5-M*n/norm([x,y,z])^3);

if (size(x,1) ~= size(y,1)) || (size(y,1) ~= size(z,1)) || (size(x,1) ~= size(z,1) || size(x,2) ~= size(y,2)) || (size(y,2) ~= size(z,2)) || (size(x,2) ~= size(z,2))
    error("ћассивы x,y,z должны иметь равные размеры.");
end

if norm(n) ~= 1
    error("Ќорма вектора нормали n должна равн€тс€ 1.");
end    
if size(x,1) == 1 || size(x,2) == 1
    if dim ~= 0
        error("јргумент ""dim"" может быть опущен или должен равн€тс€ нулю, если массивы координат одномерны.");
    end
    
    bx = zeros(size(x));
    by = zeros(size(x));
    bz = zeros(size(x));
    var = zeros(1,3);
    for i = 1:max(size(x))
        var = B(x(i),y(i),z(i));
        bx(i) = var(1);
        by(i) = var(2);
        bz(i) = var(3);
    end
    
else
    if dim > 2 || dim < 1
        error("јргумент ""dim"" должен равн€тс€ 1 или 2. ¬ случае если массивы координат одномерны, dim может быть опущено.");
    end
    
end



end

