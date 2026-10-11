function rotationMatrix = CalcRotMatrix_matr(b, ax, ay, az, bx, by, bz)
arguments
    b(3, 1, :)  double {mustBeNumeric, mustBeReal};
    ax(3, 1, :) double {mustBeNumeric, mustBeReal};
    ay(3, 1, :) double {mustBeNumeric, mustBeReal};
    az(3, 1, :) double {mustBeNumeric, mustBeReal};
    bx(3, 1, :) double {mustBeNumeric, mustBeReal};
    by(3, 1, :) double {mustBeNumeric, mustBeReal};
    bz(3, 1, :) double {mustBeNumeric, mustBeReal};
end
    sz = size(ax);
    n = 0;
    for i = 1:size(sz,2)
        
        if(sz(i) > 1)
            n = n+1;
        end
        
        if n > 1
            error("Аргументы должны быть являться одномерным вектором по одному из измерений. (b может быть скаляром)")
        end
    end
    clear sz n
    
    
        
    if ~(numel(b) == numel(ax) && numel(b) == numel(ay) && numel(b) == numel(az) && numel(b) == numel(bx) && numel(b) == numel(by) && numel(b) == numel(bz)) && ~(numel(ax) == numel(ay) && numel(ax) == numel(az) && numel(ax) == numel(bx) && numel(ax) == numel(by) && numel(ax) == numel(bz) && isscalar(b))
        error("Аргументы должны быть одинаковой размерности, или аргумент b должен быть числом, а размерность остальных аргументов должна быть одинаковой.");
    end
    
    if isscalar(b)
        bb = repmat(b, size(ax));
    else        
        bb = b;
    end
   
    r = CalcR(bb, ax, ay, az, bx, by, bz);
    %{
    x = CalcX(bb, ax, ay, az, bx, by, bz);
    y = CalcY(bb, ax, ay, az, bx, by, bz);
    
    z = CalcZ(bb, ax, ay, az, bx, by, bz);
    
    z2 = sqrt(r.^2-x.^2-y.^2);
    %}
    u = CalcX(bb, ax, ay, az, bx, by, bz)./r;
    v = CalcY(bb, ax, ay, az, bx, by, bz)./r;
    
    w = real(sqrt(1- u.^2 - v.^2));
    outvar = zeros(3,3,numel(ax));

        c = cross([ax, ay, az],[bx, by, bz]);
        
        
        cc = real(cross(...
                    ( bb./r.^3 )*[ 3*u(i)*v(i); 3*v(i)^2-1; 3*v(i)*w(i) ], ...
                    (bb(i)/r(i)^3)*[3*u(i)^2-1; 3*u(i)*v(i); 3*u(i)*w(i)] ));
                %{
                A = [[ax(i), ay(i), az(i)]', [bx(i), by(i), bz(i)]',c'];
                B = [ ...
                            (bb(i)/r(i)^3)*[3*u(i)*v(i); 3*v(i)^2-1; 3*v(i)*w(i)], ...
                            (bb(i)/r(i)^3)*[3*u(i)^2-1; 3*u(i)*v(i); 3*u(i)*w(i)], ...
                             cc];
                    %}
%         var = inv([ ...
%                 (bb(i)/r(i)^3)*[3*u(i)*v(i); 3*v(i)^2-1; 3*v(i)*w(i)], ...
%                 (bb(i)/r(i)^3)*[3*u(i)^2-1; 3*u(i)*v(i); 3*u(i)*w(i)], ...
%                  cc]);
        var = [ ...
                (bb(i)/r(i)^3)*[3*u(i)*v(i); 3*v(i)^2-1; 3*v(i)*w(i)], ...
                (bb(i)/r(i)^3)*[3*u(i)^2-1; 3*u(i)*v(i); 3*u(i)*w(i)], ...
                 cc];

        outvar(:,:,i) = [[ax(i), ay(i), az(i)]', [bx(i), by(i), bz(i)]',c']/var;
        %{
        inv([ ...
                            (bb(i)/r(i)^3)*[3*u(i)*v(i); 3*v(i)^2-1; 3*v(i)*w(i)], ...
                            (bb(i)/r(i)^3)*[3*u(i)^2-1; 3*u(i)*v(i); 3*u(i)*w(i)], ...
                             cc]);
        %}
        outvar(1,:,i) = outvar(1,:,i)./norm(outvar(1,:,i));
        outvar(2,:,i) = outvar(2,:,i)./norm(outvar(2,:,i));
        outvar(3,:,i) = outvar(3,:,i)./norm(outvar(3,:,i));
    
    rotationMatrix = outvar;
end

