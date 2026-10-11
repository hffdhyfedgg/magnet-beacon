function pos = CalcAll(b, v)

sz = [size(v),1];
vv = v(:,:,:);

pos = calcA(b, vv, sz(3:end));


end

function pos = calcA(b, v_field, sz)
    
    asq = v_field(1,1,:).^2+v_field(2,1,:).^2+v_field(3,1,:).^2;
    bsq = v_field(1,2,:).^2+v_field(2,2,:).^2+v_field(3,2,:).^2;   
    
    ab = v_field(1,1,:).*v_field(1,2,:)+v_field(2,1,:).*v_field(2,2,:)+v_field(3,1,:).*v_field(3,2,:);
    
    r = (  2*(b(:).^2)./( (asq + bsq - sqrt( ab.^2.*4+(asq-bsq).^2 )) )  ).^(1./6);
    
    temp =  r.^6./(3*b.^2);
    
    u = sqrt(temp.*asq - 1/3 );
    v = sign(ab) .* sqrt(temp.*bsq - 1/3 );
    w = real( sqrt( 5/3 - temp.*(asq+bsq) ) );
    
        

    pos.r = reshape(r,sz);
    pos.x = reshape(r.*u,sz);
    pos.y = reshape(r.*v,sz);
    pos.z = reshape(r.*w,sz);
    pos.rotmat = real(reshape(CalcRotMat(v_field(1,1,:), v_field(2,1,:), v_field(3,1,:), v_field(1,2,:), v_field(2,2,:), v_field(3,2,:), u, v, w, r, b), [3,3,sz]));
    
end