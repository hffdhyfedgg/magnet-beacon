function x = CalcX(b,v)
arguments
    b double {mustBeNumeric, mustBeReal};
    v (3,2, :) {mustBeNumeric, mustBeReal};
end

    asq = v(1,1,:).^2+v(2,1,:).^2+v(3,1,:).^2;
    bsq = v(1,2,:).^2+v(2,2,:).^2+v(3,2,:).^2;    

    sz = [size(v), 1];
%     x =  reshape( (2.0.^(1.0./6.0).*sqrt(3.0).*(1.0./(-sqrt((v(1,1,:).*v(1,2,:)+v(2,1,:).*v(2,2,:)+v(3,1,:).*v(3,2,:)).^2.*4.0+(asq-bsq).^2)+asq+bsq)).^(1.0./6.0).*(B(i).^2).^(1.0./6.0).*sqrt(((asq).*2.0)./(-sqrt((v(1,1,:).*v(1,2,:)+v(2,1,:).*v(2,2,:)+v(3,1,:).*v(3,2,:)).^2.*4.0+(asq-bsq).^2)+asq+bsq)-1.0))./3.0, sz(3:end)  );
    x =  reshape( (2.0.^(1.0./6.0).*sqrt(3.0).*(1.0./(-sqrt((v(1,1,:).*v(1,2,:)+v(2,1,:).*v(2,2,:)+v(3,1,:).*v(3,2,:)).^2.*4.0+(asq-bsq).^2)+asq+bsq)).^(1.0./6.0).*(b(:).^2).^(1.0./6.0).*sqrt(((asq).*2.0)./(-sqrt((v(1,1,:).*v(1,2,:)+v(2,1,:).*v(2,2,:)+v(3,1,:).*v(3,2,:)).^2.*4.0+(asq-bsq).^2)+asq+bsq)-1.0))./3.0, sz(3:end)  );
end
  


% function x = CalcX(b,v)
% arguments
%     b double {mustBeNumeric, mustBeReal};
%     v (3,2, :) {mustBeNumeric, mustBeReal};
% end
%     out = zeros(1,numel(bz));
%     for i = 1:numel(bz)
%         out(i) = (2.0^(1.0/6.0)*sqrt(3.0)*(1.0/(-sqrt((ax(i)*bx(i)+ay(i)*by(i)+az(i)*bz(i))^2*4.0+(ax(i)^2+ay(i)^2+az(i)^2-bx(i)^2-by(i)^2-bz(i)^2)^2)+ax(i)^2+ay(i)^2+az(i)^2+bx(i)^2+by(i)^2+bz(i)^2))^(1.0/6.0)*(B(i)^2)^(1.0/6.0)*sqrt(((ax(i)^2+ay(i)^2+az(i)^2)*2.0)/(-sqrt((ax(i)*bx(i)+ay(i)*by(i)+az(i)*bz(i))^2*4.0+(ax(i)^2+ay(i)^2+az(i)^2-bx(i)^2-by(i)^2-bz(i)^2)^2)+ax(i)^2+ay(i)^2+az(i)^2+bx(i)^2+by(i)^2+bz(i)^2)-1.0))/3.0;
%     end
%     x = reshape(out,size(bz));
% end


