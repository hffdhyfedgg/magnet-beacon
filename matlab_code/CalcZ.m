function z = CalcZ(B,ax, ay, az, bx, by, bz)
arguments
    B double {mustBeNumeric, mustBeReal};
    ax double {mustBeNumeric, mustBeReal};
    ay double {mustBeNumeric, mustBeReal};
    az double {mustBeNumeric, mustBeReal};
    bx double {mustBeNumeric, mustBeReal};
    by double {mustBeNumeric, mustBeReal};
    bz double {mustBeNumeric, mustBeReal};
end
    out = zeros(1,numel(bz));
    for i = 1:numel(bz)
        out(i) = (2.0^(1.0/6.0)*sqrt(3.0)*(1.0/(-sqrt((ax(i)*bx(i)+ay(i)*by(i)+az(i)*bz(i))^2*4.0+(ax(i)^2+ay(i)^2+az(i)^2-bx(i)^2-by(i)^2-bz(i)^2)^2)+ax(i)^2+ay(i)^2+az(i)^2+bx(i)^2+by(i)^2+bz(i)^2))^(1.0/6.0)*(B(i)^2)^(1.0/6.0)*sqrt((sqrt((ax(i)*bx(i)+ay(i)*by(i)+az(i)*bz(i))^2*4.0+(ax(i)^2+ay(i)^2+az(i)^2-bx(i)^2-by(i)^2-bz(i)^2)^2)*-5.0+ax(i)^2*3.0+ay(i)^2*3.0+az(i)^2*3.0+bx(i)^2*3.0+by(i)^2*3.0+bz(i)^2*3.0)/(-sqrt((ax(i)*bx(i)+ay(i)*by(i)+az(i)*bz(i))^2*4.0+(ax(i)^2+ay(i)^2+az(i)^2-bx(i)^2-by(i)^2-bz(i)^2)^2)+ax(i)^2+ay(i)^2+az(i)^2+bx(i)^2+by(i)^2+bz(i)^2)))/3.0;
        %{
        if ~isreal(out(i))
            out(i) = 0;
        end
        %}
    end
    z = reshape(out,size(bz));
end

