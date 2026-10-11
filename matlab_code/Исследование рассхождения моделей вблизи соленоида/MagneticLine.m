classdef MagneticLine
    %MAGNETICLINE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        B1;
        B2;
        B3;
        alf (1,1) double;
        bet (1,1) double;
        gam (1,1) double;
    end
    
    methods
        function obj = MagneticLine(alf, bet, gam, L, zeroPoint, ro)
            arguments
                alf (1, 1) double;
                bet (1, 1) double;
                gam (1, 1) double;
                L (1, 1) double;
                zeroPoint (1, 3) double;
                ro (1,1) double;
            end
            obj.alf = alf;
            obj.bet = bet;
            obj.gam = gam;
            L0 = -L/2;
            R =...
               [cos(gam ),   -sin(gam ),  0;...
                sin(gam ),   cos(gam ),   0;...
                0,          0,          1]*...
               [cos(bet ),   0,          sin(bet );...
                0,          1,          0;...
                -sin(bet ),  0,          cos(bet );]*...
               [1,          0,          0;...
                0,          cos(alf ),   -sin(alf );...
                0,          sin(alf ),   cos(alf );];
            mynorm = @(r) (sqrt(sum(r.^2)));

            syms l real positive
            m0 = 4e-7*pi;
            syms x y z real positive
            
            
            r = (R*[x-zeroPoint(1),y-zeroPoint(2),z-zeroPoint(3)]')'; % Радиус-вектор до точки
            rm = [l, 0, 0]; %Радиус-вектор до элменета катушки            
            vec_ro = [ro, 0, 0];
            syms bx by bz real
            dB = (m0/(4*pi)) * (3*(r-rm)*dot(vec_ro,(r-rm))/mynorm(r-rm)^5-vec_ro/mynorm(r-rm)^3);
            %Часть интегрирование Bx
            fn = 1/((l-x)^2+(y)^2+(z)^2)^sym(3/2);
            half_Bx1(l) = int(fn,l);
            fn = (l-x)*(l-x)/(mynorm([l-x,y,z]))^5;
            half_Bx2(l) = 3*int(fn,l);


            %Заканчиваем интегрировать Bx
            Bx = (m0*ro/(4*pi))*((half_Bx2(L0+L)-half_Bx2(L0))-(half_Bx1(L0+L)-half_Bx1(L0)));
            obj.B1 = matlabFunction(subs(simplify(Bx), m0, 4e-7*pi));

            %Интегрирование By 
            By = simplify(subs(int(dB(2), l), l, L0+L) - subs(int(dB(2), l), l, L0));
            obj.B2 = matlabFunction(simplify(By));

            %Интегрирование Bz
            Bz = simplify(subs(int(dB(3),l),l,L0+L) - subs(int(dB(3),l),l,L0));
            obj.B3 = matlabFunction(simplify(Bz));
        end
        
        function [bx, by, bz] = Field(obj, x, y, z)
            arguments 
                obj;
                x (1,:) double;
                y (1,:) double;
                z (1,:) double;
            end            
            vec_b = [obj.B1(x,y,z); obj.B2(x,y,z); obj.B3(x,y,z)];
            R =...
               [cos( obj.gam ),   -sin( obj.gam ),  0;...
                sin( obj.gam ),   cos( obj.gam ),   0;...
                0,          0,          1]*...
               [cos( obj.bet ),   0,          sin( obj.bet );...
                0,          1,          0;...
                -sin( obj.bet ),  0,          cos( obj.bet );]*...
               [1,          0,          0;...
                0,          cos( obj.alf ),   -sin( obj.alf );...
                0,          sin( obj.alf ),   cos( obj.alf );];
            for i = 1:size(vec_b, 2)
                vec_b(:, i) = vec_b(:,i);%R\vec_b(:,i); 
            end
            bx = -1.*vec_b(1, :);
            by = vec_b(2, :);
            bz = vec_b(3, :);
        end
        
    end
end

