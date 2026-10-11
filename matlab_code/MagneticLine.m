classdef MagneticLine
    %MAGNETICLINE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        B1;
        B2;
        B3;
        alf (1,1) double;
        bet (1,1) double;
    end
    
    methods
        function obj = MagneticLine(alf, bet, L, zeroPoint, ro)
            arguments
                alf (1, 1) double;
                bet (1, 1) double;
                L (1, 1) double;
                zeroPoint (1, 3) double;
                ro (1,1) double;
            end
            obj.alf = alf;
            obj.bet = bet;
            L0 = -L/2;
            R = [cos(alf)*cos(bet), sin(alf)*cos(bet), sin(bet);...
                 sin(alf),          cos(alf),          0;...
                 -cos(alf)*sin(bet),-sin(alf)*sin(bet),cos(bet)];
            mynorm = @(r) (sqrt(sum(r.^2)));

            syms l real positive
            m0 = 4e-7*pi;
            syms x y z real positive
            r = (R*[x-zeroPoint(1), y-zeroPoint(2), z-zeroPoint(3)]')'; % Радиус-вектор до точки
            rm = [l, 0, 0]; %Радиус-вектор до элменета катушки            
            vec_ro = [ro, 0, 0];
            syms bx by bz real
            dB = (m0/(4*pi)) * (3*(r-rm)*dot(vec_ro,(r-rm))/mynorm(r-rm)^5-vec_ro/mynorm(r-rm)^3);
            %Часть интегрирование Bx
            fn = 1/((l-x)^2+(y)^2+(z)^2)^sym(3/2);
            half_Bx1(l) = int(fn,l);
            fn = (l-x)*(l-x)/(mynorm([l-x,y,z]))^5;
            half_Bx2(l) = 3*int(fn,l);


            %disp(num2str(i) + " - ro = " + num2str(ro) + "[А*м]" + " | L = " + num2str(L) + "[м]") 

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
            inv_R = [ (cos(obj.alf)*cos(obj.bet))/(2*cos(obj.alf)^2 - 1), sin(obj.alf)/(2*sin(obj.alf)^2 - 1), -(cos(obj.alf)*sin(obj.bet))/(2*cos(obj.alf)^2 - 1);...
                    -(cos(obj.bet)*sin(obj.alf))/cos(2*obj.alf), cos(obj.alf)/(2*cos(obj.alf)^2 - 1), -(sin(obj.alf)*sin(obj.bet))/(2*sin(obj.alf)^2 - 1);...
                    sin(obj.bet),                       0,                            cos(obj.bet)];
            for i = 1:size(vec_b, 2)
                vec_b(:,i) = inv_R*vec_b(:,i);      
            end
            bx = vec_b(1,:);
            by = vec_b(2,:);
            bz = vec_b(3,:);
        end
        
    end
end

