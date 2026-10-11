function generator = CreateABParameterGenerator(magneticObjects)
%{
    {'l', 'cos', L, [0, 0, 0]}
    {'m', 'sin', [n1, n2, n3], [0,0,0]}
%}
arguments
    magneticObjects (:,:) cell; 
end
mynorm = @(r) (sqrt(sum(r.^2)));
resultCos = repmat(str2sym('0'), [size(magneticObjects, 1), 3]);
resultSin = repmat(str2sym('0'), [size(magneticObjects, 1), 3]);
for i = 1:size(magneticObjects, 1)
    switch magneticObjects{i, 1}
        case "l"
            
            if numel(magneticObjects{i, 4}) ~= 3 || norm(magneticObjects{i, 4}) ~= 1
                error("Не подходящий вектор направленности соленоида.")
            end
            
            l = sym('l', {'real'});           
            rm = l.*[magneticObjects{i, 4}(1), magneticObjects{i, 4}(2), magneticObjects{i, 4}(3)];%Параметрически заданная кривая, описывающая соленоид
            x = sym('x',{'real'});
            y = sym('y',{'real'});
            z = sym('z',{'real'});
            r = [x, y, z];            
            Bmr = sym('Bmr', {'real'});
            rom = [magneticObjects{i, 4}(1), magneticObjects{i, 4}(2), magneticObjects{i, 4}(3)]; % Вектор магнитной плотности. Постоянен на всей длинне сленоида. 
            B1 = Bmr * (3*(r-rm)*dot(rom, (r-rm))/mynorm(r-rm)^5 - rom/mynorm(r-rm)^3);%Bmr = (ro*m0/(4*pi))
            
            b(1) = int( B1(1)*norm(diff(rm)), l);
            b(2) = int( B1(2)*norm(diff(rm)), l);
            b(3) = int( B1(3)*norm(diff(rm)), l);
            
            L = magneticObjects{i, 3};
            if magneticObjects{i, 2} == "cos"
                resultCos(i, 1) = subs(b(1), l, L/2) - subs(b(1), l, -L/2);
                resultCos(i, 2) = subs(b(2), l, L/2) - subs(b(2), l, -L/2);
                resultCos(i, 3) = subs(b(3), l, L/2) - subs(b(3), l, -L/2);
            end
            if magneticObjects{i, 2} == "sin"
                resultSin(i, 1) = subs(b(1), l, L/2) - subs(b(1), l, -L/2);
                resultSin(i, 2) = subs(b(2), l, L/2) - subs(b(2), l, -L/2);
                resultSin(i, 3) = subs(b(3), l, L/2) - subs(b(3), l, -L/2);
            end
        otherwise
            error("Не подходящий тип объекта в строке "+num2str(i)+" - " + magneticObjects(i, 1))
    end
end

alf = sym('alf', {'real'});
bet = sym('bet', {'real'});
gam = sym('gam', {'real'});

ax = sum(resultSin(:, 1));
ay = sum(resultSin(:, 2));
az = sum(resultSin(:, 3));
bx = sum(resultCos(:, 1));
by = sum(resultCos(:, 2));
bz = sum(resultCos(:, 3));

R = ...
   [cos(gam),   -sin(gam),  0;...
    sin(gam),   cos(gam),   0;...
    0,          0,          1]*...
   [cos(bet),   0,          sin(bet);...
    0,          1,          0;...
    -sin(bet),  0,          cos(bet);]*...
   [1,          0,          0;...
    0,          cos(alf),   -sin(alf);...
    0,          sin(alf),   cos(alf);];
generator = matlabFunction(simplify(R*[ax, bx; ay, by; az, bz]));
end

