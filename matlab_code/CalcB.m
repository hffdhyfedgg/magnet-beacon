function output = CalcB(A, phi)
    arguments
        A {mustBeReal, mustBePositive, mustBeNonzero}
        phi {mustBeReal}
    end
    var = zeros(size(phi));
    for i = 1:max(size(phi))
        if abs(phi(i)) > pi
            error("Аргумент phi может быть равен только числам в диапазоне от -pi до pi");
        end
        
        if phi(i) < 0
            var(i) = A(i)*sin(phi);
        else
            var(i) = A(i)*sin(phi); 
        end
    end
    output = var;
end
