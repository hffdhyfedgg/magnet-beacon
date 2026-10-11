function [a,b] = Calcab(A,phi)
    arguments
        A (1,1) double {mustBeReal};
        phi (1,1) double {mustBeReal};
    end
    a = A*cos(phi);
    b = -A*sin(phi);
end

