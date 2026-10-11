function [nx, ny, nz, ang] = AngNormal(rotMatrix)%Рассчитывает угол и орт поворота, по полученной матрице поворота
arguments
    rotMatrix (3,3)  double {mustBeNumeric, mustBeReal};
end

trace = rotMatrix(1,1) + rotMatrix(2,2) + rotMatrix(3,3);
theta = acos((trace - 1)/2);

nx = (rotMatrix(2,3) - rotMatrix(3,2))/(2*sin(theta));
ny = (rotMatrix(3,1) - rotMatrix(1,3))/(2*sin(theta));
nz = (rotMatrix(1,2) - rotMatrix(1,3))/(2*sin(theta));
ang = theta;
end

