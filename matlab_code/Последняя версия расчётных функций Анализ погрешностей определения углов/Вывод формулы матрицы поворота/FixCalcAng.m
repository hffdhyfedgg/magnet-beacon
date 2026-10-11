function angles = FixCalcAng(rotMatrix)

arguments
    rotMatrix (3,3) {mustBeNumeric, mustBeReal};
end
    angles.alf = zeros(1, size(rotMatrix, 3));
    angles.bet = zeros(1, size(rotMatrix, 3));
    angles.gam = zeros(1, size(rotMatrix, 3));
    for i = 1:size(rotMatrix, 3)
%         angles.alf(i) = atan2(rotMatrix(3,2,i), rotMatrix(3,3,i));
%         angles.bet(i) = atan2(-rotMatrix(3,1,i), sqrt(rotMatrix(3,2,i)^2 + rotMatrix(3,3,i)^2));
%         angles.gam(i) = atan2(rotMatrix(2,1,i), rotMatrix(1,1,i));
        angles.gam(i) = atan2(rotMatrix(2,1,i), rotMatrix(1,1,i));
        angles.bet(i) = atan2(-rotMatrix(3,1,i), sqrt(rotMatrix(3,2,i)^2 + rotMatrix(3,3,i)^2));
        angles.alf(i) = atan2(rotMatrix(3,2,i), rotMatrix(3,3,i));
    end

end

