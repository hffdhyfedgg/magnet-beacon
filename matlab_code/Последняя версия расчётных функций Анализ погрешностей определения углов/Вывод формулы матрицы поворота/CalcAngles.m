function angles = CalcAngles(rotMatrix)

    arguments
        rotMatrix (3,3) {mustBeNumeric, mustBeReal};
    end
    
    % alf = atan2(r32, r33)
    % bet = atan2(-r31, sqrt(r32^2 + r33^2))
    % gam = atan2(r21, r11)
    

        angles.bet = -asin(rotMatrix(3,1));

        if abs(angles.bet) ~= pi/2

            if rotMatrix(3,3) >= 0
                angles.alf = asin(rotMatrix(3,2)/sqrt(1-rotMatrix(3,1)^2));
            else
                angles.alf = pi*sign(rotMatrix(3,2))-asin(rotMatrix(3,2)/sqrt(1-rotMatrix(3,1)^2));
            end

            if rotMatrix(1,1) >= 0
                angles.gam = asin(rotMatrix(2,1)/sqrt(1-rotMatrix(3,1)^2));
            else
                angles.gam = pi*sign(rotMatrix(2,1))-asin(rotMatrix(2,1)/sqrt(1-rotMatrix(3,1)^2));
            end
        else
            angles.alf = 0;
            if rotMatrix(2,2) >= 0
                angles.gam = -asin(rotMatrix(1,2));
            else
                angles.gam = -pi*sign(rotMatrix(1,2))+asin(rotMatrix(1,2));
            end
        end
    
end

