function [alf,bet,gam] = CalcAngles(rotMatrix)

    arguments
        rotMatrix (3,3) {mustBeNumeric, mustBeReal};
    end
    
    % alf = atan2(r32, r33)
    % bet = atan2(-r31, sqrt(r32^2 + r33^2))
    % gam = atan2(r21, r11)
    
    for i = 1:size(rotMatrix, 3)
        bet = -asin(rotMatrix(3,1));

        if abs(bet) ~= pi/2
            alf1=asin(rotMatrix(3,2)/sqrt(1-rotMatrix(3,1)^2));
            alf2=pi*sign(rotMatrix(3,2))-asin(rotMatrix(3,2)/sqrt(1-rotMatrix(3,1)^2));
            if rotMatrix(3,3) >= 0
                alf = alf1;
            else
                alf = alf2;
            end

            if rotMatrix(1,1) >= 0
                gam = asin(rotMatrix(2,1)/sqrt(1-rotMatrix(3,1)^2));
            else
                gam = pi*sign(rotMatrix(2,1))-asin(rotMatrix(2,1)/sqrt(1-rotMatrix(3,1)^2));
            end
            else
            alf = 0;
            if rotMatrix(2,2) >= 0
                gam = -asin(rotMatrix(1,2));
            else
                gam = -pi*sign(rotMatrix(1,2))+asin(rotMatrix(1,2));
            end
        end
    end
end

