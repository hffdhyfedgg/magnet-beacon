classdef MagneticMoment
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0 (1,1) double;
        y0 (1,1) double;
        z0 (1,1) double;
        mx (1,1) double;
        my (1,1) double;
        mz (1,1) double;
    end
    
    methods
        function obj = MagneticMoment(x0, y0, z0, mx, my, mz)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.x0 = x0;
            obj.y0 = y0;
            obj.z0 = z0;
            obj.mx = mx;
            obj.my = my;
            obj.mz = mz;            
            
        end
        
        function [bx, by, bz] = Field(obj, x,y,z)
            arguments
                obj;
                x (1,:) double;
                y (1,:) double;
                z (1,:) double;
            end
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            var1 = zeros(size(x));
            var2 = zeros(size(x));
            var3 = zeros(size(x));
            mynorm = @(r) (sqrt(sum(r.^2)));
            for i = 1:size(x,2)
                var1(i) = 1e-7*(3*(x(i)-obj.x0)*(obj.mx*(x(i)-obj.x0) + obj.my*(y(i)-obj.y0) + obj.mz*(z(i)-obj.z0))/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^5 - obj.mx/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^3);
                var2(i) = 1e-7*(3*(y(i)-obj.x0)*(obj.mx*(x(i)-obj.x0) + obj.my*(y(i)-obj.y0) + obj.mz*(z(i)-obj.z0))/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^5 - obj.my/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^3);
                var3(i) = 1e-7*(3*(z(i)-obj.x0)*(obj.mx*(x(i)-obj.x0) + obj.my*(y(i)-obj.y0) + obj.mz*(z(i)-obj.z0))/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^5 - obj.mz/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^3);
            end
            bx = var1;
            by = var2;
            bz = var3;
        end
        function b = AbsField(obj, x,y,z)
            arguments
                obj;
                x (1,:) double;
                y (1,:) double;
                z (1,:) double;
            end
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            mynorm = @(r) (sqrt(sum(r.^2)));
            for i = 1:size(x,2)
                var1 = 1e-7*(3*(x(i)-obj.x0)*(obj.mx*(x(i)-obj.x0) + obj.my*(y(i)-obj.y0) + obj.mz*(z(i)-obj.z0))/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^5 - obj.mx/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^3);
                var2 = 1e-7*(3*(y(i)-obj.x0)*(obj.mx*(x(i)-obj.x0) + obj.my*(y(i)-obj.y0) + obj.mz*(z(i)-obj.z0))/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^5 - obj.my/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^3);
                var3 = 1e-7*(3*(z(i)-obj.x0)*(obj.mx*(x(i)-obj.x0) + obj.my*(y(i)-obj.y0) + obj.mz*(z(i)-obj.z0))/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^5 - obj.mz/mynorm([x(i)-obj.x0, y(i)-obj.y0, z(i)-obj.z0])^3);
                var4 = sqrt(var1.^2 + var2.^2 + var3.^2);
            end
            b = var4;
        end
    end
end

