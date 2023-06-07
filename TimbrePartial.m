classdef TimbrePartial
    properties
        fy
        fx
        fim
        ay
        ax
        aim
    end 
    methods
        % Class Constructor
        function obj = TimbrePartial(fy_in,fx_in,fim_in,ay_in,ax_in,aim_in)
            if nargin == 0
                obj.fy = [1, 1];
                obj.fx = [0, 10];
                obj.fim = 'linear';
                obj.ay = [0, 0];
                obj.ax = [0, 10];
                obj.aim = 'linear';
            else
                obj.fy = fy_in;
                obj.fx = fx_in;
                obj.fim = fim_in;
                obj.ay = ay_in;
                obj.ax = ax_in;
                obj.aim = aim_in;
            end
        end
    end
end