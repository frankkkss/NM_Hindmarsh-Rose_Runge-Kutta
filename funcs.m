function [fx,fy,fz] = funcs(I,r)
    % Initialization of the Hindmarsh-Rose model equations
    % Inputs:
    % I: initial current parameter
    % r: slow ion channel efficiency parameter
    %
    % Outputs:
    % fx, fy, fz: functions of the Hindmarsh-Rose model equations
    fx = @(t,x,y,z) y + x*(3*x - x^2) - z + I;          
    fy = @(t,x,y,z) 1 - 5*x^2 - y;
    fz = @(t,x,y,z) r*(4*(x + 8/5) -z);
end