function [fx,fy,fz] = funcs(I,r)
    fx = @(t,x,y,z) y + x*(3*x - x^2) - z + I;          
    fy = @(t,x,y,z) 1 - 5*x^2 - y;
    fz = @(t,x,y,z) r*(4*(x + 8/5) -z);
end