function [x,y,z,t] = rk4sys3(fx,fy,fz,x,y,z,t,h)

    N = ceil(t(end)/h);             % Number of iterative steps

    for i = 1:N
        t(i+1) = t(i) + h;      % Definition of time vector (although it is not used)

        % First Runge-Kutta parameter
        k1x = h*fx(t(i),x(i),y(i),z(i));              
        k1y = h*fy(t(i),x(i),y(i),z(i));
        k1z = h*fz(t(i),x(i),y(i),z(i));

        % Second Runge-Kutta parameter
        k2x = h*fx(t(i) + h/2, x(i) + k1x/2, y(i) + k1y/2, z(i) + k1z/2);           
        k2y = h*fy(t(i) + h/2, x(i) + k1x/2, y(i) + k1y/2, z(i) + k1z/2);
        k2z = h*fz(t(i) + h/2, x(i) + k1x/2, y(i) + k1y/2, z(i) + k1z/2);
        
        % Third Runge-Kutta parameter
        k3x = h*fx(t(i) + h/2, x(i) + k2x/2, y(i) + k2y/2, z(i) + k2z/2);           
        k3y = h*fy(t(i) + h/2, x(i) + k2x/2, y(i) + k2y/2, z(i) + k2z/2);
        k3z = h*fz(t(i) + h/2, x(i) + k2x/2, y(i) + k2y/2, z(i) + k2z/2);

        % Fourth Runge-Kutta parameter
        k4x = h*fx(t(i) + h, x(i) + k3x, y(i) + k3y, z(i) + k3z);               
        k4y = h*fy(t(i) + h, x(i) + k3x, y(i) + k3y, z(i) + k3z);
        k4z = h*fz(t(i) + h, x(i) + k3x, y(i) + k3y, z(i) + k3z);

        % Solution vectors of the system
        x(i+1) = x(i) + 1/6*(k1x + 2*k2x + 2*k3x + k4x);             
        y(i+1) = y(i) + 1/6*(k1y + 2*k2y + 2*k3y + k4y);
        z(i+1) = z(i) + 1/6*(k1z + 2*k2z + 2*k3z + k4z);
    end

end
