close all
clear
clc


% a)
x0 = -1;                    % initial condition for x
y0 = -1;                    % Initial condition for y




h = 0.005;                                           % Change of the step-size
% d)
fprintf('Analysis of the membrane potential for 5 different initial currents, with constant step size')
t0 = [0 1000];                                          % New range for the time
r = 0.005;
Ik = (1:0.1:4);
for i = 1:length(Ik)
    if i ~= 1
        [x,y,z,t] = setO(x,y,z,t);                          % Erasing of the vectors in which the solution is stored
    end
    I = Ik(i);
    z0 = I;
    [fx,fy,fz] = funcs(I,r);                            % Call for the functions of the system
    [x,y,z,t] = rk4sys3(fx,fy,fz,x0,y0,z0,t0,h);        % Solution of the Runge-Kutta 4th order model
    figure(2)
    plot(t,x)
    title('Membrane potential for I =',I)
    if i == 1
        fprintf('\nFor a current of %.1f, the membrane potential is quiescent', I)
    elseif i == 2
        fprintf('\nFor a current of %.1f, the membrane potential consists of periodic firing of spikes', I)
    elseif i == 3
        fprintf('\nFor a current of %.1f, the membrane potential consists of periodic bursts of spikes', I)
    else
        fprintf('\nFor a current of %.1f, the membrane potential consists of periodic firing of spikes', I)
    end
    pause
end

%%

x = [1 2 3 4 10 100 10100];
y = (1:7);
loglog(x,y)