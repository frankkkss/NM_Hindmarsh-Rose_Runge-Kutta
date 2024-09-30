close all
clear
clc

 
%% Study of the effect of the step size in the solution of the system of differential equations
fprintf('Analysis of the Hindmarsh-Rose model via 4th order Runge-Kutta method for 6 different values of the step size')
I = 3.5;                    % Value of the current for the first cases of study
r = 0.005;                  % Value of the efficiency of the slow ion channels
t0 = [0 300];               % Time range for study
h = 0.001;                  % Step-size at which the system is evaluated via Runge-Kutta 4th order
x0 = -1;                    % initial condition for x
y0 = -1;                    % Initial condition for y
z0 = I;                     % Initial condition for z

[fx,fy,fz] = funcs(I,r);                            % Call for the functions of the system
n = 1;
[x,y,z,t] = rk4sys3(fx,fy,fz,x0,y0,z0,t0,h);        % Solution of the Runge-Kutta 4th order model
                                          
figure(1)
set(gcf, 'position', [0,0,1500, 1800])

hv = [0.005 0.025 0.01  0.1 0.4 0.5];
for m = 1:length(hv)
    h = hv(m);
    if m ~= 1    
        [x,y,z,t] = setO(x,y,z,t);                          % Erasing of the vectors in which the solution is stored
    end
    [x,y,z,t] = rk4sys3(fx,fy,fz,x0,y0,z0,t0,h);        % Solution of the Runge-Kutta 4th order model
    n = m + 1;
    tplots(n,x,y,z,t,hv);                                % Function for plotting the values from the solution fo the Runge-Kutta model
end

fprintf('\n\nThe first four solutions are exactly the same so they are superimposed, and only the last one is represented')
fprintf('\nWith a stepsize up to %.1f, there is no perceptible change in the solutions. \nFor %.1f there is a distortion that makes the spikes to be irregular and have no defined period. \nFor %.1f the solution blows off to infinity.\n\n', hv(3), hv(4), hv(5)')
pause
close all


h = 0.005;                                              % Change of the step-size, wuill remain constant for the rest


%% Study of the shape of membrane potential as a function of the current
fprintf('Analysis of the membrane potential for 5 different initial currents, with constant step size')
t0 = [0 1000];                                          % New range for the time
r = 0.005;
Ik = [1.2 1.4 2.5 3.1 3.5];
for i = 1:length(Ik)
    [x,y,z,t] = setO(x,y,z,t);                          % Erasing of the vectors in which the solution is stored
    I = Ik(i);
    z0 = I;
    [fx,fy,fz] = funcs(I,r);                            % Call for the functions of the system
    [x,y,z,t] = rk4sys3(fx,fy,fz,x0,y0,z0,t0,h);        % Solution of the Runge-Kutta 4th order model
    figure(2)
    set(gcf, 'position', [0,0,1500, 1800])
    plot(t,x)
    title('Membrane potential for I =',I)
    if i == 1
        fprintf('\nFor a current of %.1f, the membrane potential is quiescent as there is just one spike', I)
    elseif i == 2
        fprintf('\nFor a current of %.1f, the membrane potential consists of periodic firing of spikes', I)
    elseif i == 3
        fprintf('\nFor a current of %.1f, the membrane potential consists of periodic bursts of spikes', I)
    elseif i == 4
        fprintf('\nFor a current of %.1f, the membrane potential consists of irregular firings of bursts of spikes', I')
    else
        fprintf('\nFor a current of %.1f, the membrane potential consists of periodic firing of spikes', I)
    end
    pause
end

close all



%% Detection of peaks in the membrane potential for different values of r
fprintf('\n\nDetection of the peaks in the membrane potential for 5 different values of r')
rk = [0.005 0.01 0.025 0.05 0.1];
I = 3.1;
z0 = I;                                                 % Setting the initial condition for z
x_th = -0.4;                                            % Threshold value for the finding of peaks
for j = 1:length(rk)
    [x,y,z,t] = setO(x,y,z,t);                          % Erasing of the vectors in which the solution is stored
    r = rk(j);                                          % Setting the value of r from the proposed values to study
    [fx,fy,fz] = funcs(I,r);                            % Call for the functions of the system to include the newest values of the parameters
    [x,y,z,t] = rk4sys3(fx,fy,fz,x0,y0,z0,t0,h);        % Solution of the Runge-Kutta 4th order model
    [t_spi, spi] = findspikes(x_th, x, h);              % Function for finding the spikes in the membrane potential from a previously stated value x_th
    figure(3)
    set(gcf, 'position', [0,0,1500, 1800])
    plot(t,x)                                           % Plot of the membrane potential
    hold on
    plotspike(t_spi, spi);                              % Plot for the spikes found in the membrane potential
    title('Membrane potential for r =',r)
    ylabel('Membrane potential')
    xlabel('t')
    pause
    close all
end



%% Study of the bifurcation diagram for t = [0 1000] and gamma = 1.05
gamma = 1.05;
p = floor(log(0.1/0.005)/log(gamma));                   % Calculation of the length of the r vector
b = 1000;                                               % Time boundary
t0 = [0 b];                                             % Time range for study

fprintf('\n\nBifurcation diagram for values of r in the sequence 0.005*gamma^p in the time interval [0 %d]\n', b)
figure(4)
set(gcf, 'position', [0,0,1500, 1800])
hold on
for k = 0:p
    r = 0.005*gamma^k;                                  % Value of r at each iteration
    [x,y,z,t] =setO(x,y,z,t);                           % Erasing of the vectors in which the solution is stored
    [fx,fy,fz] = funcs(I,r);                            % Call for the functions of the system
    [x,y,z,t] = rk4sys3(fx,fy,fz,x0,y0,z0,t0,h);        % Solution of the Runge-Kutta 4th order model
    [t_spi, spi] = findspikes(x_th, x, h);              % Function for finding the spikes in the membrane potential from a previously stated value x_th
    isv = [];
    isv = diff(t_spi);
    loglog(r,isv,'b.')
end

pause
close all
fprintf('\n\nThis process lasts approxiamtely 15 seconds\n')
%% Bifurcation diagram for longer computation
fprintf('For reaching a time of 2h or 7200s, whe have 15s*f_t*f_p, where f_t = b/1000 and f_p = p/61 (as the previous value of p is 61)\n')
fprintf('Appropriate values for a computation time near 2h would be a time interval of t = [0 10000] and around 2900 values of r, which means a value of gamma = 1.00103')
fprintf('\nMy computer cannot support such amounts of data for representation so the initial conditions selected for study are b = 10000 and gamma = 1.002')

tic
gamma = 1.001;
p = floor(log(0.1/0.005)/log(gamma));                   % Calculation of the length of the r vector
b = 10000;                                              % Time boundary
t0 = [0 b];                                             % Time range for study

%fprintf('\n\nBifurcation diagram for extended time of study and lower value of gamma (more values of r and peaks are looked upon)')
figure(5)
set(gcf, 'position', [0,0,1500, 1800])
hold on
for k = 0:p
    r = 0.005*gamma^k;
    [x,y,z,t] =setO(x,y,z,t);                           % Erasing of the vectors in which the solution is stored
    [fx,fy,fz] = funcs(I,r);                            % Call for the functions of the system
    [x,y,z,t] = rk4sys3(fx,fy,fz,x0,y0,z0,t0,h);        % Solution of the Runge-Kutta 4th order model
    [t_spi, spi] = findspikes(x_th, x, h);              % Function for finding the spikes in the membrane potential from a previously stated value x_th
    isv = [];
    isv = diff(t_spi);
    loglog(r,isv,'b.')
end
toc
pause
close all

fprintf('\n\nAnalysis of the spikes in the membrane potential for 3 values of r')

%% Plot of the membrane potential of 3 values of r for t= [0 1000]
t0 = [0 1000];
rk = [0.005 0.025 0.88];
figure(6)
set(gcf, 'position', [0,0,1500, 1800])
hold on 
for o = 1:length(rk)
    r = rk(o);
    [x,y,z,t] =setO(x,y,z,t);                           % Erasing of the vectors in which the solution is stored
    [fx,fy,fz] = funcs(I,r);                            % Call for the functions of the system
    [x,y,z,t] = rk4sys3(fx,fy,fz,x0,y0,z0,t0,h);        % Solution of the Runge-Kutta 4th order model
    [t_spi, spi] = findspikes(x_th, x, h);              % Function for finding the spikes in the membrane potential from a previously stated value x_th
    subplot(3,1,o)
    plot(t,x)
    xlim(t0);
end
fprintf('\nIn the bifurcation diagram, many peaks could be seen for r=0.005, although much more inter spike intervals were plotted as the interval for study was t=[0,2000]')
fprintf('\nFor r=0.025, the peaks are very regular following two frequencies, as it could be seen in the bifurcation diagram (if a line r=0.0025 was drawn, there would be just two intersections)')
fprintf('\nFor the bigger values of r, the was little separation in the inter spike distance, and there is just one spike in the selected interval.')


%%
function [] = tplots(n,x,y,z,t,h)
    subplot(3,1,1)
    hold on
    plot(t,x)
    ylabel('Membrane potential')
    xlabel('Time')
    ylim([-1 2])
    xlim([0 300])
    for u = 1:length(h)
        leg(u) = 'h = ' + string(h(u));
    end       
    legend(leg(:), 'Orientation', 'horizontal', 'Location', 'best')
    
    subplot(3,1,2)
    hold on
    plot(t,y)
    ylabel('Fast ion channels exchange')
    xlabel('Time')
    ylim([-6 1])
    xlim([0 300])
    legend(leg(:), 'Orientation', 'horizontal', 'Location', 'best')
    
    subplot(3,1,3)
    hold on
    plot(t,z)
    ylabel('Slow ion channels exchange')
    xlabel('Time')
    ylim([3.4 4])
    xlim([0 300])
    legend(leg(:), 'Orientation', 'horizontal', 'Location', 'best')
end

%%
function [x,y,z,t] = setO(x,y,z,t)
    x = [];
    y = [];
    z = [];
    t = [];
end

%% 
function [t_spi, spi] = findspikes(x_th, x, h)
    t_spi = [];
    pos_spi = [];
    spi = [];

    for j = 3:(length(x)-2)
        if (x(j) > x(j-1)) && (x(j) > x(j-2)) && (x(j) > x(j+1)) && (x(j) > x(j+2)) && (x(j) > x_th)
            pos_spi = [pos_spi j];
            spi = [spi x(j)];
            t_spi = [t_spi h*j];
        end
    end
end

%% 
function [] = plotspike(t_spi, spi)
    hold on
    for e = 1:length(spi)
        plot(t_spi(e), spi(e),'xg')

    end
end

%%
function [fx,fy,fz] = funcs(I,r)
    fx = @(t,x,y,z) y + x*(3*x - x^2) - z + I;          
    fy = @(t,x,y,z) 1 - 5*x^2 - y;
    fz = @(t,x,y,z) r*(4*(x + 8/5) -z);
end