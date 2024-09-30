function [] = tplots(x,y,z,t)
    subplot(3,1,1)
    hold on
    plot(t,x)
    ylabel('Membrane potential')
    xlabel('Time')
    ylim([-1 2])
    
    subplot(3,1,2)
    hold on
    plot(t,y)
    ylabel('Fast ion channels exchange')
    xlabel('Time')
    ylim([-6 1])    
    
    subplot(3,1,3)
    hold on
    plot(t,z)
    ylabel('Slow ion channels exchange')
    xlabel('Time')
    ylim([3.4 4])
end