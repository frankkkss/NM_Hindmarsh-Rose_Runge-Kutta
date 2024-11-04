function [t_spi, spi] = findspikes(x_th, x, h)
    % Finding spikes along the electrical current
    % Inputs:
    % x_th: threshold value for spike detection
    % x: electrical current
    % h: time step
    %
    % Outputs:
    % t_spi: vector containng the timestamps of the spikes
    % spi: vector containing the height of the spikes
    

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