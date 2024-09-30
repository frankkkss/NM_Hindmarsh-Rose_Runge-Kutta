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