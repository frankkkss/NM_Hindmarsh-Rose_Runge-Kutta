function [] = plotspike(t_spi, spi)
    hold on
    for e = 1:length(spi)
        plot(t_spi(e), spi(e),'xg')
    end
end
