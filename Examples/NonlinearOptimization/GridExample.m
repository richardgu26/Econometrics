function GridExample(points, doprint=false)
        close all;

        % plot the line
        x = linspace(-pi,pi/2,1000);
        obj = @(theta) 0.5*theta.*sin(theta.^2);
        y = obj(x);
        plot(x, y);
        hold on;

        % plot the grid points
        x = linspace(-pi,pi/2,points);
        y = obj(x);
        plot(x, y, 'ob', 'Markersize', 4);

        % plot the best point found
        [junk ind] = min(y);
        plot(x(ind), y(ind), 'or', 'Markersize', 8);
        hold off;

        if doprint print('gridsearch.svg', '-dsvg'); endif
endfunction

