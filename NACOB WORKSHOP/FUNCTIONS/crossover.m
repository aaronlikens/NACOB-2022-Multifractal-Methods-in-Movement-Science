function [ashort,along] = crossover(logn, logF, breakpoint, doplot)
x = logn;
y = logF;
y1 = y(1:breakpoint);
y2 = y(breakpoint+1:end);
m1 = polyfit(x(1:breakpoint), y1, 1);
p1 = polyval(m1, x(1:breakpoint));
m2 = polyfit(x(breakpoint+1:end), y2, 1);
p2 = polyval(m2, x(breakpoint+1:end));

if doplot
    plot(x,y,'o'); hold on;
    plot(x(1:breakpoint), p1,'r', 'LineWidth', 3);
    plot(x(breakpoint+1:end), p2,'r', 'LineWidth', 3);
    hold off;
    
end

% output slopes
ashort = m1(1);
along = m2(1);
end


