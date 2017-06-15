function plotlog(x, y, alpha, xrange)
    numx = 2;
    xv = logspace(log10(x), log10(x) + xrange, numx);
    C = y / (x.^alpha);
    yv = C * xv.^(alpha);
    plot(xv, yv,'k--', 'linewidth', 2)
end