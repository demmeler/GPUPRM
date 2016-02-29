x = (0:10)';
y = sin(x);
xi = (0:.25:10)';
yi = interp1q(1:3',2:4',xi);
plot(x,y,'o',xi,yi)
