function [f, g, l] = fit_int_autocorr(t,h)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, h );

% Set up fittype and options.
ft = fittype( 'a*exp(-x/b1)+(1-a)*exp(-x/b2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = [0.25 .25 3];

% Fit model to data.
[f, g] = fit( xData, yData, ft, opts );

% Integrate to get hbond lifetimes
l = (1-f.a).*(f.b2) + (f.a).*(f.b1);
end




