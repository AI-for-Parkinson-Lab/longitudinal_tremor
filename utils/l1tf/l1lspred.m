function ypred = l1lspred(xpred,x,yhat)
% Out-of-sample prediction for L1 regularized linear spline regression
% smoothing, using interpolation/extrapolation.
%
% Input arguments:
%  x: original covariate values
%  xpred: new covariate values
%  yhat: linear spline smoothed solution
%
% Output arguments:
%  ypred: new predicted, smoothed output values

% Input values should be in ascending order of covariates
[x,i] = sort(x);
yhat = yhat(i);

% Matlab's interp1 with linear interpolation is ready-made here
ypred = interp1(x,yhat,xpred,'linear','extrap');
