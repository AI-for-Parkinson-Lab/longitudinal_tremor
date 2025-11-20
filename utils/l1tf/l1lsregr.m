function [yhat,status] = l1lsregr(x,y,lambda)
% L1 regularized linear spline regression smoothing, solved using
% interior-point quadratic programming.
%
% Input arguments:
%  x: covariate values
%  y: output variable values
%  lambda: positive L1 regularization parameter
%
% Output arguments:
%  yhat: linear spline smoothed solution
%  status: 'solved' or 'maxiter exceeded'

N = size(x,1);

% Input values must be in ascending order of covariates
[x,i] = sort(x);
y = y(i);

% Uniform grid spacing correction factor
hn = ((N-2)*N)/(max(x)-min(x))^2;

% Input value grid spacings
h0 = x(2:end-1)-x(1:end-2);
h1 = x(3:end)-x(2:end-1);

% 2nd derivative, non-uniform grid finite difference coefficients
a = 2./(h0.*(h0+h1));
b = -2./(h0.*h1);
c = 2./(h1.*(h0+h1));

% Construct 2nd derivative differencing matrix
A = spdiags(a,0,N-2,N-2);
B = spdiags(b,0,N-2,N-2);
C = spdiags(c,0,N-2,N-2);
O2 = zeros(N-2,1);
D = ([A O2 O2] + [O2 B O2] + [O2 O2 C])/hn;

% Perform L1 trend filtering with non-uniform differencing matrix
[yhat,status] = l1tf(y,lambda,D);
