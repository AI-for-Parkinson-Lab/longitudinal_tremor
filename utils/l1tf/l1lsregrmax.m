function lambdamax = l1lsregrmax(x,y)
% Calculate the value of lambda so that if lambda >= lambdamax, the L1
% regression smoothing is minimized to the trivial least-squares line fit
% to the data. This can then be used to determine a useful range
% of values of lambda, for example.
%
% Usage:
% lambdamax = l1lsregrmax(x,y)
%
% Input arguments:
% - x          Covariates
% - y          Data to smooth
%
% Output arguments:
% - lambdamax  Maximum meaningful value of lambda
%
% (c) Max Little, 2012. If you use this code for your research, please cite:
% M.A. Little, Nick S. Jones (2012)
% "Generalized Methods and Solvers for Noise Removal from Piecewise
% Constant Signals: Part I - Background Theory"
% Proceedings of the Royal Society A (in press)

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

DDT = D*D';
Dy  = D*y;
lambdamax = max(abs(DDT\Dy));
