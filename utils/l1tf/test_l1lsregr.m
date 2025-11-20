clear all;

% Generate underlying PWL signal
Ne = 7*365;
te = linspace(0,2,Ne);
sigma = 5.3;
Nb = floor(0.7*Ne);
xe = [linspace(15,50,Nb)'; linspace(50,25,Ne-Nb)'];

% Generate noisy signal
ue = xe + sigma*randn(Ne,1);

% Sample data from noisy signal
N = 150;
i = sort(randsample(Ne,N));
t = te(i)';
u = ue(i);

% Range of possible regularization constants
lambda = logspace(-2,6,50);

% CV train fraction
ftrncv = 0.5;

% Cross validation repeats
% Ncv = 100;
Ncv = 50;

% lmax = l1lsregrmax(t,u)/2

[uhat,lambdacv,Ecv] = l1lscv(t,u,lambda,ftrncv,Ncv);

%%
close all;
figure;
subplot(2,1,1);
hold on;
plot(te,ue,'-','Color',[0.8 0.8 0.8]);
plot(t,u,'k.');
plot(te,xe,'b-');
plot(t,uhat,'r.-');
subplot(2,1,2);
semilogx(lambda,Ecv,'k.-');
