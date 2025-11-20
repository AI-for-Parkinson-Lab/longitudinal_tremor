clear all;
Ne = 15*365;
te = linspace(0,15,Ne);
sigma = 5.3;
Nb = floor(0.7*Ne);
xe = [linspace(15,50,Nb)'; linspace(50,25,Ne-Nb)'];

tt = [];
ee = [];
for i = 1:500
    ue = xe + sigma*randn(Ne,1);
    N = 75;
    % N = floor(0.1*Ne);
    i = sort(randsample(Ne,N));
    t = te(i)';
    u = ue(i);

    [us,status] = l1lsregr(t,u,200);
    ee = [ee; us - u];
    tt = [tt; t];
end

close all;
figure;
plot(tt,ee,'k.');


% Test for heteroscedasticity - blocks of time
tb = 0:(1/52):max(te);
tw = 1.0;
Nb = length(tb)-1;
lb = zeros(Nb,3);
for i = 1:Nb
    t0 = tb(i);
    t1 = t0 + tw;
    j = (tt >= t0) & (tt <= t1);
    [l,ci] = expfit(abs(ee(j)));
    lb(i,1) = l;
    lb(i,2) = ci(1);
    lb(i,3) = ci(2);
end

%%
% Errors seem to be skewed and heteroscedastic
figure;
hold on;
plot(tt,abs(ee),'.','Color',[0.8 0.8 0.8]);
plot(tb(1:end-1),lb,'k-');
axis tight;
