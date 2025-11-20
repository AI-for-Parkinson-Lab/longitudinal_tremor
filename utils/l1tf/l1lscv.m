% function [yhat,lambdacv,Ecv] = l1lscv(x,y,lambdarng,ftrn,Ncv)
function Ecv = l1lscv(x,y,lambdarng)%,ftrn,Ncv)

N = size(x,1);
% Ntrn = floor(N*ftrn);
Nl = length(lambdarng);

% Etrn = zeros(Nl,Ncv);
% Etst = zeros(Nl,Ncv);
Etst = zeros(Nl,N);

% Do train/test cross-validation
for n = 1:N%Ncv
    for j = 1:Nl

        % i = sort(randsample(N,Ntrn));
        % xtrn = x(i);
        % ytrn = y(i);
        % 
        % ii = setdiff(1:N,i);
        % xtst = x(ii);
        % ytst = y(ii);

        xtst = x(n);
        ytst = y(n);
        ii = setdiff(1:N,n);
        xtrn = x(ii);
        ytrn = y(ii);
  
        yhat = l1lsregr(xtrn,ytrn,lambdarng(j));
        % Etrn(j,n) = mean(abs(ytrn-yhat));

        yprd = l1lspred(xtst,xtrn,yhat);
        % Etst(j,n) = mean(abs(yprd-ytst));
        Etst(j,n) = abs(yprd-ytst);
    end
        % close all;
        % figure(1); hold on
        % scatter(xtrn,ytrn,'k')
        % plot(xtrn,yhat,'k')
        % scatter(xtst,ytst,'r','filled')
        % title(num2str(Etst(j,n)))
        % ylim([0 1])
end

% Find value of lambda at smallest mean test error
Ecv = mean(Etst,2);
% % Ecv = filtfilt(ones(15,1)/15,1,Ecv);
% % Ecv = medfilt1(Ecv,5);
% [v,i] = min(Ecv);
% lambdacv = lambdarng(i);
% 
% % Do smoothing on all data
% yhat = l1lsregr(x,y,lambdacv);
end