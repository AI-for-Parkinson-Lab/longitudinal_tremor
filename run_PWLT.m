%% Estimate the optimal lambda
lambdarng = logspace(-8,2,2000); % Range of possible regularization constants

N = length(FileListPD_updated);
Etst = NaN(N,length(lambdarng));

measure = tremor_amplitude_90th_rest; % weekly aggregated measure
week_vector = 0:2:104;

for i = 1:N
    tremor = measure(i,~isnan(measure(i,:)))';
    weeks = week_vector(~isnan(measure(i,:)))';
    if length(tremor)>3 % at least 3 datapoints are necessary
        Etst(i,:) = l1lscv(weeks, tremor, lambdarng);
    end
end

mean_Etst = mean(Etst,'omitnan');
[~,j] = min(mean_Etst);
lambdacv = lambdarng(j)

%% Fit the PWLT
uhat = NaN(N,length(week_vector));

for i = 1:N
        tremor = measure(i,~isnan(measure(i,:)))';
        weeks = week_vector(~isnan(measure(i,:)))';
        if length(tremor)>3
            uhat(i,~isnan(measure(i,:))) = l1lsregr(weeks,tremor,lambdacv);
        end
end