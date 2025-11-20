clear all; close all;

%% Load sensor data
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_L1trend\trends2.mat')
trend_tremor_time_unmedicated_filled_above_threshold = trend_tremor_time_unmedicated_filled;
trend_tremor_time_unmedicated_filled_above_threshold(isnan(trend_median_tremor_power_unmedicated_filled))=NaN;
trend_tremor_time_medicated_filled_above_threshold = trend_tremor_time_medicated_filled;
trend_tremor_time_medicated_filled_above_threshold(isnan(trend_median_tremor_power_medicated_filled))=NaN;


%% Load IDs
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\IDs.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250410_matching_info\IDs_BaselineMedicated_matched.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_l1trend\IDs_BaselineUnmedicated_L1trend.mat')

IDs_BaselineMedicated_complete = IDs_BaselineMedicated(all(~isnan(trend_tremor_time_medicated_filled(1:length(IDs_BaselineMedicated),[2 51])),2));
IDs_BaselineMedicated_matched_complete = intersect(IDs_BaselineMedicated_complete,IDs_BaselineMedicated_matched);

IDs_BaselineMedicated_tremor = IDs_BaselineMedicated(all(~isnan(trend_median_tremor_power_medicated_filled(1:length(IDs_BaselineMedicated),[2 51])),2));
IDs_BaselineMedicated_matched_tremor = intersect(IDs_BaselineMedicated_tremor,IDs_BaselineMedicated_matched);
IDs_BaselineMedicated_matched_notremor = setxor(IDs_BaselineMedicated_tremor,IDs_BaselineMedicated_matched);

IDs_BaselineUnmedicated = [IDs_StartMedication; IDs_AllUnmedicated];
IDs_BaselineUnmedicated_complete = IDs_BaselineUnmedicated(all(~isnan(trend_tremor_time_unmedicated_filled(:,[2 51])),2));

IDs_BaselineUnmedicated_tremor = IDs_BaselineUnmedicated(all(~isnan(trend_median_tremor_power_unmedicated_filled(:,[2 51])),2));

%% Plot distribution over time (medicated group)
close all;

idx = contains([IDs_BaselineMedicated; IDs_StartMedication],IDs_BaselineMedicated);
measure = 100*trend_tremor_time_medicated_filled_above_threshold(idx,2:51);
% measure = trend_perc90_tremor_power_medicated_filled(idx,2:51);
N = sum(~isnan(measure))

figure(); hold on;
Xband = [0:2:98 98:-2:0];
plot(0:2:98,median(measure,'omitnan'),'k','LineWidth',2)
upperband = prctile(measure,75);
fill(Xband,[prctile(measure,25) upperband(end:-1:1)],'k','FaceAlpha',0.3,'EdgeColor','none')
upperband = prctile(measure,90);
fill(Xband,[prctile(measure,10) upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
xlabel('Weeks since baseline')
ylabel('Tremor time (% of inactive time)')
% ylabel('Modal tremor power (log values)')
ylim([0 100])
% ylim([0 4])
xlim([0 100])
% legend('median', '25th-75th percentiles','10th-90th percentiles')
title('Medicated group')

%% Plot delta progression (medicated group)
% close all;

idx = contains([IDs_BaselineMedicated; IDs_StartMedication],IDs_BaselineMedicated);
measure = logit(trend_tremor_time_medicated_filled_above_threshold(idx,3:51)) - logit(trend_tremor_time_medicated_filled_above_threshold(idx,2));
% measure = trend_perc90_tremor_power_medicated_filled(idx,3:51) - trend_perc90_tremor_power_medicated_filled(idx,2);
N = sum(~isnan(measure))

figure(); hold on;
Xband = [2:2:98 98:-2:2];
plot(2:2:98,median(measure,'omitnan'),'.-k')
upperband = prctile(measure,75);
fill(Xband,[prctile(measure,25) upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
xlabel('Weeks since baseline')
yline(0,'--')
% ylabel('\Delta 90th percentile of tremor power (log values)')
ylabel('\Delta tremor time (log odds-ratio)')
ylim([-1 1.5])
xlim([0 100])
legend('median', 'IQR')
title('Medicated group')

%% Plot distribution over time  (unmedicated group)

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\v3_survival\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(77,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100

weights = weights(contains(IDs_BaselineUnmedicated_L1trend,IDs_BaselineUnmedicated_L1trend),:);

idx = contains([IDs_StartMedication; IDs_AllUnmedicated],IDs_BaselineUnmedicated_L1trend);
measure = 100*trend_tremor_time_unmedicated_filled_above_threshold(idx,2:51);
% measure = trend_perc90_tremor_power_unmedicated_filled(idx,2:51);

n = size(measure,2);
medians = zeros(n, 1);
p25 = zeros(n, 1);
p75 = zeros(n, 1);
p10 = zeros(n, 1);
p90 = zeros(n, 1);

for i = 1:n

    outcome_t = measure(:,i);
    weights_t = weights(:,i);

    valid = ~isnan(outcome_t) & ~isnan(weights_t);
    outcome_t = outcome_t(valid);
    weights_t = weights_t(valid);
    
    % percentiles = weighted_percentiles(outcome_t, weights_t, [10, 25, 50, 75, 90]);
    percentiles = wprctile(outcome_t,[10, 25, 50, 75, 90],weights_t);
    
    p10(i) = percentiles(1);
    p25(i) = percentiles(2);
    medians(i) = percentiles(3);
    p75(i) = percentiles(4);
    p90(i) = percentiles(5);

end

N = sum(~isnan(measure))

figure(); hold on;
Xband = [0:2:98 98:-2:0];
plot(0:2:98,medians,'k','LineWidth',2)
upperband = p75;
fill(Xband,[p25; upperband(end:-1:1)],'k','FaceAlpha',0.3,'EdgeColor','none')
upperband = p90;
fill(Xband,[p10; upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
xlabel('Weeks since baseline')
ylabel('Tremor time (% of inactive time)')
% ylabel('90th percentile of tremor power (log values)')
ylim([0 100])
xlim([0 100])
legend('median', '25th-75th percentiles','10th-90th percentiles')
title('Unmedicated group')

%% Plot delta progression (unmedicated group)

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\v3_survival\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(77,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100

weights = weights(contains(IDs_BaselineUnmedicated_L1trend,IDs_BaselineUnmedicated_L1trend),:);

idx = contains([IDs_StartMedication; IDs_AllUnmedicated],IDs_BaselineUnmedicated_L1trend);
% measure = trend_modal_tremor_power_unmedicated_filled(idx,3:end) -  trend_modal_tremor_power_unmedicated_filled(idx,2);
measure = logit(trend_tremor_time_unmedicated_filled_above_threshold(idx,3:end)) - logit(trend_tremor_time_unmedicated_filled_above_threshold(idx,2));

n = size(measure,2)-2;
medians = zeros(n, 1);
p25 = zeros(n, 1);
p75 = zeros(n, 1);
p10 = zeros(n, 1);
p90 = zeros(n, 1);

for i = 1:n

    outcome_t = measure(:,i);
    weights_t = weights(:,i+2);

    valid = ~isnan(outcome_t) & ~isnan(weights_t);
    outcome_t = outcome_t(valid);
    weights_t = weights_t(valid);
    
    % percentiles = weighted_percentiles(outcome_t, weights_t, [10, 25, 50, 75, 90]);
    percentiles = wprctile(outcome_t,[10, 25, 50, 75, 90],weights_t);
    
    p10(i) = percentiles(1);
    p25(i) = percentiles(2);
    medians(i) = percentiles(3);
    p75(i) = percentiles(4);
    p90(i) = percentiles(5);

end

N = sum(~isnan(measure))

figure(); hold on;
Xband = [2:2:98 98:-2:2];
plot(2:2:98,medians,'.-k')
upperband = p75;
fill(Xband,[p25; upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
xlabel('Weeks since baseline')
yline(0,'--')
ylabel('\Delta tremor time (log odds-ratio)')
% ylabel('\Delta modal of tremor power (log values)')
ylim([-1 1.5])
xlim([0 100])
legend('median', 'IQR')
title('Unmedicated group')

%% Plot SRM over time (medicated group)
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))
SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(trend_tremor_time_medicated_filled));
% data = tremor_time_logit;
data = trend_modal_tremor_power_medicated_filled;
first_week_idx = 2;
last_week_idxs = 3:51;

idx = contains([IDs_BaselineMedicated; IDs_StartMedication],IDs_BaselineMedicated);

SRM_sensor_medicated = [];
CI_sensor_medicated = [];
N_sensor_medicated = [];

for m = 1:length(last_week_idxs)
    last_week_idx = last_week_idxs(m);
    
    SRM_sensor_medicated(m) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    CI_sensor_medicated(m,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));
    N_sensor_medicated(m) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx))));

end

week_vector = 0:2:104;
Xband = [2:2:98 98:-2:2];

figure(); hold on;
plot(week_vector(last_week_idxs)-2,SRM_sensor_medicated,'.-',Color='k')
upperband = CI_sensor_medicated(:,2)';
fill(Xband,[CI_sensor_medicated(:,1)' upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
yline(0,'--')
ylabel('SRM with 95%-CI')
ylim([-1 2])
xlabel('Weeks since baseline')
% title('Medicated group')

%% Plot SRM over time (unmedicated group)
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))

weighted_SRM_function = @(change, weights) ...
    (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan')) / ...
    sqrt(sum(weights .* (change - (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan'))).^2, 'omitnan') / sum(weights, 'omitnan'));

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\v3_survival\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(77,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100
weights = weights(contains(IDs_BaselineUnmedicated_L1trend,IDs_BaselineUnmedicated_L1trend),:);

idx = contains([IDs_StartMedication; IDs_AllUnmedicated],IDs_BaselineUnmedicated_L1trend);

tremor_time_logit = real(logit(trend_tremor_time_unmedicated_filled));

data = tremor_time_logit;

% data = trend_modal_tremor_power_unmedicated_filled;
first_week_idx = 2;
last_week_idxs = 3:51;

SRM_sensor_unmedicated = [];
CI_sensor_unmedicated = [];
N_sensor_unmedicated = [];

for m = 1:length(last_week_idxs)

    last_week_idx = last_week_idxs(m);

    change = data(idx,last_week_idx) - data(idx,first_week_idx);
    
    weighted_SRM_sensor_unmedicated(m) = weighted_SRM_function(change,weights(:,m+2));
    weighted_CI_sensor_unmedicated(m,1:2) = bootci(10000,{weighted_SRM_function,change,weights(:,m+2)});

    % SRM_sensor_unmedicated(m) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    % CI_sensor_unmedicated(m,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));

    N_sensor_unmedicated(m) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx))));
end

week_vector = 0:2:104;
Xband = [2:2:98 98:-2:2];

figure(); hold on;
plot(week_vector(last_week_idxs)-2,weighted_SRM_sensor_unmedicated,'.-',Color='k')
upperband = weighted_CI_sensor_unmedicated(:,2)';
fill(Xband,[weighted_CI_sensor_unmedicated(:,1)' upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
yline(0,'--')
ylabel('SRM with 95%-CI')
ylim([-1 2])
xlabel('Weeks since baseline')
% title('Unmedicated group')
