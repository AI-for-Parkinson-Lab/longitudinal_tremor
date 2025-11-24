%% Step 6b: Create plots of SRM over time
% Based on smoothed and raw weekly sensor-derived tremor measures of the two-year
% unmedicated (tremor) group

clear all; close all;

%% Load sensor data and IDs
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Trends_filled.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IDs_selected.mat'); 

%% Select two-year unmedicated participants (with tremor)
IDs_BaselineUnmedicated_complete = IDs_BaselineUnmedicated(all(~isnan(trend_tremor_time_unmedicated_filled(:,[2 51])),2));
IDs_BaselineUnmedicated_tremor = IDs_BaselineUnmedicated(all(~isnan(trend_modal_tremor_power_unmedicated_filled(:,[2 51])),2)); 

%% Plot SRM over time based on smoothed data

weighted_SRM_function = @(change, weights) ...
    (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan')) / ...
    sqrt(sum(weights .* (change - (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan'))).^2, 'omitnan') / sum(weights, 'omitnan'));

load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(78,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100, duplicate to increase availability of weights

idx = contains(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated_tremor); % Select full group or tremor group 

tremor_time_logit = real(logit(trend_tremor_time_unmedicated_filled));

% Select measure
% data = tremor_time_logit;
data = trend_perc90_tremor_power_unmedicated_filled;

first_week_idx = 2;
last_week_idxs = 3:51;

SRM_sensor_unmedicated = [];
CI_sensor_unmedicated = [];
N_sensor_unmedicated = [];

for m = 1:length(last_week_idxs)

    last_week_idx = last_week_idxs(m);
    change = data(idx,last_week_idx) - data(idx,first_week_idx);
    
    weighted_SRM_sensor_unmedicated(m) = weighted_SRM_function(change,weights(idx,m+2));
    weighted_CI_sensor_unmedicated(m,1:2) = bootci(10000,{weighted_SRM_function,change,weights(idx,m+2)});
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
title('Unmedicated tremor group')

%% Load raw weekly sensor data

addpath(genpath('utils\jsonlab'))
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Inclusion.mat');

filename = 'Tremor_aggregates.json';
weeks = 0:2:104;
tremor_time = [];
modal_tremor_power = [];
perc90_tremor_power = [];

for i = 1:length(IDs_BaselineUnmedicated)
    id = IDs_BaselineUnmedicated{i};
    for k = 1:length(weeks)
        week = weeks(k);
        if contains(Inclusion.Group(ismember(Inclusion.ID,id)),{'PPP'})
            file = ['C:\Users\z835211\Documents\Data\PPP\aggregated_output_191125\ppp\' num2str(week) '\' id '\' filename];
        else
            file = ['C:\Users\z835211\Documents\Data\DeNovo\aggregated_output_201125\denovo\' num2str(week) '\' id '\' filename];
        end
        if isfile(file)
            Tremor_aggregates = loadjson(file);
            if Tremor_aggregates.metadata.nr_valid_days >= 3
                tremor_time(i,k) = Tremor_aggregates.aggregated_tremor_measures.perc_windows_tremor/100; % store as proportion instead of percentage
                if tremor_time(i,k) >= 0.035
                    modal_tremor_power(i,k) = Tremor_aggregates.aggregated_tremor_measures.modal_tremor_power;
                    perc90_tremor_power(i,k) = Tremor_aggregates.aggregated_tremor_measures.x0x39_0p_tremor_power;
                else
                    modal_tremor_power(i,k) = NaN;
                    perc90_tremor_power(i,k) = NaN;
                end
            else
                tremor_time(i,k) = NaN;
                modal_tremor_power(i,k) = NaN;
                perc90_tremor_power(i,k) = NaN;
            end
        else
            tremor_time(i,k) = NaN;
            modal_tremor_power(i,k) = NaN;
            perc90_tremor_power(i,k) = NaN;
        end
    end
end

%% Plot SRM over time based on raw data

idx = contains(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated_tremor); % Select full group or tremor group 

tremor_time_logit = real(logit(tremor_time));

% Select measure
% data = tremor_time_logit;
data = modal_tremor_power;

first_week_idx = 2;
last_week_idxs = 3:51;

SRM_sensor_unmedicated = [];
CI_sensor_unmedicated = [];
N_sensor_unmedicated = [];

for m = 1:length(last_week_idxs)

    last_week_idx = last_week_idxs(m);
    change = data(idx,last_week_idx) - data(idx,first_week_idx);
    
    weighted_SRM_sensor_unmedicated(m) = weighted_SRM_function(change,weights(idx,m+2));
    weighted_CI_sensor_unmedicated(m,1:2) = bootci(10000,{weighted_SRM_function,change,weights(idx,m+2)});
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
title('Unmedicated tremor group')
