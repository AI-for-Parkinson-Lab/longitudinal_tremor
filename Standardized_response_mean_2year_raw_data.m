%% 
clear all; close all;

%% Load IDs
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_L1trend\trends2.mat')

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\IDs.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250410_matching_info\IDs_BaselineMedicated_matched.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_l1trend\IDs_BaselineUnmedicated_L1trend.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\Inclusion.mat')

for i = 1:length(IDs_BaselineUnmedicated_L1trend)
    id = IDs_BaselineUnmedicated_L1trend(i);
    StartWeek(i) = Inclusion.StartWeek(ismember(Inclusion.ID,id));
end

IDs_BaselineMedicated_tremor = IDs_BaselineMedicated(all(~isnan(trend_median_tremor_power_medicated_filled(1:length(IDs_BaselineMedicated),[2 51])),2));
IDs_BaselineMedicated_matched_tremor = intersect(IDs_BaselineMedicated_tremor,IDs_BaselineMedicated_matched,'stable');

IDs_BaselineUnmedicated = [IDs_StartMedication; IDs_AllUnmedicated];
IDs_BaselineUnmedicated_complete = IDs_BaselineUnmedicated(all(~isnan(trend_tremor_time_unmedicated_filled(:,[2 51])),2));
IDs_BaselineUnmedicated_tremor = IDs_BaselineUnmedicated(all(~isnan(trend_median_tremor_power_unmedicated_filled(:,[2 51])),2));

%% Load sensor data
addpath(genpath('C:\Users\z835211\Documents\GitHub\tsdf4matlab'))
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\jsonlab'))

filename = 'Tremor_aggregates.json';

week_vector = 0:2:104;

tremor_time = [];
median_tremor_power = [];
modal_tremor_power = [];
perc90_tremor_power = [];

IDs_inclusion = [IDs_BaselineMedicated; IDs_StartMedication; IDs_AllUnmedicated];

for i = 1:length(IDs_inclusion)

    id = IDs_inclusion{i};

    folder = strcat('C:\Users\z835211\OneDrive - Radboudumc\Weekly aggregates\',id);

    if isfolder(folder)

        [metadata, data] = load_tsdf_metadata_from_path(strcat(folder,'\Tremor_weekly_aggregates_meta.json'));

        for k=1:length(week_vector)
            if data{1,1}(k,3)>=3
                idx = find(data{1,1}(:,1)==week_vector(k));

                tremor_time(i,k) = data{1,2}(idx,5);

                if tremor_time(i,k) >= 0.035
                    median_tremor_power(i,k) = data{1,2}(idx,13);
                    modal_tremor_power(i,k) = data{1,2}(idx,14);
                    perc90_tremor_power(i,k) = data{1,2}(idx,15);
                else
                    median_tremor_power(i,k) = NaN;
                    modal_tremor_power(i,k) = NaN;
                    perc90_tremor_power(i,k) = NaN;
                end
            else
                tremor_time(i,k) = NaN;
                median_tremor_power(i,k) = NaN;
                modal_tremor_power(i,k) = NaN;
                perc90_tremor_power(i,k) = NaN;
            end
        end

    else

        for k = 1:length(week_vector)
            week = week_vector(k);
            file = ['C:\Users\z835211\Documents\Data\DeNovo\aggregated_output\' num2str(week) '\' id '\' filename];
            if isfile(file)
                Tremor_aggregates = loadjson(file);

                if Tremor_aggregates.metadata.nr_windows_total >= 3*9000

                    tremor_time(i,k) = Tremor_aggregates.aggregated_tremor_measures.perc_windows_tremor/100;
                    if tremor_time(i,k) >= 0.035
                        median_tremor_power(i,k) = Tremor_aggregates.aggregated_tremor_measures.median_tremor_power;
                        modal_tremor_power(i,k) = Tremor_aggregates.aggregated_tremor_measures.modal_tremor_power;
                        perc90_tremor_power(i,k) = Tremor_aggregates.aggregated_tremor_measures.x0x39_0p_tremor_power;
                    else
                        median_tremor_power(i,k) = NaN;
                        modal_tremor_power(i,k) = NaN;
                        perc90_tremor_power(i,k) = NaN;
                    end
                else
                    tremor_time(i,k) = NaN;
                    median_tremor_power(i,k) = NaN;
                    modal_tremor_power(i,k) = NaN;
                    perc90_tremor_power(i,k) = NaN;
                end

            else
                tremor_time(i,k) = NaN;
                median_tremor_power(i,k) = NaN;
                modal_tremor_power(i,k) = NaN;
                perc90_tremor_power(i,k) = NaN;
            end
        end
    end
end


%% Plot SRM over time (medicated group)
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))
SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(tremor_time));
data = tremor_time_logit;
% data = trend_perc90_tremor_power_medicated_filled;
first_week_idx = 2;
last_week_idxs = 3:51;

idx = contains(IDs_inclusion,IDs_BaselineMedicated_tremor);

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
xlabel('Weeks since week 2')
title('Medicated group')

%% Plot SRM over time (unmedicated group)
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))

weighted_SRM_function = @(change, weights) ...
    (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan')) / ...
    sqrt(sum(weights .* (change - (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan'))).^2, 'omitnan') / sum(weights, 'omitnan'));

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\v3_survival\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(77,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100
weights = weights(contains(IDs_BaselineUnmedicated_L1trend,IDs_BaselineUnmedicated_tremor),:);

idx = contains(IDs_inclusion,IDs_BaselineUnmedicated_tremor);

tremor_time_logit = real(logit(tremor_time));

% data = tremor_time_logit;

data = perc90_tremor_power;
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
xlabel('Weeks since week 2')
title('Unmedicated group')
