%% Step 6: Assess the sensitivity to change
% Compute one- and two-year standardized response means of sensor-derived
% tremor measures and clinical scores

clear all; close all;

%% Load sensor data and IDs
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Trends_filled.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IDs_selected.mat'); 
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Inclusion.mat');
start_week = Inclusion.StartWeek(ismember(Inclusion.ID,IDs_BaselineUnmedicated));

%% Select participants with tremor
IDs_BaselineMedicated_tremor = IDs_BaselineMedicated(all(~isnan(trend_modal_tremor_power_medicated_filled(1:length(IDs_BaselineMedicated),[2 51])),2)); % Change index for one- or two-year group
IDs_BaselineUnmedicated_tremor = IDs_BaselineUnmedicated(all(~isnan(trend_modal_tremor_power_unmedicated_filled(:,[2 51])),2)); % Change index for one- or two-year group

%% Load clinical data
Visit1DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit1_DeNovo.csv");
Visit2DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit2_DeNovo.csv");
Visit3DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit3_DeNovo.csv");

Visit1PPP = readtable("C:\Users\z835211\Documents\Data\PPP\csv_files\General_visit1.csv");
Visit2PPP = readtable("C:\Users\z835211\Documents\Data\PPP\csv_files\General_visit2.csv");
Visit3PPP = readtable("C:\Users\z835211\Documents\Data\PPP\csv_files\General_visit3.csv");

% Convert IDs of PPP to same format
ids = char(Visit1PPP.id);
new_ids = ids(:,1:16);
Visit1PPP.id  = cellstr(strcat(repmat('POMU',length(new_ids),1),new_ids));
ids = char(Visit2PPP.id);
new_ids = ids(:,1:16);
Visit2PPP.id  = cellstr(strcat(repmat('POMU',length(new_ids),1),new_ids));
ids = char(Visit3PPP.id);
new_ids = ids(:,1:16);
Visit3PPP.id  = cellstr(strcat(repmat('POMU',length(new_ids),1),new_ids));

%% Extract UPDRS scores baseline medicated group

UPDRS_317OFF_1 = [];
UPDRS_318OFF_1 = [];
UPDRS_317ON_1 = [];
UPDRS_318ON_1 = [];
UPDRS_210_1 = [];
UPDRS_317OFF_2 = [];
UPDRS_318OFF_2 = [];
UPDRS_317ON_2 = [];
UPDRS_318ON_2 = [];
UPDRS_210_2 = [];
UPDRS_317OFF_3 = [];
UPDRS_318OFF_3 = [];
UPDRS_317ON_3 = [];
UPDRS_318ON_3 = [];
UPDRS_210_3 = [];

load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\visit_week_numbers.mat') % load visit week numbers

for i = 1:length(IDs_BaselineMedicated)
    id = IDs_BaselineMedicated{i};

    UPDRS_317OFF_1(i) = Visit1PPP.Up3OfRAmpArmYesDev(contains(Visit1PPP.id,id));
    UPDRS_318OFF_1(i) = Visit1PPP.Up3OfConstan(contains(Visit1PPP.id,id));
    UPDRS_317ON_1(i) = Visit1PPP.Up3OnRAmpArmYesDev(contains(Visit1PPP.id,id));
    UPDRS_318ON_1(i) = Visit1PPP.Up3OnConstan(contains(Visit1PPP.id,id));
    UPDRS_210_1(i) = Visit1PPP.Updrs2It23(contains(Visit1PPP.id,id));

    if ~isempty(find(contains(Visit2PPP.id,id)))
        visit2_week_number = cell2mat(visit_week_numbers.visit2(contains(visit_week_numbers.ID,id)));
        if visit2_week_number <= 60
            UPDRS_317OFF_2(i) = Visit2PPP.Up3OfRAmpArmYesDev(contains(Visit2PPP.id,id));
            UPDRS_318OFF_2(i) = Visit2PPP.Up3OfConstan(contains(Visit2PPP.id,id));
            UPDRS_317ON_2(i) = Visit2PPP.Up3OnRAmpArmYesDev(contains(Visit2PPP.id,id));
            UPDRS_318ON_2(i) = Visit2PPP.Up3OnConstan(contains(Visit2PPP.id,id));
            UPDRS_210_2(i) = Visit2PPP.Updrs2It23(contains(Visit2PPP.id,id));
        else
            UPDRS_317OFF_2(i) = NaN;
            UPDRS_318OFF_2(i) = NaN;
            UPDRS_317ON_2(i) = NaN;
            UPDRS_318ON_2(i) = NaN;
            UPDRS_210_2(i) = NaN;
        end
    else
        UPDRS_317OFF_2(i) = NaN;
        UPDRS_318OFF_2(i) = NaN;
        UPDRS_317ON_2(i) = NaN;
        UPDRS_318ON_2(i) = NaN;
        UPDRS_210_2(i) = NaN;
    end

    if ~isempty(find(contains(Visit3PPP.id,id)))
        visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));
        if visit3_week_number <= 120
            UPDRS_317OFF_3(i) = Visit3PPP.Up3OfRAmpArmYesDev(contains(Visit3PPP.id,id));
            UPDRS_318OFF_3(i) = Visit3PPP.Up3OfConstan(contains(Visit3PPP.id,id));
            UPDRS_317ON_3(i) = Visit3PPP.Up3OnRAmpArmYesDev(contains(Visit3PPP.id,id));
            UPDRS_318ON_3(i) = Visit3PPP.Up3OnConstan(contains(Visit3PPP.id,id));
            UPDRS_210_3(i) = Visit3PPP.Updrs2It23(contains(Visit3PPP.id,id));
        else
            UPDRS_317OFF_3(i) = NaN;
            UPDRS_318OFF_3(i) = NaN;
            UPDRS_317ON_3(i) = NaN;
            UPDRS_318ON_3(i) = NaN;
            UPDRS_210_3(i) = NaN;
        end
    else
        UPDRS_317OFF_3(i) = NaN;
        UPDRS_318OFF_3(i) = NaN;
        UPDRS_317ON_3(i) = NaN;
        UPDRS_318ON_3(i) = NaN;
        UPDRS_210_3(i) = NaN;
    end
end

%% Determine SRM of UPDRS scores
SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');
UPDRS_names = {'UPDRS_317OFF','UPDRS_318OFF','UPDRS_317ON','UPDRS_318ON','UPDRS_210'};

SRM_UPDRS_medicated = [];
CI_UPDRS_medicated = [];
N_UPDRS_medicated = [];

for k = 1:length(UPDRS_names)
    UPDRS_score3 = eval([UPDRS_names{k} '_2']); % score 2 or 3
    UPDRS_score1 = eval([UPDRS_names{k} '_1']);

    SRM_UPDRS_medicated(k) = mean(UPDRS_score3-UPDRS_score1,'omitnan')/std(UPDRS_score3-UPDRS_score1,'omitnan');
    CI_UPDRS_medicated(k,1:2) = bootci(10000,SRM_function,UPDRS_score3 - UPDRS_score1);

    N_UPDRS_medicated(k) = length(find(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1)));
end

%% Determine SRM of sensor data
addpath(genpath('utils'))
tremor_time_logit = real(logit(trend_tremor_time_medicated_filled));

% Select measures (tremor time for full group, all measures for tremor group)
% sensor_names = {'tremor_time_logit','trend_modal_tremor_power_medicated_filled','trend_perc90_tremor_power_medicated_filled'};
sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 26; % Modify for one- or two-year follow-up

idx = contains([IDs_BaselineMedicated; IDs_BaselineUnmedicated],IDs_BaselineMedicated);

SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

SRM_sensor_medicated = [];
CI_sensor_medicated = [];
N_sensor_medicated = [];

for k = 1:length(sensor_names)    

    data = eval(sensor_names{k});
    
    SRM_sensor_medicated(k) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    CI_sensor_medicated(k,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));

    N_sensor_medicated(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx))));
end

%% Extract UPDRS scores unmedicated group (no ON assessment)
UPDRS_317OFF_1 = [];
UPDRS_318OFF_1 = [];
UPDRS_210_1 = [];
UPDRS_317OFF_2 = [];
UPDRS_318OFF_2 = [];
UPDRS_210_2 = [];
UPDRS_317OFF_3 = [];
UPDRS_318OFF_3 = [];
UPDRS_210_3 = [];

for i = 1:length(IDs_BaselineUnmedicated)
    id = IDs_BaselineUnmedicated{i};

    if ~isempty(find(contains(Visit1PPP.id,id)))

        UPDRS_317OFF_1(i) = Visit1PPP.Up3OfRAmpArmYesDev(contains(Visit1PPP.id,id));
        UPDRS_318OFF_1(i) = Visit1PPP.Up3OfConstan(contains(Visit1PPP.id,id));
        UPDRS_210_1(i) = Visit1PPP.Updrs2It23(contains(Visit1PPP.id,id));

        if ~isempty(find(contains(Visit2PPP.id,id)))
            visit2_week_number = cell2mat(visit_week_numbers.visit2(contains(visit_week_numbers.ID,id)));
            if visit2_week_number <= 60
                UPDRS_317OFF_2(i) = Visit2PPP.Up3OfRAmpArmYesDev(contains(Visit2PPP.id,id));
                UPDRS_318OFF_2(i) = Visit2PPP.Up3OfConstan(contains(Visit2PPP.id,id));
                UPDRS_210_2(i) = Visit2PPP.Updrs2It23(contains(Visit2PPP.id,id));
            else
                UPDRS_317OFF_2(i) = NaN;
                UPDRS_318OFF_2(i) = NaN;
                UPDRS_210_2(i) = NaN;
            end
        else
            UPDRS_317OFF_2(i) = NaN;
            UPDRS_318OFF_2(i) = NaN;
            UPDRS_210_2(i) = NaN;
        end

        if ~isempty(find(contains(Visit3PPP.id,id)))
            visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));
            if visit3_week_number <= 120
                UPDRS_317OFF_3(i) = Visit3PPP.Up3OfRAmpArmYesDev(contains(Visit3PPP.id,id));
                UPDRS_318OFF_3(i) = Visit3PPP.Up3OfConstan(contains(Visit3PPP.id,id));
                UPDRS_210_3(i) = Visit3PPP.Updrs2It23(contains(Visit3PPP.id,id));
            else
                UPDRS_317OFF_3(i) = NaN;
                UPDRS_318OFF_3(i) = NaN;
                UPDRS_210_3(i) = NaN;
            end
        else
            UPDRS_317OFF_3(i) = NaN;
            UPDRS_318OFF_3(i) = NaN;
            UPDRS_210_3(i) = NaN;
        end

    else
        UPDRS_317OFF_1(i) = Visit1DeNovo.Up3OfRAmpArmYesDev(contains(Visit1DeNovo.id,id));
        UPDRS_318OFF_1(i) = Visit1DeNovo.Up3OfConstan(contains(Visit1DeNovo.id,id));
        UPDRS_210_1(i) = Visit1DeNovo.Updrs2It23(contains(Visit1DeNovo.id,id));

        if contains(id,'POMU600C11F136E6FB4D')

            UPDRS_317OFF_2(i) = NaN;
            UPDRS_318OFF_2(i) = NaN;
            UPDRS_210_2(i) = NaN;

            UPDRS_317OFF_3(i) = Visit2DeNovo.Up3OfRAmpArmYesDev(contains(Visit2DeNovo.id,id));
            UPDRS_318OFF_3(i) = Visit2DeNovo.Up3OfConstan(contains(Visit2DeNovo.id,id));
            UPDRS_210_3(i) = Visit2DeNovo.Updrs2It23(contains(Visit2DeNovo.id,id));

        else
            if ~isempty(find(contains(Visit2DeNovo.id,id)))
                visit2_week_number = cell2mat(visit_week_numbers.visit2(contains(visit_week_numbers.ID,id)));
                if visit2_week_number <= 60
                    UPDRS_317OFF_2(i) = Visit2DeNovo.Up3OfRAmpArmYesDev(contains(Visit2DeNovo.id,id));
                    UPDRS_318OFF_2(i) = Visit2DeNovo.Up3OfConstan(contains(Visit2DeNovo.id,id));
                    UPDRS_210_2(i) = Visit2DeNovo.Updrs2It23(contains(Visit2DeNovo.id,id));
                else
                    UPDRS_317OFF_2(i) = NaN;
                    UPDRS_318OFF_2(i) = NaN;
                    UPDRS_210_2(i) = NaN;
                end
            else
                UPDRS_317OFF_2(i) = NaN;
                UPDRS_318OFF_2(i) = NaN;
                UPDRS_210_2(i) = NaN;
            end
            if ~isempty(find(contains(Visit3DeNovo.id,id)))
                visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));
                if visit3_week_number <= 120
                    UPDRS_317OFF_3(i) = Visit3DeNovo.Up3OfRAmpArmYesDev(contains(Visit3DeNovo.id,id));
                    UPDRS_318OFF_3(i) = Visit3DeNovo.Up3OfConstan(contains(Visit3DeNovo.id,id));
                    UPDRS_210_3(i) = Visit3DeNovo.Updrs2It23(contains(Visit3DeNovo.id,id));
                else
                    UPDRS_317OFF_3(i) = NaN;
                    UPDRS_318OFF_3(i) = NaN;
                    UPDRS_210_3(i) = NaN;
                end
            else
                UPDRS_317OFF_3(i) = NaN;
                UPDRS_318OFF_3(i) = NaN;
                UPDRS_210_3(i) = NaN;
            end
        end
    end
end

%% Determine SRM of UPDRS scores

load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IPCW.mat')

weighted_SRM_function = @(change, weights) ...
    (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan')) / ...
    sqrt(sum(weights .* (change - (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan'))).^2, 'omitnan') / sum(weights, 'omitnan'));

UPDRS_names = {'UPDRS_317OFF','UPDRS_318OFF','UPDRS_210'};

weighted_SRM_UPDRS_unmedicated = [];
weighted_CI_UPDRS_unmedicated = [];
N_UPDRS_unmedicated = [];

weights = IPCW(:,2:end).Variables;
weights = [ones(78,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100, duplicate to increase availability of weights

idx_2year = contains(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated); % Select full group or tremor group
last_week_idx = 26;
weights_2year = weights(idx_2year,last_week_idx);

for k = 1:length(UPDRS_names)

    UPDRS_score3 = eval([UPDRS_names{k} '_2']);
    UPDRS_score1 = eval([UPDRS_names{k} '_1']);

    change = (UPDRS_score3-UPDRS_score1)';

    weighted_SRM_UPDRS_unmedicated(k) = weighted_SRM_function(change,weights_2year);
    weighted_CI_UPDRS_unmedicated(k,1:2)= bootci(10000,{weighted_SRM_function,change,weights_2year});

    N_UPDRS_unmedicated(k) = length(find(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1) & ~isnan(weights_2year')));
end

%% Determine SRM of sensor data

SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(trend_tremor_time_unmedicated_filled));

% Select measures (tremor time for full group, all measures for tremor group)
% sensor_names = {'tremor_time_logit','trend_modal_tremor_power_unmedicated_filled','trend_perc90_tremor_power_unmedicated_filled'};
sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 26;

idx = contains(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated);

weighted_SRM_sensor_unmedicated = [];
weighted_CI_sensor_unmedicated = [];
N_sensor_unmedicated = [];

for k = 1:length(sensor_names)

    data = eval(sensor_names{k});

    change = data(idx,last_week_idx) - data(idx,first_week_idx);

    weights_boostrap = weights_2year;

    weighted_SRM_sensor_unmedicated(k) = weighted_SRM_function(change,weights_boostrap);
    weighted_CI_sensor_unmedicated(k,1:2) = bootci(10000,{weighted_SRM_function,change,weights_boostrap});

    N_sensor_unmedicated(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx)) & ~isnan(weights_boostrap)));
end

%% Create forest plot
% Names = {'rest tremor severity in OFF', 'rest tremor constancy in OFF', 'rest tremor severity in ON', 'rest tremor constancy in ON', 'patient-reported tremor','tremor time', 'modal tremor power', '90th percentile of tremor power',...
%     'rest tremor severity in OFF', 'rest tremor constancy in OFF', 'patient-reported tremor', 'tremor time', 'modal tremor power', '90th percentile of tremor power'}; % Study names]
Names = {'rest tremor severity in OFF', 'rest tremor constancy in OFF', 'rest tremor severity in ON', 'rest tremor constancy in ON', 'patient-reported tremor','tremor time',...
    'rest tremor severity in OFF', 'rest tremor constancy in OFF', 'patient-reported tremor', 'tremor time'}; % Study names]
Names = strcat(Names,' (n =',{' '},string([N_UPDRS_medicated N_sensor_medicated N_UPDRS_unmedicated N_sensor_unmedicated]),')');

SRM = [SRM_UPDRS_medicated SRM_sensor_medicated weighted_SRM_UPDRS_unmedicated weighted_SRM_sensor_unmedicated]; % Standardized Response Means
CI = [CI_UPDRS_medicated; CI_sensor_medicated; weighted_CI_UPDRS_unmedicated; weighted_CI_sensor_unmedicated]; % Lower bounds of 95% confidence intervals

% Plotting
figure;
hold on;
C = colororder('glow12');

% Number of studies
n = length(SRM);

% Plot the confidence intervals
for i = 1:n
    if i == 6 || i == 10
    % if i == 6 || i == 7 || i == 8 || i == 12 || i == 13 || i == 14
        line([CI(i,1), CI(i,2)], [n-i+1, n-i+1], 'Color', C(2,:), 'LineWidth', 1.5); % CI line
        plot(SRM(i), n-i+1, 'ko', 'MarkerFaceColor', C(2,:), 'MarkerEdgeColor', C(2,:), 'MarkerSize', 6); % SRM point
    else
        line([CI(i,1), CI(i,2)], [n-i+1, n-i+1], 'Color', 'k', 'LineWidth', 1.5); % CI line
        plot(SRM(i), n-i+1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % SRM point
    end
end

% Add labels
yticks(1:n)
yticklabels(Names(end:-1:1))
xlabel('SRM with 95%-CI');
title('1-year change')

xlim([-0.5 1.2]);
ylim([0 n+1])

% Reference line at zero effect
line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

hold off;

%%  Select participants 
% IDs_BaselineUnmedicated_statistics = IDs_AllUnmedicated(~isnan(trend_tremor_time_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_AllUnmedicated),2)) ...
%     & ~isnan(trend_tremor_time_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_AllUnmedicated),51)));
% % IDs_BaselineUnmedicated_statistics = IDs_1year_unmedicated(~isnan(trend_tremor_time_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_1year_unmedicated),2)) ...
% %     & ~isnan(trend_tremor_time_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_1year_unmedicated),26)));
% IDs_BaselineUnmedicated_statistics_tremor = intersect(IDs_BaselineUnmedicated_statistics,IDs_BaselineUnmedicated_tremor,'stable');
% IDs_BaselineMedicated_statistics = IDs_BaselineMedicated(~isnan(trend_tremor_time_medicated_filled(1:length(IDs_BaselineMedicated),2)) ...
%     & ~isnan(trend_tremor_time_medicated_filled(1:length(IDs_BaselineMedicated),26)));
% IDs_BaselineMedicated_statistics_tremor =  intersect(IDs_BaselineMedicated_statistics,IDs_BaselineMedicated_tremor,'stable');

%% Extract UPDRS scores baseline medicated group

UPDRS_317OFF_1 = [];
UPDRS_318OFF_1 = [];
UPDRS_317ON_1 = [];
UPDRS_318ON_1 = [];
UPDRS_210_1 = [];
UPDRS_317OFF_2 = [];
UPDRS_318OFF_2 = [];
UPDRS_317ON_2 = [];
UPDRS_318ON_2 = [];
UPDRS_210_2 = [];
UPDRS_317OFF_3 = [];
UPDRS_318OFF_3 = [];
UPDRS_317ON_3 = [];
UPDRS_318ON_3 = [];
UPDRS_210_3 = [];

load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\visit_week_numbers.mat') % load visit week numbers


for i = 1:length(IDs_BaselineMedicated_tremor)
    id = IDs_BaselineMedicated_tremor{i};

    UPDRS_317OFF_1(i) = Visit1PPP.Up3OfRAmpArmYesDev(contains(Visit1PPP.id,id));
    UPDRS_318OFF_1(i) = Visit1PPP.Up3OfConstan(contains(Visit1PPP.id,id));
    UPDRS_317ON_1(i) = Visit1PPP.Up3OnRAmpArmYesDev(contains(Visit1PPP.id,id));
    UPDRS_318ON_1(i) = Visit1PPP.Up3OnConstan(contains(Visit1PPP.id,id));
    UPDRS_210_1(i) = Visit1PPP.Updrs2It23(contains(Visit1PPP.id,id));

    if ~isempty(find(contains(Visit2PPP.id,id)))
        visit2_week_number = cell2mat(visit_week_numbers.visit2(contains(visit_week_numbers.ID,id)));
        if visit2_week_number <= 60
            UPDRS_317OFF_2(i) = Visit2PPP.Up3OfRAmpArmYesDev(contains(Visit2PPP.id,id));
            UPDRS_318OFF_2(i) = Visit2PPP.Up3OfConstan(contains(Visit2PPP.id,id));
            UPDRS_317ON_2(i) = Visit2PPP.Up3OnRAmpArmYesDev(contains(Visit2PPP.id,id));
            UPDRS_318ON_2(i) = Visit2PPP.Up3OnConstan(contains(Visit2PPP.id,id));
            UPDRS_210_2(i) = Visit2PPP.Updrs2It23(contains(Visit2PPP.id,id));
        else
            UPDRS_317OFF_2(i) = NaN;
            UPDRS_318OFF_2(i) = NaN;
            UPDRS_317ON_2(i) = NaN;
            UPDRS_318ON_2(i) = NaN;
            UPDRS_210_2(i) = NaN;
        end
    else
        UPDRS_317OFF_2(i) = NaN;
        UPDRS_318OFF_2(i) = NaN;
        UPDRS_317ON_2(i) = NaN;
        UPDRS_318ON_2(i) = NaN;
        UPDRS_210_2(i) = NaN;
    end

    if ~isempty(find(contains(Visit3PPP.id,id)))
        visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));
        if visit3_week_number <= 120
            UPDRS_317OFF_3(i) = Visit3PPP.Up3OfRAmpArmYesDev(contains(Visit3PPP.id,id));
            UPDRS_318OFF_3(i) = Visit3PPP.Up3OfConstan(contains(Visit3PPP.id,id));
            UPDRS_317ON_3(i) = Visit3PPP.Up3OnRAmpArmYesDev(contains(Visit3PPP.id,id));
            UPDRS_318ON_3(i) = Visit3PPP.Up3OnConstan(contains(Visit3PPP.id,id));
            UPDRS_210_3(i) = Visit3PPP.Updrs2It23(contains(Visit3PPP.id,id));
        else
            UPDRS_317OFF_3(i) = NaN;
            UPDRS_318OFF_3(i) = NaN;
            UPDRS_317ON_3(i) = NaN;
            UPDRS_318ON_3(i) = NaN;
            UPDRS_210_3(i) = NaN;
        end
    else
        UPDRS_317OFF_3(i) = NaN;
        UPDRS_318OFF_3(i) = NaN;
        UPDRS_317ON_3(i) = NaN;
        UPDRS_318ON_3(i) = NaN;
        UPDRS_210_3(i) = NaN;
    end
end

%% Determine SRM of UPDRS scores
SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');
UPDRS_names = {'UPDRS_317OFF','UPDRS_318OFF','UPDRS_317ON','UPDRS_318ON','UPDRS_210'};

SRM_UPDRS_medicated = [];
CI_UPDRS_medicated = [];
N_UPDRS_medicated = [];

bootstat_UPDRS_medicated = [];
Bootstrap_IDs = cell(1, length(UPDRS_names));  

for k = 1:length(UPDRS_names)
    UPDRS_score3 = eval([UPDRS_names{k} '_2']); % score 2 or 3
    UPDRS_score1 = eval([UPDRS_names{k} '_1']);

    SRM_UPDRS_medicated(k) = mean(UPDRS_score3-UPDRS_score1,'omitnan')/std(UPDRS_score3-UPDRS_score1,'omitnan');
    [CI_UPDRS_medicated(k,1:2), bootstat_UPDRS_medicated(:,k)] = bootci(10000,SRM_function,UPDRS_score3 - UPDRS_score1);

    N_UPDRS_medicated(k) = length(find(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1)));
    % Bootstrap_IDs{k} = IDs_BaselineMedicated_statistics(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1)); 

end

%% Determine SRM of sensor data
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))
SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(trend_tremor_time_medicated_filled));
sensor_names = {'tremor_time_logit','trend_modal_tremor_power_medicated_filled','trend_perc90_tremor_power_medicated_filled'};
% sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 26; % Modify for one- or two-year follow-up

idx = contains([IDs_BaselineMedicated; IDs_StartMedication],IDs_AllUnmedicated);

SRM_sensor_medicated = [];
CI_sensor_medicated = [];
N_sensor_medicated = [];
bootstat_sensor_medicated = [];

for k = 1:length(sensor_names)
% for k = 1:length(UPDRS_names)
    
    % idx = contains([IDs_BaselineMedicated; IDs_StartMedication],Bootstrap_IDs{k});

    data = eval(sensor_names{1});
    
    SRM_sensor_medicated(k) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    [CI_sensor_medicated(k,1:2) bootstat_sensor_medicated(:,k)] = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));

    % SE_mean_diff_data_PPP = std(data(idx,last_week_idx) - data(idx,1),'omitnan')/sqrt(length(find(~isnan(data(idx,last_week_idx))==1 & ~isnan(data(idx,1))==1)));
    % CI_sensor_data_low_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')-1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');
    % CI_sensor_data_high_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')+1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');

    N_sensor_medicated(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx))));
end

%% Extract UPDRS scores unmedicated group (no ON assessment) and disease duration
UPDRS_317OFF_1 = [];
UPDRS_318OFF_1 = [];
UPDRS_210_1 = [];
UPDRS_317OFF_2 = [];
UPDRS_318OFF_2 = [];
UPDRS_210_2 = [];
UPDRS_317OFF_3 = [];
UPDRS_318OFF_3 = [];
UPDRS_210_3 = [];

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\visit_week_numbers.mat')

for i = 1:length(IDs_BaselineUnmedicated_statistics)
    id = IDs_BaselineUnmedicated_statistics{i};

    if ~isempty(find(contains(Visit1PPP.id,id)))

        UPDRS_317OFF_1(i) = Visit1PPP.Up3OfRAmpArmYesDev(contains(Visit1PPP.id,id));
        UPDRS_318OFF_1(i) = Visit1PPP.Up3OfConstan(contains(Visit1PPP.id,id));
        UPDRS_210_1(i) = Visit1PPP.Updrs2It23(contains(Visit1PPP.id,id));

        if ~isempty(find(contains(Visit2PPP.id,id)))
            visit2_week_number = cell2mat(visit_week_numbers.visit2(contains(visit_week_numbers.ID,id)));
            if visit2_week_number <= 60
                UPDRS_317OFF_2(i) = Visit2PPP.Up3OfRAmpArmYesDev(contains(Visit2PPP.id,id));
                UPDRS_318OFF_2(i) = Visit2PPP.Up3OfConstan(contains(Visit2PPP.id,id));
                UPDRS_210_2(i) = Visit2PPP.Updrs2It23(contains(Visit2PPP.id,id));
            else
                UPDRS_317OFF_2(i) = NaN;
                UPDRS_318OFF_2(i) = NaN;
                UPDRS_210_2(i) = NaN;
            end
        else
            UPDRS_317OFF_2(i) = NaN;
            UPDRS_318OFF_2(i) = NaN;
            UPDRS_210_2(i) = NaN;
        end

        if ~isempty(find(contains(Visit3PPP.id,id)))
            visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));
            if visit3_week_number <= 120
                UPDRS_317OFF_3(i) = Visit3PPP.Up3OfRAmpArmYesDev(contains(Visit3PPP.id,id));
                UPDRS_318OFF_3(i) = Visit3PPP.Up3OfConstan(contains(Visit3PPP.id,id));
                UPDRS_210_3(i) = Visit3PPP.Updrs2It23(contains(Visit3PPP.id,id));
            else
                UPDRS_317OFF_3(i) = NaN;
                UPDRS_318OFF_3(i) = NaN;
                UPDRS_210_3(i) = NaN;
            end
        else
            UPDRS_317OFF_3(i) = NaN;
            UPDRS_318OFF_3(i) = NaN;
            UPDRS_210_3(i) = NaN;
        end

    else
        UPDRS_317OFF_1(i) = Visit1DeNovo.Up3OfRAmpArmYesDev(contains(Visit1DeNovo.id,id));
        UPDRS_318OFF_1(i) = Visit1DeNovo.Up3OfConstan(contains(Visit1DeNovo.id,id));
        UPDRS_210_1(i) = Visit1DeNovo.Updrs2It23(contains(Visit1DeNovo.id,id));

        if contains(id,'POMU600C11F136E6FB4D')

            UPDRS_317OFF_2(i) = NaN;
            UPDRS_318OFF_2(i) = NaN;
            UPDRS_210_2(i) = NaN;

            UPDRS_317OFF_3(i) = Visit2DeNovo.Up3OfRAmpArmYesDev(contains(Visit2DeNovo.id,id));
            UPDRS_318OFF_3(i) = Visit2DeNovo.Up3OfConstan(contains(Visit2DeNovo.id,id));
            UPDRS_210_3(i) = Visit2DeNovo.Updrs2It23(contains(Visit2DeNovo.id,id));

        else
            if ~isempty(find(contains(Visit2DeNovo.id,id)))
                visit2_week_number = cell2mat(visit_week_numbers.visit2(contains(visit_week_numbers.ID,id)));
                if visit2_week_number <= 60
                    UPDRS_317OFF_2(i) = Visit2DeNovo.Up3OfRAmpArmYesDev(contains(Visit2DeNovo.id,id));
                    UPDRS_318OFF_2(i) = Visit2DeNovo.Up3OfConstan(contains(Visit2DeNovo.id,id));
                    UPDRS_210_2(i) = Visit2DeNovo.Updrs2It23(contains(Visit2DeNovo.id,id));
                else
                    UPDRS_317OFF_2(i) = NaN;
                    UPDRS_318OFF_2(i) = NaN;
                    UPDRS_210_2(i) = NaN;
                end
            else
                UPDRS_317OFF_2(i) = NaN;
                UPDRS_318OFF_2(i) = NaN;
                UPDRS_210_2(i) = NaN;
            end
            if ~isempty(find(contains(Visit3DeNovo.id,id)))
                visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));
                if visit3_week_number <= 120
                    UPDRS_317OFF_3(i) = Visit3DeNovo.Up3OfRAmpArmYesDev(contains(Visit3DeNovo.id,id));
                    UPDRS_318OFF_3(i) = Visit3DeNovo.Up3OfConstan(contains(Visit3DeNovo.id,id));
                    UPDRS_210_3(i) = Visit3DeNovo.Updrs2It23(contains(Visit3DeNovo.id,id));
                else
                    UPDRS_317OFF_3(i) = NaN;
                    UPDRS_318OFF_3(i) = NaN;
                    UPDRS_210_3(i) = NaN;
                end
            else
                UPDRS_317OFF_3(i) = NaN;
                UPDRS_318OFF_3(i) = NaN;
                UPDRS_210_3(i) = NaN;
            end
        end
    end
end

%% Determine SRM of UPDRS scores

weighted_SRM_function = @(change, weights) ...
    (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan')) / ...
    sqrt(sum(weights .* (change - (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan'))).^2, 'omitnan') / sum(weights, 'omitnan'));

UPDRS_names = {'UPDRS_317OFF','UPDRS_318OFF','UPDRS_210'};

SRM_UPDRS_unmedicated = [];
CI_UPDRS_unmedicated = [];
N_UPDRS_unmedicated = [];
weighted_SRM_UPDRS_unmedicated = [];
weighted_CI_UPDRS_unmedicated = [];
bootstat_UPDRS_unmedicated = [];
Bootstrap_IDs = cell(1, length(UPDRS_names));  

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\v3_survival\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(77,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100

idx_2year = contains(IDs_BaselineUnmedicated_L1trend,IDs_BaselineUnmedicated_statistics);
last_week_idx = 51;
weights_2year = weights(idx_2year,last_week_idx);

for k = 1:length(UPDRS_names)

    UPDRS_score3 = eval([UPDRS_names{k} '_3']);
    UPDRS_score1 = eval([UPDRS_names{k} '_1']);

    change = (UPDRS_score3-UPDRS_score1)';

    weighted_SRM_UPDRS_unmedicated(k) = weighted_SRM_function(change,weights_2year);
    [weighted_CI_UPDRS_unmedicated(k,1:2), bootstat_UPDRS_unmedicated(:,k)] = bootci(10000,{weighted_SRM_function,change,weights_2year});

    % SRM_UPDRS_unmedicated(k) = mean(UPDRS_score3-UPDRS_score1,'omitnan')/std(UPDRS_score3-UPDRS_score1,'omitnan');
    % CI_UPDRS_unmedicated(k,1:2) = bootci(10000,SRM_function,UPDRS_score3 - UPDRS_score1);

    N_UPDRS_unmedicated(k) = length(find(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1) & ~isnan(weights_2year')));
    Bootstrap_IDs{k} = IDs_BaselineUnmedicated_statistics(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1) & ~isnan(weights_2year'));

end

%% Determine SRM of sensor data
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))

SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(trend_tremor_time_unmedicated_filled));
sensor_names = {'tremor_time_logit','trend_modal_tremor_power_unmedicated_filled','trend_perc90_tremor_power_unmedicated_filled'};
% sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 51;

% idx = contains([IDs_StartMedication; IDs_AllUnmedicated],IDs_BaselineUnmedicated_statistics);

weighted_SRM_sensor_unmedicated = [];
weighted_CI_sensor_unmedicated = [];
SRM_sensor_unmedicated = [];
CI_sensor_unmedicated = [];
N_sensor_unmedicated = [];
bootstat_sensor_unmedicated = [];

for k = 1:length(UPDRS_names)
    % data = eval(sensor_names{k});
    data = eval(sensor_names{1});

    idx = contains([IDs_StartMedication; IDs_AllUnmedicated],Bootstrap_IDs{k});

    change = data(idx,last_week_idx) - data(idx,first_week_idx);

    weights_boostrap = weights_2year(ismember(IDs_BaselineUnmedicated_statistics,Bootstrap_IDs{k}));
    % weights_boostrap = weights_2year;

    weighted_SRM_sensor_unmedicated(k) = weighted_SRM_function(change,weights_boostrap);
    [weighted_CI_sensor_unmedicated(k,1:2), bootstat_sensor_unmedicated(:,k)] = bootci(10000,{weighted_SRM_function,change,weights_boostrap});

    % SRM_sensor_unmedicated(k) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    % CI_sensor_unmedicated(k,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));

    % SE_mean_diff_data_PPP = std(data(idx,last_week_idx) - data(idx,1),'omitnan')/sqrt(length(find(~isnan(data(idx,last_week_idx))==1 & ~isnan(data(idx,1))==1)));
    % CI_sensor_data_low_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')-1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');
    % CI_sensor_data_high_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')+1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');

    N_sensor_unmedicated(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx)) & ~isnan(weights_boostrap)));
end

%% Testing significant difference
% Step 1: Bootstrap difference
diff_boot = bootstat_sensor_unmedicated - bootstat_UPDRS_unmedicated;

p_val = [];

% Step 2: Compute one-tailed bootstrap p-value
for i = 1:length(UPDRS_names)
    p_val(i) = length(find(diff_boot(:,i) >= 0))/10000;
end

% Step 4: Display results
fprintf('Bootstrap p-value: %.4f\n', p_val);

p_val = p_val';

effect_size = mean(diff_boot,'omitnan')
CI_effect_size = prctile(diff_boot, [2.5, 97.5])

