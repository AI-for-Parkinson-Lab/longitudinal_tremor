%% Step 6c: Sensitivity analyses for SRM
% Assess the SRMs in the matched medicated group and unweighted unmedicated
% group

clear all; close all;

%% Load sensor data, IDs and descriptives
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Trends_filled.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IDs_selected.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Descriptives.mat')

%% Perform propensity score matching

treatment = [zeros(1,length(IDs_BaselineMedicated)) ones(1,length(IDs_BaselineUnmedicated))]';
X  = [BaselineMedicated(:,[2,4,5]).Variables; BaselineUnmedicated(:,[2,4,5]).Variables];

% Step 1: Estimate propensity scores using logistic regression
mdl = fitglm(X, treatment, 'Distribution', 'binomial');
propensity_scores = predict(mdl, X);

% Step 2: Nearest neighbor matching (1:1 without replacement)
treated_idx = find(treatment == 1);
control_idx = find(treatment == 0);

matched_pairs = [];
used_controls = false(length(control_idx), 1);

for i = 1:length(treated_idx)
    t_idx = treated_idx(i);
    ps_treated = propensity_scores(t_idx);
    
    % Get available controls
    available_controls = find(~used_controls);
    ps_controls = propensity_scores(control_idx(available_controls));
    
    % Find nearest neighbor
    [~, min_idx] = min(abs(ps_controls - ps_treated));
    c_idx = control_idx(available_controls(min_idx));
    
    % Store match
    matched_pairs = [matched_pairs; t_idx, c_idx];
    used_controls(available_controls(min_idx)) = true;
end

IDs_BaselineMedicated_matched = IDs_BaselineMedicated(matched_pairs(:,2))

save('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IDs_BaselineMedicated_matched.mat',"IDs_BaselineMedicated_matched")

%% Extract descriptives
BaselineMedicated_matched = BaselineMedicated(ismember(IDs_BaselineMedicated,IDs_BaselineMedicated_matched),:);

%% Age
median(BaselineMedicated_matched.Age)
prctile(BaselineMedicated_matched.Age,25)
prctile(BaselineMedicated_matched.Age,75)

%% Gender
length(find(BaselineMedicated_matched.Gender==1))

%% Disease duration
median(BaselineMedicated_matched.Disease_duration,'omitnan')
prctile(BaselineMedicated_matched.Disease_duration,25)
prctile(BaselineMedicated_matched.Disease_duration,75)

%% Most-affected
n1 = length(find(BaselineMedicated_matched.Most_affected==1))
n2 = length(find(BaselineMedicated_matched.Most_affected==0))

%% UPDRS 1
median(BaselineMedicated_matched.UPDRS1,'omitnan')
prctile(BaselineMedicated_matched.UPDRS1,25)
prctile(BaselineMedicated_matched.UPDRS1,75)

%% UPDRS 2
median(BaselineMedicated_matched.UPDRS2,'omitnan')
prctile(BaselineMedicated_matched.UPDRS2,25)
prctile(BaselineMedicated_matched.UPDRS2,75)

%% UPDRS 3
median(BaselineMedicated_matched.UPDRS3_OFF,'omitnan')
prctile(BaselineMedicated_matched.UPDRS3_OFF,25)
prctile(BaselineMedicated_matched.UPDRS3_OFF,75)

%% UPDRS 4
median(BaselineMedicated_matched.UPDRS4,'omitnan')
prctile(BaselineMedicated_matched.UPDRS4,25)
prctile(BaselineMedicated_matched.UPDRS4,75)

%% UPDRS 3 tremor subscore
median(BaselineMedicated_matched.UPDRS_tremor_OFF,'omitnan')
prctile(BaselineMedicated_matched.UPDRS_tremor_OFF,25)
prctile(BaselineMedicated_matched.UPDRS_tremor_OFF,75)

%% Rest tremor device sided arm
median(BaselineMedicated_matched.UPDRS_devside_tremor_OFF,'omitnan')
prctile(BaselineMedicated_matched.UPDRS_devside_tremor_OFF,25)
prctile(BaselineMedicated_matched.UPDRS_devside_tremor_OFF,75)

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

%% Select participants with tremor
IDs_BaselineMedicated_tremor = IDs_BaselineMedicated(all(~isnan(trend_modal_tremor_power_medicated_filled(1:length(IDs_BaselineMedicated),[2 26])),2));
IDs_BaselineMedicated_matched_tremor = intersect(IDs_BaselineMedicated_tremor,IDs_BaselineMedicated_matched,'stable');

%% Extract UPDRS scores matched medicated group

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

load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\visit_week_numbers.mat')

for i = 1:length(IDs_BaselineMedicated_matched_tremor)
    id = IDs_BaselineMedicated_matched_tremor{i};

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

SRM_UPDRS_medicated_matched = [];
CI_UPDRS_medicated_matched = [];
N_UPDRS_medicated_matched = [];

for k = 1:length(UPDRS_names)
    UPDRS_score_followup = eval([UPDRS_names{k} '_2']); % score 2 or 3
    UPDRS_score_baseline = eval([UPDRS_names{k} '_1']);

    SRM_UPDRS_medicated_matched(k) = mean(UPDRS_score_followup-UPDRS_score_baseline,'omitnan')/std(UPDRS_score_followup-UPDRS_score_baseline,'omitnan');
    CI_UPDRS_medicated_matched(k,1:2) = bootci(10000,SRM_function,UPDRS_score_followup - UPDRS_score_baseline);
    N_UPDRS_medicated_matched(k) = length(find(~isnan(UPDRS_score_followup) & ~isnan(UPDRS_score_baseline)));

end

%% Determine SRM of sensor data

tremor_time_logit = real(logit(trend_tremor_time_medicated_filled));

% Select measures (tremor time for full group, all measures for tremor group)
sensor_names = {'tremor_time_logit','trend_modal_tremor_power_medicated_filled','trend_perc90_tremor_power_medicated_filled'};
% sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 26; % Modify for one- or two-year follow-up

idx = contains([IDs_BaselineMedicated; IDs_BaselineUnmedicated],IDs_BaselineMedicated_matched_tremor); % Select full group or tremor group

SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

SRM_sensor_medicated_matched = [];
CI_sensor_medicated_matched = [];
N_sensor_medicated_matched = [];

for k = 1:length(sensor_names)    

    data = eval(sensor_names{k});
    
    SRM_sensor_medicated_matched(k) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    CI_sensor_medicated_matched(k,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));
    N_sensor_medicated_matched(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx))));
end

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

for k = 1:length(UPDRS_names)
    UPDRS_score_followup = eval([UPDRS_names{k} '_2']); % score 2 or 3
    UPDRS_score_baseline = eval([UPDRS_names{k} '_1']);

    SRM_UPDRS_medicated(k) = mean(UPDRS_score_followup-UPDRS_score_baseline,'omitnan')/std(UPDRS_score_followup-UPDRS_score_baseline,'omitnan');
    CI_UPDRS_medicated(k,1:2) = bootci(10000,SRM_function,UPDRS_score_followup - UPDRS_score_baseline);

    N_UPDRS_medicated(k) = length(find(~isnan(UPDRS_score_followup) & ~isnan(UPDRS_score_baseline)));
end

%% Determine SRM of sensor data
addpath(genpath('utils'))
tremor_time_logit = real(logit(trend_tremor_time_medicated_filled));

% Select measures (tremor time for full group, all measures for tremor group)
sensor_names = {'tremor_time_logit','trend_modal_tremor_power_medicated_filled','trend_perc90_tremor_power_medicated_filled'};
% sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 26; % Modify for one- or two-year follow-up

idx = contains([IDs_BaselineMedicated; IDs_BaselineUnmedicated],IDs_BaselineMedicated_tremor);

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

%% Create forest plot
Names = {'rest tremor severity in OFF (matched', 'rest tremor severity in OFF (unmatched','rest tremor constancy in OFF (matched','rest tremor constancy in OFF (unmatched',...
    'rest tremor severity in ON (matched','rest tremor severity in ON (unmatched', 'rest tremor constancy in ON (matched', 'rest tremor constancy in ON (unmatched',...
    'patient-reported tremor (matched','patient-reported tremor (unmatched','tremor time (matched','tremor time (unmatched','modal tremor power (matched',...
    'modal tremor power (unmatched','90th percentile of tremor power (matched','90th percentile of tremor power (unmatched',}; % Study names]
% Names = {'rest tremor severity in OFF (matched', 'rest tremor severity in OFF (unmatched','rest tremor constancy in OFF (matched','rest tremor constancy in OFF (unmatched',...
%     'rest tremor severity in ON (matched','rest tremor severity in ON (unmatched', 'rest tremor constancy in ON (matched', 'rest tremor constancy in ON (unmatched',...
%     'patient-reported tremor (matched','patient-reported tremor (unmatched','tremor time (matched','tremor time (unmatched'}; % Study names]
counts = [ ...
    reshape([N_UPDRS_medicated_matched;   N_UPDRS_medicated], 1, []), ...
    reshape([N_sensor_medicated_matched;  N_sensor_medicated], 1, []) ...
];
Names = strcat(Names, ", n = ", string(counts), ")");
SRM = [reshape([SRM_UPDRS_medicated_matched; SRM_UPDRS_medicated],1,[]) reshape([SRM_sensor_medicated_matched; SRM_sensor_medicated],1,[])]; % Standardized Response Means
CI = [reshape([CI_UPDRS_medicated_matched.'; CI_UPDRS_medicated.'],2,[]).'; reshape([CI_sensor_medicated_matched.'; CI_sensor_medicated.'],2,[]).']; % Lower bounds of 95% confidence intervals

% Plotting
figure;
hold on;
C = colororder('glow12');

% Number of studies
n = length(SRM);

% Plot the confidence intervals
for i = 1:n
    % if i == 11 || i == 12
    if i == 11 || i == 12 || i == 13 || i == 14 || i == 15 || i == 16
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

xlim([-1 1]);
ylim([0 n+1])

% Reference line at zero effect
line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

hold off;

%% Extract UPDRS scores unmedicated group (no ON assessment) 

IDs_BaselineUnmedicated_tremor = IDs_BaselineUnmedicated(all(~isnan(trend_modal_tremor_power_unmedicated_filled(:,[2 26])),2)); % Change index for one- or two-year group


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
%% Determine SRM of UPDRS scores (weighted)

load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IPCW.mat')

weighted_SRM_function = @(change, weights) ...
    (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan')) / ...
    sqrt(sum(weights .* (change - (sum(weights .* change, 'omitnan') / sum(weights, 'omitnan'))).^2, 'omitnan') / sum(weights, 'omitnan'));

UPDRS_names = {'UPDRS_317OFF','UPDRS_318OFF','UPDRS_210'};

weighted_SRM_UPDRS_unmedicated = [];
weighted_CI_UPDRS_unmedicated = [];
weighted_N_UPDRS_unmedicated = [];

weights = IPCW(:,2:end).Variables;
weights = [ones(78,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100, duplicate to increase availability of weights

idx_followup = contains(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated); % Select full group or tremor group
last_week_idx = 26;
weights_followup = weights(idx_followup,last_week_idx);

for k = 1:length(UPDRS_names)

    UPDRS_score_followup = eval([UPDRS_names{k} '_2']);
    UPDRS_score_baseline = eval([UPDRS_names{k} '_1']);

    change = (UPDRS_score_followup-UPDRS_score_baseline)';

    weighted_SRM_UPDRS_unmedicated(k) = weighted_SRM_function(change,weights_followup);
    weighted_CI_UPDRS_unmedicated(k,1:2)= bootci(10000,{weighted_SRM_function,change,weights_followup});
    weighted_N_UPDRS_unmedicated(k) = length(find(~isnan(UPDRS_score_followup) & ~isnan(UPDRS_score_baseline) & ~isnan(weights_followup')));
end

%% Determine SRM of sensor data (weighted)

tremor_time_logit = real(logit(trend_tremor_time_unmedicated_filled));

% Select measures (tremor time for full group, all measures for tremor group)
% sensor_names = {'tremor_time_logit','trend_modal_tremor_power_unmedicated_filled','trend_perc90_tremor_power_unmedicated_filled'};
sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 26;

idx = contains(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated);

weighted_SRM_sensor_unmedicated = [];
weighted_CI_sensor_unmedicated = [];
weighted_N_sensor_unmedicated = [];

for k = 1:length(sensor_names)

    data = eval(sensor_names{k});

    change = data(idx,last_week_idx) - data(idx,first_week_idx);

    weighted_SRM_sensor_unmedicated(k) = weighted_SRM_function(change,weights_followup);
    weighted_CI_sensor_unmedicated(k,1:2) = bootci(10000,{weighted_SRM_function,change,weights_followup});
    weighted_N_sensor_unmedicated(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx)) & ~isnan(weights_followup)));
end

%% Determine SRM of UPDRS scores (unweighted)

weights(~isnan(weights)) = 1; % Show unweighted results

UPDRS_names = {'UPDRS_317OFF','UPDRS_318OFF','UPDRS_210'};

unweighted_SRM_UPDRS_unmedicated = [];
unweighted_CI_UPDRS_unmedicated = [];
unweighted_N_UPDRS_unmedicated = [];

idx_followup = contains(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated);
last_week_idx = 26;
weights_followup = weights(idx_followup,last_week_idx);

for k = 1:length(UPDRS_names)
    UPDRS_score_followup = eval([UPDRS_names{k} '_2']);
    UPDRS_score_baseline = eval([UPDRS_names{k} '_1']);

    change = (UPDRS_score_followup-UPDRS_score_baseline)';

    unweighted_SRM_UPDRS_unmedicated(k) = weighted_SRM_function(change,weights_followup);
    unweighted_CI_UPDRS_unmedicated(k,1:2) = bootci(10000,{weighted_SRM_function,change,weights_followup});
    unweighted_N_UPDRS_unmedicated(k) = length(find(~isnan(UPDRS_score_followup) & ~isnan(UPDRS_score_baseline) & ~isnan(weights_followup')));
end

%% Determine SRM of sensor data (unweighted)

tremor_time_logit = real(logit(trend_tremor_time_unmedicated_filled));
% sensor_names = {'tremor_time_logit','trend_modal_tremor_power_unmedicated_filled','trend_perc90_tremor_power_unmedicated_filled'};
sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 26;

idx = contains(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated);

unweighted_SRM_sensor_unmedicated = [];
unweighted_CI_sensor_unmedicated = [];
unweighted_N_sensor_unmedicated = [];

for k = 1:length(sensor_names)
    data = eval(sensor_names{k});

    change = data(idx,last_week_idx) - data(idx,first_week_idx);

    unweighted_SRM_sensor_unmedicated(k) = weighted_SRM_function(change,weights_followup);
    unweighted_CI_sensor_unmedicated(k,1:2) = bootci(10000,{weighted_SRM_function,change,weights_followup});
    unweighted_N_sensor_unmedicated(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx)) & ~isnan(weights_followup)));
end

%% Create forest plot
% Names = {'rest tremor severity in OFF (unweighted', 'rest tremor severity in OFF (weighted','rest tremor constancy in OFF (unweighted','rest tremor constancy in OFF (weighted',...
%     'patient-reported tremor (unweighted','patient-reported tremor (weighted','tremor time (unweighted','tremor time (weighted',...
%     'modal tremor power (unweighted','modal tremor power (weighted','90th percentile of tremor power (unweighted','90th percentile of tremor power (weighted'}; % Study names]
Names = {'rest tremor severity in OFF (unweighted', 'rest tremor severity in OFF (weighted','rest tremor constancy in OFF (unweighted','rest tremor constancy in OFF (weighted',...
    'patient-reported tremor (unweighted','patient-reported tremor (weighted','tremor time (unweighted','tremor time (weighted'}; % Study names]
counts = [reshape([unweighted_N_UPDRS_unmedicated; weighted_N_UPDRS_unmedicated],1,[]) reshape([unweighted_N_sensor_unmedicated; weighted_N_sensor_unmedicated],1,[])];
Names = strcat(Names, ", n = ", string(counts), ")");
SRM = [reshape([unweighted_SRM_UPDRS_unmedicated; weighted_SRM_UPDRS_unmedicated],1,[]) reshape([unweighted_SRM_sensor_unmedicated; weighted_SRM_sensor_unmedicated],1,[])]; % Standardized Response Means
CI = [reshape([unweighted_CI_UPDRS_unmedicated.'; weighted_CI_UPDRS_unmedicated.'],2,[]).'; reshape([unweighted_CI_sensor_unmedicated.'; weighted_CI_sensor_unmedicated.'],2,[]).']; % Lower bounds of 95% confidence intervals

% Plotting
figure;
hold on;
C = colororder('glow12');

% Number of studies
n = length(SRM);

% Plot the confidence intervals
for i = 1:n
    if i == 7 || i == 8
    % if i == 7 || i == 8 || i == 9 || i == 10 || i == 11 || i == 12
        line([CI(i,1), CI(i,2)], [n-i+1, n-i+1], 'Color', 'r', 'LineWidth', 1.5); % CI line
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

xlim([-1 1.7]);
ylim([0 n+1])

% Reference line at zero effect
line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

hold off;