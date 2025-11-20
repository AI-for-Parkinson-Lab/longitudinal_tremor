%% 
clear all; close all;
%% Load sensor data
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_L1trend\trends2.mat')

%% Load IDs
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\IDs.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250410_matching_info\IDs_BaselineMedicated_matched.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_l1trend\IDs_BaselineUnmedicated_L1trend.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\Inclusion.mat')

for i = 1:length(IDs_BaselineUnmedicated_L1trend)
    id = IDs_BaselineUnmedicated_L1trend(i);
    StartWeek(i) = Inclusion.StartWeek(ismember(Inclusion.ID,id));
end

% IDs_AllUnmedicated_complete = IDs_AllUnmedicated(all(~isnan(trend_tremor_time_unmedicated_filled(1+length(IDs_StartMedication):end,[2 26])),2));

% IDs_BaselineMedicated_tremor = IDs_BaselineMedicated(all(~isnan(trend_median_tremor_power_medicated_filled(1:length(IDs_BaselineMedicated),[2 51])),2));
% IDs_BaselineMedicated_matched_tremor = intersect(IDs_BaselineMedicated_tremor,IDs_BaselineMedicated_matched,'stable');

IDs_BaselineUnmedicated = [IDs_StartMedication; IDs_AllUnmedicated];
IDs_1year_unmedicated = IDs_BaselineUnmedicated_L1trend(StartWeek>50);
IDs_BaselineUnmedicated_tremor = IDs_BaselineUnmedicated(all(~isnan(trend_median_tremor_power_unmedicated_filled(:,[2 51])),2));


%% Load clinical data
Visit1DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit1_DeNovo.csv");
Visit2DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit2_DeNovo.csv");
Visit3DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit3_DeNovo.csv");

Visit1PPP = readtable("C:\Users\z835211\Documents\Data\CSV files clinical and dempographic data\General_visit1.csv");
Visit2PPP = readtable("C:\Users\z835211\Documents\Data\CSV files clinical and dempographic data\General_visit2.csv");
Visit3PPP = readtable("C:\Users\z835211\Documents\Data\CSV files clinical and dempographic data\General_visit3.csv");

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

%% Extract UPDRS scores baseline medicated group (and disease duration)

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

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\visit_week_numbers.mat')


for i = 1:length(IDs_BaselineMedicated_matched)
    id = IDs_BaselineMedicated_matched{i};

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
    UPDRS_score3 = eval([UPDRS_names{k} '_3']);
    UPDRS_score1 = eval([UPDRS_names{k} '_1']);

    SRM_UPDRS_medicated_matched(k) = mean(UPDRS_score3-UPDRS_score1,'omitnan')/std(UPDRS_score3-UPDRS_score1,'omitnan');
    CI_UPDRS_medicated_matched(k,1:2) = bootci(10000,SRM_function,UPDRS_score3 - UPDRS_score1);

    % SE_UPDRS = std(UPDRS_score3-UPDRS_score1,'omitnan')/sqrt(length(find(~isnan(UPDRS_score1)==1 & ~isnan(UPDRS_score3)==1)));
    % CI_low_UPDRS(k) = (mean(UPDRS_score3-UPDRS_score1,'omitnan')-1.96*SE_UPDRS)/std(UPDRS_score3-UPDRS_score1,'omitnan');
    % CI_high_UPDRS(k) = (mean(UPDRS_score3-UPDRS_score1,'omitnan')+1.96*SE_UPDRS)/std(UPDRS_score3-UPDRS_score1,'omitnan');

    N_UPDRS_medicated_matched(k) = length(find(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1)));
end

%% Determine SRM of sensor data
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))
SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(trend_tremor_time_medicated_filled));
% sensor_names = {'tremor_time_logit','trend_modal_tremor_power_medicated_filled','trend_perc90_tremor_power_medicated_filled'};
sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 51;

idx = contains([IDs_BaselineMedicated; IDs_StartMedication],IDs_BaselineMedicated_matched);

SRM_sensor_medicated_matched = [];
CI_sensor_medicated_matched = [];
N_sensor_medicated_matched = [];

for k = 1:length(sensor_names)
    data = eval(sensor_names{k});
    SRM_sensor_medicated_matched(k) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    CI_sensor_medicated_matched(k,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));

    % SE_mean_diff_data_PPP = std(data(idx,last_week_idx) - data(idx,1),'omitnan')/sqrt(length(find(~isnan(data(idx,last_week_idx))==1 & ~isnan(data(idx,1))==1)));
    % CI_sensor_data_low_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')-1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');
    % CI_sensor_data_high_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')+1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');

    N_sensor_medicated_matched(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx))));
end

%% Extract UPDRS scores baseline medicated group (and disease duration)

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

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\visit_week_numbers.mat')


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

SRM_UPDRS_medicated_unmatched = [];
CI_UPDRS_medicated_unmatched = [];
N_UPDRS_medicated_unmatched = [];

for k = 1:length(UPDRS_names)
    UPDRS_score3 = eval([UPDRS_names{k} '_3']);
    UPDRS_score1 = eval([UPDRS_names{k} '_1']);

    SRM_UPDRS_medicated_unmatched(k) = mean(UPDRS_score3-UPDRS_score1,'omitnan')/std(UPDRS_score3-UPDRS_score1,'omitnan');
    CI_UPDRS_medicated_unmatched(k,1:2) = bootci(10000,SRM_function,UPDRS_score3 - UPDRS_score1);

    % SE_UPDRS = std(UPDRS_score3-UPDRS_score1,'omitnan')/sqrt(length(find(~isnan(UPDRS_score1)==1 & ~isnan(UPDRS_score3)==1)));
    % CI_low_UPDRS(k) = (mean(UPDRS_score3-UPDRS_score1,'omitnan')-1.96*SE_UPDRS)/std(UPDRS_score3-UPDRS_score1,'omitnan');
    % CI_high_UPDRS(k) = (mean(UPDRS_score3-UPDRS_score1,'omitnan')+1.96*SE_UPDRS)/std(UPDRS_score3-UPDRS_score1,'omitnan');

    N_UPDRS_medicated_unmatched(k) = length(find(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1)));
end

%% Determine SRM of sensor data
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))
SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(trend_tremor_time_medicated_filled));
% sensor_names = {'tremor_time_logit','trend_modal_tremor_power_medicated_filled','trend_perc90_tremor_power_medicated_filled'};
sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 51;

idx = contains([IDs_BaselineMedicated; IDs_StartMedication],IDs_BaselineMedicated);

SRM_sensor_medicated_unmatched = [];
CI_sensor_medicated_unmatched = [];
N_sensor_medicated_unmatched = [];

for k = 1:length(sensor_names)
    data = eval(sensor_names{k});
    SRM_sensor_medicated_unmatched(k) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    CI_sensor_medicated_unmatched(k,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));

    % SE_mean_diff_data_PPP = std(data(idx,last_week_idx) - data(idx,1),'omitnan')/sqrt(length(find(~isnan(data(idx,last_week_idx))==1 & ~isnan(data(idx,1))==1)));
    % CI_sensor_data_low_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')-1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');
    % CI_sensor_data_high_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')+1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');

    N_sensor_medicated_unmatched(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx))));
end

%% Create forest plot
% Names = {'rest tremor severity in OFF (matched', 'rest tremor severity in OFF (unmatched','rest tremor constancy in OFF (matched','rest tremor constancy in OFF (unmatched',...
%     'rest tremor severity in ON (matched','rest tremor severity in ON (unmatched', 'rest tremor constancy in ON (matched', 'rest tremor constancy in ON (unmatched',...
%     'patient-reported tremor (matched','patient-reported tremor (unmatched','tremor time (matched','tremor time (unmatched','modal tremor power (matched',...
%     'modal tremor power (unmatched','90th percentile of tremor power (matched','90th percentile of tremor power (unmatched',}; % Study names]
Names = {'rest tremor severity in OFF (matched', 'rest tremor severity in OFF (unmatched','rest tremor constancy in OFF (matched','rest tremor constancy in OFF (unmatched',...
    'rest tremor severity in ON (matched','rest tremor severity in ON (unmatched', 'rest tremor constancy in ON (matched', 'rest tremor constancy in ON (unmatched',...
    'patient-reported tremor (matched','patient-reported tremor (unmatched','tremor time (matched','tremor time (unmatched'}; % Study names]
counts = [ ...
    reshape([N_UPDRS_medicated_matched;   N_UPDRS_medicated_unmatched], 1, []), ...
    reshape([N_sensor_medicated_matched;  N_sensor_medicated_unmatched], 1, []) ...
];
Names = strcat(Names, ", n = ", string(counts), ")");
SRM = [reshape([SRM_UPDRS_medicated_matched; SRM_UPDRS_medicated_unmatched],1,[]) reshape([SRM_sensor_medicated_matched; SRM_sensor_medicated_unmatched],1,[])]; % Standardized Response Means
CI = [reshape([CI_UPDRS_medicated_matched.'; CI_UPDRS_medicated_unmatched.'],2,[]).'; reshape([CI_sensor_medicated_matched.'; CI_sensor_medicated_unmatched.'],2,[]).']; % Lower bounds of 95% confidence intervals

% Plotting
figure;
hold on;


% Number of studies
n = length(SRM);

% Plot the confidence intervals
for i = 1:n
    if i == 6 || i == 10
    % if i == 6 || i == 7 || i == 8 || i == 12 || i == 13 || i == 14
        line([CI(i,1), CI(i,2)], [n-i+1, n-i+1], 'Color', 'r', 'LineWidth', 1.5); % CI line
        plot(SRM(i), n-i+1, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % SRM point
    else
        line([CI(i,1), CI(i,2)], [n-i+1, n-i+1], 'Color', 'k', 'LineWidth', 1.5); % CI line
        plot(SRM(i), n-i+1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % SRM point
    end
end

% Add labels
yticks(1:n)
yticklabels(Names(end:-1:1))
xlabel('SRM with 95%-CI');
title('2-year change')

xlim([-1 1]);
ylim([0 n+1])

% Reference line at zero effect
line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

hold off;

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


for i = 1:length(IDs_BaselineUnmedicated_L1trend)
    id = IDs_BaselineUnmedicated_L1trend{i};

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
unweighted_SRM_UPDRS_unmedicated = [];
unweighted_CI_UPDRS_unmedicated = [];

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\v3_survival\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(77,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100

weights(~isnan(weights)) = 1; % Show unweighted results

idx_2year = contains(IDs_BaselineUnmedicated_L1trend,IDs_BaselineUnmedicated_L1trend);
last_week_idx = 51;
weights_2year = weights(idx_2year,last_week_idx);

for k = 1:length(UPDRS_names)
    UPDRS_score3 = eval([UPDRS_names{k} '_3']);
    UPDRS_score1 = eval([UPDRS_names{k} '_1']);

    change = (UPDRS_score3-UPDRS_score1)';

    unweighted_SRM_UPDRS_unmedicated(k) = weighted_SRM_function(change,weights_2year);
    unweighted_CI_UPDRS_unmedicated(k,1:2) = bootci(10000,{weighted_SRM_function,change,weights_2year});

    % SRM_UPDRS_unmedicated(k) = mean(UPDRS_score3-UPDRS_score1,'omitnan')/std(UPDRS_score3-UPDRS_score1,'omitnan');
    % CI_UPDRS_unmedicated(k,1:2) = bootci(10000,SRM_function,UPDRS_score3 - UPDRS_score1);

    N_UPDRS_unmedicated(k) = length(find(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1) & ~isnan(weights_2year')));
end

%% Determine SRM of sensor data
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))

SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(trend_tremor_time_unmedicated_filled));
% sensor_names = {'tremor_time_logit','trend_modal_tremor_power_unmedicated_filled','trend_perc90_tremor_power_unmedicated_filled'};
sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 51;

idx = contains([IDs_StartMedication; IDs_AllUnmedicated],IDs_BaselineUnmedicated_L1trend);

unweighted_SRM_sensor_unmedicated = [];
unweighted_CI_sensor_unmedicated = [];
N_sensor_unmedicated = [];

for k = 1:length(sensor_names)
    data = eval(sensor_names{k});

    change = data(idx,last_week_idx) - data(idx,first_week_idx);

    unweighted_SRM_sensor_unmedicated(k) = weighted_SRM_function(change,weights_2year);
    unweighted_CI_sensor_unmedicated(k,1:2) = bootci(10000,{weighted_SRM_function,change,weights_2year});

    % SRM_sensor_unmedicated(k) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    % CI_sensor_unmedicated(k,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));

    % SE_mean_diff_data_PPP = std(data(idx,last_week_idx) - data(idx,1),'omitnan')/sqrt(length(find(~isnan(data(idx,last_week_idx))==1 & ~isnan(data(idx,1))==1)));
    % CI_sensor_data_low_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')-1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');
    % CI_sensor_data_high_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')+1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');

    N_sensor_unmedicated(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx)) & ~isnan(weights_2year)));
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

load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\v3_survival\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(77,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100

idx_2year = contains(IDs_BaselineUnmedicated_L1trend,IDs_BaselineUnmedicated_L1trend);
last_week_idx = 51;
weights_2year = weights(idx_2year,last_week_idx);

for k = 1:length(UPDRS_names)
    UPDRS_score3 = eval([UPDRS_names{k} '_3']);
    UPDRS_score1 = eval([UPDRS_names{k} '_1']);

    change = (UPDRS_score3-UPDRS_score1)';

    weighted_SRM_UPDRS_unmedicated(k) = weighted_SRM_function(change,weights_2year);
    weighted_CI_UPDRS_unmedicated(k,1:2) = bootci(10000,{weighted_SRM_function,change,weights_2year});

    % SRM_UPDRS_unmedicated(k) = mean(UPDRS_score3-UPDRS_score1,'omitnan')/std(UPDRS_score3-UPDRS_score1,'omitnan');
    % CI_UPDRS_unmedicated(k,1:2) = bootci(10000,SRM_function,UPDRS_score3 - UPDRS_score1);

    N_UPDRS_unmedicated(k) = length(find(~isnan(UPDRS_score3) & ~isnan(UPDRS_score1) & ~isnan(weights_2year')));
end

%% Determine SRM of sensor data
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))

SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

tremor_time_logit = real(logit(trend_tremor_time_unmedicated_filled));
% sensor_names = {'tremor_time_logit','trend_modal_tremor_power_unmedicated_filled','trend_perc90_tremor_power_unmedicated_filled'};
sensor_names = {'tremor_time_logit'};
first_week_idx = 2;
last_week_idx = 51;

idx = contains([IDs_StartMedication; IDs_AllUnmedicated],IDs_BaselineUnmedicated_L1trend);

weighted_SRM_sensor_unmedicated = [];
weighted_CI_sensor_unmedicated = [];
N_sensor_unmedicated = [];

for k = 1:length(sensor_names)
    data = eval(sensor_names{k});

    change = data(idx,last_week_idx) - data(idx,first_week_idx);

    weighted_SRM_sensor_unmedicated(k) = weighted_SRM_function(change,weights_2year);
    weighted_CI_sensor_unmedicated(k,1:2) = bootci(10000,{weighted_SRM_function,change,weights_2year});

    % SRM_sensor_unmedicated(k) = mean(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan')/std(data(idx,last_week_idx) - data(idx,first_week_idx),'omitnan');
    % CI_sensor_unmedicated(k,1:2) = bootci(10000,SRM_function,data(idx,last_week_idx) - data(idx,first_week_idx));

    % SE_mean_diff_data_PPP = std(data(idx,last_week_idx) - data(idx,1),'omitnan')/sqrt(length(find(~isnan(data(idx,last_week_idx))==1 & ~isnan(data(idx,1))==1)));
    % CI_sensor_data_low_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')-1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');
    % CI_sensor_data_high_PPP(k) = (mean(data(idx,last_week_idx) - data(idx,1),'omitnan')+1.96*SE_mean_diff_data_PPP)/std(data(idx,last_week_idx) - data(idx,1),'omitnan');

    N_sensor_unmedicated(k) = length(find(~isnan(data(idx,last_week_idx)) & ~isnan(data(idx,first_week_idx)) & ~isnan(weights_2year)));
end

%% Create forest plot
% Names = {'rest tremor severity in OFF (unweighted', 'rest tremor severity in OFF (weighted','rest tremor constancy in OFF (unweighted','rest tremor constancy in OFF (weighted',...
%     'patient-reported tremor (unweighted','patient-reported tremor (weighted','tremor time (unweighted','tremor time (weighted',...
%     'modal tremor power (unweighted','modal tremor power (weighted','90th percentile of tremor power (unweighted','90th percentile of tremor power (weighted'}; % Study names]
Names = {'rest tremor severity in OFF (unweighted', 'rest tremor severity in OFF (weighted','rest tremor constancy in OFF (unweighted','rest tremor constancy in OFF (weighted',...
    'patient-reported tremor (unweighted','patient-reported tremor (weighted','tremor time (unweighted','tremor time (weighted'}; % Study names]
counts = [reshape([N_UPDRS_unmedicated; N_UPDRS_unmedicated],1,[]) reshape([N_sensor_unmedicated; N_sensor_unmedicated],1,[])];
Names = strcat(Names, ", n = ", string(counts), ")");
SRM = [reshape([unweighted_SRM_UPDRS_unmedicated; weighted_SRM_UPDRS_unmedicated],1,[]) reshape([unweighted_SRM_sensor_unmedicated; weighted_SRM_sensor_unmedicated],1,[])]; % Standardized Response Means
CI = [reshape([unweighted_CI_UPDRS_unmedicated.'; weighted_CI_UPDRS_unmedicated.'],2,[]).'; reshape([unweighted_CI_sensor_unmedicated.'; weighted_CI_sensor_unmedicated.'],2,[]).']; % Lower bounds of 95% confidence intervals

% Plotting
figure;
hold on;


% Number of studies
n = length(SRM);

% Plot the confidence intervals
for i = 1:n
    if i == 6 || i == 10
    % if i == 6 || i == 7 || i == 8 || i == 12 || i == 13 || i == 14
        line([CI(i,1), CI(i,2)], [n-i+1, n-i+1], 'Color', 'r', 'LineWidth', 1.5); % CI line
        plot(SRM(i), n-i+1, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % SRM point
    else
        line([CI(i,1), CI(i,2)], [n-i+1, n-i+1], 'Color', 'k', 'LineWidth', 1.5); % CI line
        plot(SRM(i), n-i+1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % SRM point
    end
end

% Add labels
yticks(1:n)
yticklabels(Names(end:-1:1))
xlabel('SRM with 95%-CI');
title('2-year change')

xlim([-1 1.7]);
ylim([0 n+1])

% Reference line at zero effect
line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

hold off;