%% 
clear all, close all;

%% Load sensor data
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_L1trend\trends2.mat')

%% Load IDs
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\IDs.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250410_matching_info\IDs_BaselineMedicated_matched.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_l1trend\IDs_BaselineUnmedicated_L1trend.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\Inclusion.mat')

%% Extract start week, most-affected side and disease-duration
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250410_matching_info\matching_info.mat')

for i = 1:length(IDs_BaselineUnmedicated_L1trend)
    id = IDs_BaselineUnmedicated_L1trend(i);
    StartWeek(i) = Inclusion.StartWeek(ismember(Inclusion.ID,id));
    Most_affected(i) = BaselineUnmedicated.Most_affected(ismember(BaselineUnmedicated.ID,id));
    Disease_duration(i) = BaselineUnmedicated.Disease_duration(ismember(BaselineUnmedicated.ID,id));
end

%% Calculate moving average of 4 weeks
IDs_BaselineUnmedicated = [IDs_StartMedication; IDs_AllUnmedicated];

tremor_time = trend_tremor_time_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_BaselineUnmedicated_L1trend),:);

numSubjects = size(tremor_time,1);
times = (0:2:104)';
data = [];

for i = 1:numSubjects
    event = (times < StartWeek(i) - 1)';
    X = tremor_time(i,event)'; 
    data = [data; repmat(i, length(X), 1), repmat(Most_affected(i), length(X), 1), repmat(Disease_duration(i), length(X), 1), times(1:length(X)), X, event(1:length(X))'];
    if length(X)<53
        data = [data; i, Most_affected(i), Disease_duration(i), times(length(X)+1), tremor_time(i,length(X)+1), 0];
    end
end

% Convert to table for easier handling
dataTable = array2table(data, 'VariableNames', {'ID','Most-affected','Disease duration', 'Time', 'X', 'Event'});
%% Implement time lags
windowSize = 2; % 4-week moving average
for i = 1:numSubjects
    movavg = movmean(dataTable.X(dataTable.ID==i), windowSize, 'Endpoints', 'fill'); 
    dataTable.X_MA(dataTable.ID==i) = movavg;
    if length(movavg)>9
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = [NaN(2,1); movavg(1:end-2)];
        dataTable.X_MA_Lagged3(dataTable.ID==i) = [NaN(3,1); movavg(1:end-3)];
        dataTable.X_MA_Lagged4(dataTable.ID==i) = [NaN(4,1); movavg(1:end-4)];
        dataTable.X_MA_Lagged5(dataTable.ID==i) = [NaN(5,1); movavg(1:end-5)];
        dataTable.X_MA_Lagged6(dataTable.ID==i) = [NaN(6,1); movavg(1:end-6)];
        dataTable.X_MA_Lagged7(dataTable.ID==i) = [NaN(7,1); movavg(1:end-7)];
        dataTable.X_MA_Lagged8(dataTable.ID==i) = [NaN(8,1); movavg(1:end-8)];
        dataTable.X_MA_Lagged9(dataTable.ID==i) = [NaN(9,1); movavg(1:end-9)];
    elseif length(movavg)>8
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = [NaN(2,1); movavg(1:end-2)];
        dataTable.X_MA_Lagged3(dataTable.ID==i) = [NaN(3,1); movavg(1:end-3)];
        dataTable.X_MA_Lagged4(dataTable.ID==i) = [NaN(4,1); movavg(1:end-4)];
        dataTable.X_MA_Lagged5(dataTable.ID==i) = [NaN(5,1); movavg(1:end-5)];
        dataTable.X_MA_Lagged6(dataTable.ID==i) = [NaN(6,1); movavg(1:end-6)];
        dataTable.X_MA_Lagged7(dataTable.ID==i) = [NaN(7,1); movavg(1:end-7)];
        dataTable.X_MA_Lagged8(dataTable.ID==i) = [NaN(8,1); movavg(1:end-8)];
        dataTable.X_MA_Lagged9(dataTable.ID==i) = NaN(9,1);
    elseif length(movavg)>7
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = [NaN(2,1); movavg(1:end-2)];
        dataTable.X_MA_Lagged3(dataTable.ID==i) = [NaN(3,1); movavg(1:end-3)];
        dataTable.X_MA_Lagged4(dataTable.ID==i) = [NaN(4,1); movavg(1:end-4)];
        dataTable.X_MA_Lagged5(dataTable.ID==i) = [NaN(5,1); movavg(1:end-5)];
        dataTable.X_MA_Lagged6(dataTable.ID==i) = [NaN(6,1); movavg(1:end-6)];
        dataTable.X_MA_Lagged7(dataTable.ID==i) = [NaN(7,1); movavg(1:end-7)];
        dataTable.X_MA_Lagged8(dataTable.ID==i) = NaN(8,1);
        dataTable.X_MA_Lagged9(dataTable.ID==i) = NaN(8,1);
    elseif length(movavg)>6
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = [NaN(2,1); movavg(1:end-2)];
        dataTable.X_MA_Lagged3(dataTable.ID==i) = [NaN(3,1); movavg(1:end-3)];
        dataTable.X_MA_Lagged4(dataTable.ID==i) = [NaN(4,1); movavg(1:end-4)];
        dataTable.X_MA_Lagged5(dataTable.ID==i) = [NaN(5,1); movavg(1:end-5)];
        dataTable.X_MA_Lagged6(dataTable.ID==i) = [NaN(6,1); movavg(1:end-6)];
        dataTable.X_MA_Lagged7(dataTable.ID==i) = NaN(7,1);
        dataTable.X_MA_Lagged8(dataTable.ID==i) = NaN(7,1);
        dataTable.X_MA_Lagged9(dataTable.ID==i) = NaN(7,1);
    elseif length(movavg)>5
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = [NaN(2,1); movavg(1:end-2)];
        dataTable.X_MA_Lagged3(dataTable.ID==i) = [NaN(3,1); movavg(1:end-3)];
        dataTable.X_MA_Lagged4(dataTable.ID==i) = [NaN(4,1); movavg(1:end-4)];
        dataTable.X_MA_Lagged5(dataTable.ID==i) = [NaN(5,1); movavg(1:end-5)];
        dataTable.X_MA_Lagged6(dataTable.ID==i) = NaN(6,1);
        dataTable.X_MA_Lagged7(dataTable.ID==i) = NaN(6,1);
        dataTable.X_MA_Lagged8(dataTable.ID==i) = NaN(6,1);
        dataTable.X_MA_Lagged9(dataTable.ID==i) = NaN(6,1);
    elseif length(movavg)>4
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = [NaN(2,1); movavg(1:end-2)];
        dataTable.X_MA_Lagged3(dataTable.ID==i) = [NaN(3,1); movavg(1:end-3)];
        dataTable.X_MA_Lagged4(dataTable.ID==i) = [NaN(4,1); movavg(1:end-4)];
        dataTable.X_MA_Lagged5(dataTable.ID==i) = NaN(5,1);
        dataTable.X_MA_Lagged6(dataTable.ID==i) = NaN(5,1);
        dataTable.X_MA_Lagged7(dataTable.ID==i) = NaN(5,1);
        dataTable.X_MA_Lagged8(dataTable.ID==i) = NaN(5,1);
        dataTable.X_MA_Lagged9(dataTable.ID==i) = NaN(5,1);
    elseif length(movavg)>3
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = [NaN(2,1); movavg(1:end-2)];
        dataTable.X_MA_Lagged3(dataTable.ID==i) = [NaN(3,1); movavg(1:end-3)];
        dataTable.X_MA_Lagged4(dataTable.ID==i) = NaN(4,1);
        dataTable.X_MA_Lagged5(dataTable.ID==i) = NaN(4,1);
        dataTable.X_MA_Lagged6(dataTable.ID==i) = NaN(4,1);
        dataTable.X_MA_Lagged7(dataTable.ID==i) = NaN(4,1);
        dataTable.X_MA_Lagged8(dataTable.ID==i) = NaN(4,1);
        dataTable.X_MA_Lagged9(dataTable.ID==i) = NaN(4,1);
    elseif length(movavg)>2
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = [NaN(2,1); movavg(1:end-2)];
        dataTable.X_MA_Lagged3(dataTable.ID==i) = NaN(3,1);
        dataTable.X_MA_Lagged4(dataTable.ID==i) = NaN(3,1);
        dataTable.X_MA_Lagged5(dataTable.ID==i) = NaN(3,1);
        dataTable.X_MA_Lagged6(dataTable.ID==i) = NaN(3,1);
        dataTable.X_MA_Lagged7(dataTable.ID==i) = NaN(3,1);
        dataTable.X_MA_Lagged8(dataTable.ID==i) = NaN(3,1);
        dataTable.X_MA_Lagged9(dataTable.ID==i) = NaN(3,1);
    elseif length(movavg)>1
        dataTable.X_MA_Lagged1(dataTable.ID==i) = [NaN; movavg(1:end-1)];
        dataTable.X_MA_Lagged2(dataTable.ID==i) = NaN(2,1);
        dataTable.X_MA_Lagged3(dataTable.ID==i) = NaN(2,1);
        dataTable.X_MA_Lagged4(dataTable.ID==i) = NaN(2,1);
        dataTable.X_MA_Lagged5(dataTable.ID==i) = NaN(2,1);
        dataTable.X_MA_Lagged6(dataTable.ID==i) = NaN(2,1);
        dataTable.X_MA_Lagged7(dataTable.ID==i) = NaN(2,1);
        dataTable.X_MA_Lagged8(dataTable.ID==i) = NaN(2,1);
        dataTable.X_MA_Lagged9(dataTable.ID==i) = NaN(2,1);
    end
end

%% Load and modify survival probabilities from R

survivalprobabilities3(48,8:18) =  survivalprobabilities3(48,5:15);
survivalprobabilities3.time_14(48) = 1;
survivalprobabilities3.time_16(48) = 1;
survivalprobabilities3.time_18(48) = 1;

%%
Survival = survivalprobabilities3;

weeks = {};
for i = 6:2:106
    weeks = [weeks; ['Week',num2str(i)]]
end
Survival = removevars(Survival, "VarName1");
Survival = renamevars(Survival,Survival.Properties.VariableNames,['ID'; weeks]);
Survival.ID = IDs_BaselineUnmedicated_L1trend;

%%
IPCW = Survival;
IPCW(:,2:end) = 1./IPCW(:,2:end);

%% Analysis of start medication group
IDs_BaselineUnmedicated = [IDs_StartMedication; IDs_AllUnmedicated];

IDs_StartMedication_L1trend = intersect(IDs_StartMedication,IDs_BaselineUnmedicated_L1trend,'stable');
tremor_time_before_start = trend_tremor_time_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_StartMedication_L1trend),:);
tremor_time_after_start = trend_tremor_time_medicated_filled(ismember([IDs_BaselineMedicated; IDs_StartMedication] ,IDs_StartMedication_L1trend),:);
median_tremor_power_before_start = trend_median_tremor_power_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_StartMedication_L1trend),:);
median_tremor_power_after_start = trend_median_tremor_power_medicated_filled(ismember([IDs_BaselineMedicated; IDs_StartMedication] ,IDs_StartMedication_L1trend),:);
modal_tremor_power_before_start = trend_modal_tremor_power_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_StartMedication_L1trend),:);
modal_tremor_power_after_start = trend_modal_tremor_power_medicated_filled(ismember([IDs_BaselineMedicated; IDs_StartMedication] ,IDs_StartMedication_L1trend),:);
perc90_tremor_power_before_start = trend_perc90_tremor_power_unmedicated_filled(ismember(IDs_BaselineUnmedicated,IDs_StartMedication_L1trend),:);
perc90_tremor_power_after_start = trend_perc90_tremor_power_medicated_filled(ismember([IDs_BaselineMedicated; IDs_StartMedication] ,IDs_StartMedication_L1trend),:);

include_idx = find(~all(isnan(tremor_time_after_start),2));

IDs_StartMedication_L1trend_include = IDs_StartMedication_L1trend(include_idx);
StartWeek_include = StartWeek(ismember(IDs_BaselineUnmedicated_L1trend,IDs_StartMedication_L1trend_include));
StartWeek_include(mod(StartWeek_include,2)>0) = StartWeek_include(mod(StartWeek_include,2)>0) + 1;

tremor_time_before_start = tremor_time_before_start(include_idx,:);
tremor_time_after_start = tremor_time_after_start(include_idx,:);
median_tremor_power_before_start = median_tremor_power_before_start(include_idx,:);
median_tremor_power_after_start = median_tremor_power_after_start(include_idx,:);
modal_tremor_power_before_start = modal_tremor_power_before_start(include_idx,:);
modal_tremor_power_after_start = modal_tremor_power_after_start(include_idx,:);
perc90_tremor_power_before_start = perc90_tremor_power_before_start(include_idx,:);
perc90_tremor_power_after_start = perc90_tremor_power_after_start(include_idx,:);

week_vector = 0:2:104;

tremor_time_before = [];
tremor_time_after = [];
median_tremor_power_before = [];
median_tremor_power_after = [];
modal_tremor_power_before = [];
modal_tremor_power_after = [];
perc90_tremor_power_before = [];
perc90_tremor_power_after = [];


IDs_StartMedication_L1trend_include(38) = []; % participant is using anticholinergic medication

for i = 1:length(IDs_StartMedication_L1trend_include)
    
    weeks_before = find(week_vector<StartWeek_include(i));
    tremor_time_before(i) = tremor_time_before_start(i,weeks_before(end));
    median_tremor_power_before(i) = median_tremor_power_before_start(i,weeks_before(end));
    modal_tremor_power_before(i) = modal_tremor_power_before_start(i,weeks_before(end));
    perc90_tremor_power_before(i) = perc90_tremor_power_before_start(i,weeks_before(end));

    weeks_after = find(week_vector>=StartWeek_include(i));
    tremor_time_after(i,1:8) = tremor_time_after_start(i,weeks_after(1:8));
    median_tremor_power_after(i,1:8) = median_tremor_power_after_start(i,weeks_after(1:8));
    modal_tremor_power_after(i,1:8) = modal_tremor_power_after_start(i,weeks_after(1:8));
    perc90_tremor_power_after(i,1:8) = perc90_tremor_power_after_start(i,weeks_after(1:8));

end

tremor_time_before = logit(tremor_time_before);
tremor_time_after = logit(tremor_time_after);

%% Calculate SRM
addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression'))
SRM_function = @(x)mean(x,'omitnan')/std(x,'omitnan');

SRM_tremor_time = [];
CI_SRM_tremor_time = [];
SRM_median_tremor_power = [];
CI_SRM_median_tremor_power = [];
SRM_modal_tremor_power = [];
CI_SRM_modal_tremor_power = [];
SRM_perc90_tremor_power = [];
CI_SRM_perc90_tremor_power = [];

for i = 1:8
    SRM_tremor_time(i) = mean(tremor_time_after(:,i) - tremor_time_before','omitnan')/(std(tremor_time_after(:,i) - tremor_time_before','omitnan'));
    CI_SRM_tremor_time(:,i) = bootci(1000,SRM_function,tremor_time_after(:,i) - tremor_time_before');
    N_tremor_time(i) = length(find(~isnan(tremor_time_after(:,i)) & ~isnan(tremor_time_before')));
    SRM_median_tremor_power(i) = mean(median_tremor_power_after(:,i) - median_tremor_power_before','omitnan')/(std(median_tremor_power_after(:,i) - median_tremor_power_before','omitnan'));
    CI_SRM_median_tremor_power(:,i) = bootci(1000,SRM_function,median_tremor_power_after(:,i) - median_tremor_power_before');
    N_median_tremor_power(i) = length(find(~isnan(median_tremor_power_after(:,i)) & ~isnan(median_tremor_power_before')));
    SRM_modal_tremor_power(i) = mean(modal_tremor_power_after(:,i) - modal_tremor_power_before','omitnan')/(std(modal_tremor_power_after(:,i) - modal_tremor_power_before','omitnan'));
    CI_SRM_modal_tremor_power(:,i) = bootci(1000,SRM_function,modal_tremor_power_after(:,i) - modal_tremor_power_before');
    N_modal_tremor_power(i) = length(find(~isnan(modal_tremor_power_after(:,i)) & ~isnan(modal_tremor_power_before')));
    SRM_perc90_tremor_power(i) = mean(perc90_tremor_power_after(:,i) - perc90_tremor_power_before','omitnan')/(std(perc90_tremor_power_after(:,i) - perc90_tremor_power_before','omitnan'));
    CI_SRM_perc90_tremor_power(:,i) = bootci(1000,SRM_function,perc90_tremor_power_after(:,i) - perc90_tremor_power_before');
    N_perc90_tremor_power(i) = length(find(~isnan(perc90_tremor_power_after(:,i)) & ~isnan(perc90_tremor_power_before')));
end

%% Plot delta's
close all;
figure(); hold on;
measure = tremor_time_after - tremor_time_before';
Xband = [0:2:14 14:-2:0];
plot(0:2:14,median(measure,'omitnan'),'.-k')
upperband = prctile(measure,75);
fill(Xband,[prctile(measure,25) upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
yline(0,'--')
ylabel('\Delta tremor time (log odds-ratio)')
xlabel('Weeks since start of dopaminergic treatment')
legend('median', 'IQR')
ylim([-0.6 0.3])

measure = median_tremor_power_after - median_tremor_power_before';
figure(); hold on;
Xband = [0:2:14 14:-2:0];
plot(0:2:14,median(measure,'omitnan'),'.-k')
upperband = prctile(measure,75);
fill(Xband,[prctile(measure,25) upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
yline(0,'--')
ylabel('\Delta median tremor power')
xlabel('Weeks since start of dopaminergic treatment')
legend('median', 'IQR')
ylim([-0.6 0.3])

measure = modal_tremor_power_after - modal_tremor_power_before';
figure(); hold on;
Xband = [0:2:14 14:-2:0];
plot(0:2:14,median(measure,'omitnan'),'.-k')
upperband = prctile(measure,75);
fill(Xband,[prctile(measure,25) upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
yline(0,'--')
ylabel('\Delta modal tremor power')
xlabel('Weeks since start of dopaminergic treatment')
legend('median', 'IQR')
ylim([-0.6 0.3])

measure = perc90_tremor_power_after - perc90_tremor_power_before';
figure(); hold on;
Xband = [0:2:14 14:-2:0];
plot(0:2:14,median(measure,'omitnan'),'.-k')
upperband = prctile(measure,75);
fill(Xband,[prctile(measure,25) upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
yline(0,'--')
ylabel('\Delta 90th percentile of tremor power')
xlabel('Weeks since start of dopaminergic treatment')
legend('median', 'IQR')
ylim([-0.6 0.3])


%% Plot SRM start medication
figure(); hold on;
errorbar(0:2:14,SRM_tremor_time,(CI_SRM_tremor_time(2,:)-CI_SRM_tremor_time(1,:))/2,'k')
xlabel('Weeks since start of dopaminergic treatment')
ylabel('SRM (with 95%-CI)')
xlim([ -1 15])
yline(0,'--')
ylim([-1.2 1.2])
title('tremor time')

figure(); hold on;
errorbar(0:2:14,SRM_median_tremor_power,(CI_SRM_median_tremor_power(2,:)-CI_SRM_median_tremor_power(1,:))/2,'k')
xlabel('Weeks since start of dopaminergic treatment')
ylabel('SRM (with 95%-CI)')
xlim([ -1 15])
yline(0,'--')
ylim([-1.2 1.2])
title('median tremor power')

figure(); hold on;
errorbar(0:2:14,SRM_modal_tremor_power,(CI_SRM_modal_tremor_power(2,:)-CI_SRM_modal_tremor_power(1,:))/2,'k')
xlabel('Weeks since start of dopaminergic treatment')
ylabel('SRM (with 95%-CI)')
xlim([ -1 15])
yline(0,'--')
ylim([-1.2 1.2])
title('modal tremor power')

figure(); hold on;
errorbar(0:2:14,SRM_perc90_tremor_power,(CI_SRM_perc90_tremor_power(2,:)-CI_SRM_perc90_tremor_power(1,:))/2,'k')
xlabel('Weeks since start of dopaminergic treatment')
ylabel('SRM (with 95%-CI)')
xlim([ -1 15])
yline(0,'--')
ylim([-1.2 1.2])
title('90th percentile of tremor power')
