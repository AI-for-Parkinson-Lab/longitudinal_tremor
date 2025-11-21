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