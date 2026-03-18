%% Step 5: Create group level figures 
% Show the distribution of sensor-derived tremor measures, and the
% distribution of individual changes in these measures over time

clear all; close all;

%% Load sensor data and IDs
load('\\umcn.nl\nas\RBS\NEURO_AI4P\Users\Nienke Timmermans\Tremor progression\Derived_data\Trends_filled.mat')
load('\\umcn.nl\nas\RBS\NEURO_AI4P\Users\Nienke Timmermans\Tremor progression\Derived_data\IDs_selected.mat'); 

%% Select tremor time data above the false positive threshold
trend_tremor_time_unmedicated_filled_above_threshold = trend_tremor_time_unmedicated_filled;
trend_tremor_time_unmedicated_filled_above_threshold(isnan(trend_modal_tremor_power_unmedicated_filled))=NaN;
trend_tremor_time_medicated_filled_above_threshold = trend_tremor_time_medicated_filled;
trend_tremor_time_medicated_filled_above_threshold(isnan(trend_modal_tremor_power_medicated_filled))=NaN;

%% Plot distribution over time for medicated group

idx = contains([IDs_BaselineMedicated; IDs_BaselineUnmedicated],IDs_BaselineMedicated);

% Select measure (and change axes labels and limits accordingly)
measure = 100*trend_tremor_time_medicated_filled(idx,2:51);
% measure = trend_modal_tremor_power_medicated_filled(idx,2:51);
N = sum(~isnan(measure));

figure(); hold on;
Xband = [0:2:98 98:-2:0];
yyaxis left
plot(0:2:98,median(measure,'omitnan'),'k','LineWidth',2)
upperband = prctile(measure,75);
fill(Xband,[prctile(measure,25) upperband(end:-1:1)],'k','FaceAlpha',0.3,'EdgeColor','none')
upperband = prctile(measure,90);
fill(Xband,[prctile(measure,10) upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
xlabel('Weeks since baseline')
ylabel('Tremor time (% of inactive time)')
% ylabel('Modal percentile of tremor power (log values)')
ylim([0 100])
% ylim([0 4])
xlim([0 98])
yyaxis right
plot(0:2:98,N)
ylim([0 500])
ylabel('Number of participants')
% legend('median', '25th-75th percentiles','10th-90th percentiles')
ax = gca;
ax.YAxis(1).Color = 'k'; % Make left side black
ax.YAxis(2).Color = [0.8500 0.3250 0.0980]; % Make right side orange
title('Medicated tremor group')

%% Plot distribution over time for unmedicated group

addpath(genpath('utils'))

load('\\umcn.nl\nas\RBS\NEURO_AI4P\Users\Nienke Timmermans\Tremor progression\Derived_data\IPCW.mat')
weights = IPCW(:,2:end).Variables;
weights = [ones(78,3) weights]; % The weights start at week 6, so add 3 columns of ones
weights(:,51) = weights(:,50); % The weights don't change between week 98 and 100, duplicate to increase availability of weights

% Select measure (and change axes labels and limits accordingly)
measure = 100*trend_tremor_time_unmedicated_filled_above_threshold(:,2:51);
% measure = trend_modal_tremor_power_unmedicated_filled(:,2:51);

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
    
    percentiles = wprctile(outcome_t,[10, 25, 50, 75, 90],weights_t);
    
    p10(i) = percentiles(1);
    p25(i) = percentiles(2);
    medians(i) = percentiles(3);
    p75(i) = percentiles(4);
    p90(i) = percentiles(5);

end

N = sum(~isnan(measure));

figure(); hold on;
yyaxis left
Xband = [0:2:98 98:-2:0];
plot(0:2:98,medians,'k','LineWidth',2)
upperband = p75;
fill(Xband,[p25; upperband(end:-1:1)],'k','FaceAlpha',0.3,'EdgeColor','none')
upperband = p90;
fill(Xband,[p10; upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
xlabel('Weeks since baseline')
ylabel('Tremor time (% of inactive time)')
% ylabel('Modal tremor power (log values)')
ylim([0 100])
% ylim([0 4])
xlim([0 98])
yyaxis right
plot(0:2:98,N)
ylabel('Number of participants')
ylim([0 50])
ax = gca;
ax.YAxis(1).Color = 'k'; % Make left side black
ax.YAxis(2).Color = [0.8500 0.3250 0.0980]; % Make right side orange
title('Unmedicated tremor group')
% legend('median', '25th-75th percentiles','10th-90th percentiles')

%% Plot change over time for medicated group

idx = contains([IDs_BaselineMedicated; IDs_BaselineUnmedicated],IDs_BaselineMedicated);

% Select measure (and change axes labels and limits accordingly)
% measure = logit(trend_tremor_time_medicated_filled_above_threshold(idx,3:51)) - logit(trend_tremor_time_medicated_filled_above_threshold(idx,2));
measure = trend_modal_tremor_power_medicated_filled(idx,3:51) - trend_modal_tremor_power_medicated_filled(idx,2);
N = sum(~isnan(measure))

figure(); hold on;
Xband = [2:2:98 98:-2:2];
plot(2:2:98,median(measure,'omitnan'),'.-k')
upperband = prctile(measure,75);
fill(Xband,[prctile(measure,25) upperband(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none')
xlabel('Weeks since baseline')
yline(0,'--')
ylabel('\Delta modal tremor power (log values)')
% ylabel('\Delta tremor time (log odds-ratio)')
ylim([-0.5 1])
xlim([0 100])
legend('median', 'IQR')
title('Medicated tremor group')

%% Plot change over time for unmedicated group

% Select measure (and change axes labels and limits accordingly)
measure = trend_perc90_tremor_power_unmedicated_filled(:,3:end) -  trend_perc90_tremor_power_unmedicated_filled(:,2);
% measure = logit(trend_tremor_time_unmedicated_filled_above_threshold(:,3:end)) - logit(trend_tremor_time_unmedicated_filled_above_threshold(:,2));

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
% ylabel('\Delta tremor time (log odds-ratio)')
ylabel('\Delta 90th percentile of tremor power (log values)')
ylim([-0.5 1])
xlim([0 100])
legend('median', 'IQR')
title('Unmedicated tremor group')

%% Dominant side versus non-dominant side 

close all;
n1 = length(find(BaselineMedicated.Dominant_side==1));
n2 = length(find(BaselineMedicated.Dominant_side==0));
g1 = repmat({'Dominant side (n=38) '},n1,1);
g2 = repmat({'Non-dominant side (n=162)'},n2,1);
g = [g1; g2;];

x1 = trend_tremor_time_medicated_filled_above_threshold(BaselineMedicated.Dominant_side==1,2);
y1 = trend_tremor_time_medicated_filled_above_threshold(BaselineMedicated.Dominant_side==0,2);
length(find(~isnan(x1)))
length(find(~isnan(y1)))

figure(); hold on;
boxplot([100*x1;100*y1],g,'Colors','k','Symbol','.k','OutlierSize',10);
xlim([0.5 2.5])
ylabel('Tremor time (% of inactive time)')
title('Medicated tremor group')

n1 = length(find(BaselineUnmedicated.Dominant_side==1));
n2 = length(find(BaselineUnmedicated.Dominant_side==0));
g1 = repmat({'Dominant side (n=19) '},n1,1);
g2 = repmat({'Non-dominant side (n=28)'},n2,1);
g = [g1; g2;];

x2 = trend_tremor_time_unmedicated_filled_above_threshold(BaselineUnmedicated.Dominant_side==1,2);
y2 = trend_tremor_time_unmedicated_filled_above_threshold(BaselineUnmedicated.Dominant_side==0,2);
length(find(~isnan(x2)))
length(find(~isnan(y2)))

figure(); hold on;
boxplot([100*x2;100*y2],g,'Colors','k','Symbol','.k','OutlierSize',10);
xlim([0.5 2.5])
ylabel('Tremor time (% of inactive time)')
title('Unmedicated tremor group')


%%

p = ranksum(x1,y1) 

p = ranksum(x2,y2) 