%% 
clear all; close all;
%%
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\IDs.mat')

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
%% Split sensor data in medicated and unmedicated group
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\Inclusion.mat');

tremor_time_medicated = tremor_time(ismember(IDs_inclusion,[IDs_BaselineMedicated; IDs_StartMedication]),:);
tremor_time_unmedicated = tremor_time(ismember(IDs_inclusion,[IDs_StartMedication; IDs_AllUnmedicated]),:);
median_tremor_power_medicated = median_tremor_power(ismember(IDs_inclusion,[IDs_BaselineMedicated; IDs_StartMedication]),:); 
median_tremor_power_unmedicated = median_tremor_power(ismember(IDs_inclusion,[IDs_StartMedication; IDs_AllUnmedicated]),:);
modal_tremor_power_medicated = modal_tremor_power(ismember(IDs_inclusion,[IDs_BaselineMedicated; IDs_StartMedication]),:); 
modal_tremor_power_unmedicated = modal_tremor_power(ismember(IDs_inclusion,[IDs_StartMedication; IDs_AllUnmedicated]),:);
perc90_tremor_power_medicated = perc90_tremor_power(ismember(IDs_inclusion,[IDs_BaselineMedicated; IDs_StartMedication]),:); 
perc90_tremor_power_unmedicated = perc90_tremor_power(ismember(IDs_inclusion,[IDs_StartMedication; IDs_AllUnmedicated]),:);

start_week = Inclusion.StartWeek(ismember(Inclusion.ID,IDs_StartMedication));
start_week(mod(start_week,2)>0) = start_week(mod(start_week,2)>0) + 1; % Change odd start week number to even start week number

for i = 1:length(IDs_StartMedication)
    
    tremor_time_medicated(i+length(IDs_BaselineMedicated),week_vector<start_week(i)) = NaN(1,length(find(week_vector<start_week(i))));
    median_tremor_power_medicated(i+length(IDs_BaselineMedicated),week_vector<start_week(i)) = NaN(1,length(find(week_vector<start_week(i))));
    modal_tremor_power_medicated(i+length(IDs_BaselineMedicated),week_vector<start_week(i)) = NaN(1,length(find(week_vector<start_week(i))));
    perc90_tremor_power_medicated(i+length(IDs_BaselineMedicated),week_vector<start_week(i)) = NaN(1,length(find(week_vector<start_week(i))));

    tremor_time_unmedicated(i,week_vector>=start_week(i)) = NaN(1,length(find(week_vector>=start_week(i))));
    median_tremor_power_unmedicated(i,week_vector>=start_week(i)) = NaN(1,length(find(week_vector>=start_week(i))));
    modal_tremor_power_unmedicated(i,week_vector>=start_week(i)) = NaN(1,length(find(week_vector>=start_week(i))));
    perc90_tremor_power_unmedicated(i,week_vector>=start_week(i)) = NaN(1,length(find(week_vector>=start_week(i))));

end

%% Determine lambda per measure (concatenate medicated and unmedicated groups)

addpath(genpath('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\l1tf'))

measure = [tremor_time_medicated; tremor_time_unmedicated]; % Change to correct measure
lambdarng = logspace(-3,1,500); % Range of possible regularization constants
N = length(measure);
Etst = NaN(N,length(lambdarng));

%%
for i = 1:N
    tremor = measure(i,~isnan(measure(i,:)))';
    weeks = week_vector(~isnan(measure(i,:)))';
    if length(tremor)>=8
        Etst(i,:) = l1lscv(weeks, tremor, lambdarng);
    else
        Etst(i,:) = NaN(1,length(lambdarng));
    end
end

%% Find optimal lambda based on weighted error curve

number_of_weeks_available = [];
for i = 1:length(measure)
    number_of_weeks_available(i) = length(find(~isnan(measure(i,:))));
end
weights_unmedicated = number_of_weeks_available/mean(number_of_weeks_available);

mean_Etst = mean(Etst.*weights_unmedicated','omitnan'); % Calculate the mean error function
[~,j] = min(mean_Etst); % Determine the minimum
lambda = lambdarng(j) % Optimal lambda

figure();
plot(lambdarng, mean_Etst)

%% Calculate SNR (splitted)

measure = tremor_time_medicated; % Modify to run for the other measures
trend_tremor_time_medicated = NaN(length(measure),length(week_vector));
var_signal_medicated = [];
var_noise_medicated = [];
SNR_medicated = [];

for k = 1:length(measure)
    tremor = measure(k,~isnan(measure(k,:)))';
    weeks = week_vector(~isnan(measure(k,:)))';

    if length(tremor)>=8
        uhat = l1lsregr(weeks,tremor,lambda);
        trend_tremor_time_medicated(k,~isnan(measure(k,:))) = uhat;
        error = tremor-uhat;
        var_signal_medicated(k) = var(uhat,'omitnan');
        var_noise_medicated(k) = var(error,'omitnan');
        SNR_medicated(k) = var_signal_medicated(k)/var_noise_medicated(k);
    else
        var_signal_medicated(k) =  NaN;
        var_noise_medicated(k) =  NaN;
        SNR_medicated(k) = NaN;
    end
end

measure = tremor_time_unmedicated; % Modify to run for the other measures
trend_tremor_time_unmedicated = NaN(length(measure),length(week_vector));
var_signal_unmedicated = [];
var_noise_unmedicated = [];
SNR_unmedicated = [];

for k = 1:length(measure)
    tremor = measure(k,~isnan(measure(k,:)))';
    weeks = week_vector(~isnan(measure(k,:)))';

    if length(tremor)>=8
        uhat = l1lsregr(weeks,tremor,lambda);
        trend_tremor_time_unmedicated(k,~isnan(measure(k,:))) = uhat;
        error = tremor-uhat;
        var_signal_unmedicated(k) = var(uhat,'omitnan');
        var_noise_unmedicated(k) = var(error,'omitnan');
        SNR_unmedicated(k) = var_signal_unmedicated(k)/var_noise_unmedicated(k);
    else
        var_signal_unmedicated(k) =  NaN;
        var_noise_unmedicated(k) =  NaN;
        SNR_unmedicated(k) = NaN;
    end
end
%% Determine median and IQR of SNR

% median_SNR = median(SNR_medicated,'omitnan')
% IQR_SNR = [prctile(SNR_medicated,25) prctile(SNR_medicated,75)]
median_var_signal = median(var_signal_medicated./(var_signal_medicated + var_noise_medicated),'omitnan')
IQR_var_signal = [prctile(var_signal_medicated./(var_signal_medicated + var_noise_medicated),25) prctile(var_signal_medicated./(var_signal_medicated+var_noise_medicated),75)]
% median_var_noise = median(var_noise_medicated,'omitnan')
% IQR_var_noise = [prctile(var_noise_medicated,25) prctile(var_noise_medicated,75)]

% median_SNR = median(SNR_unmedicated,'omitnan')
% IQR_SNR = [prctile(SNR_unmedicated,25) prctile(SNR_unmedicated,75)]
median_var_signal = median(var_signal_unmedicated./(var_signal_unmedicated + var_noise_unmedicated),'omitnan')
IQR_var_signal = [prctile(var_signal_unmedicated./(var_signal_unmedicated + var_noise_unmedicated),25) prctile(var_signal_unmedicated./(var_signal_unmedicated + var_noise_unmedicated),75)]
% median_var_noise = median(var_noise_unmedicated,'omitnan')
% IQR_var_noise = [prctile(var_noise_unmedicated,25) prctile(var_noise_unmedicated,75)]

%% Interpolation and extrapolation

data = trend_modal_tremor_power_unmedicated; % Modify to run for the other measures
% Get the size of the matrix
[rows, cols] = size(data);

% Loop through each row to interpolate and extrapolate
filledData = data; % Initialize filledData as a copy of the original
for row = 1:rows
    % Extract the row data
    rowData = data(row,:);
    
    % Find the indices of non-NaN values
    validIdx = find(~isnan(rowData));
    
    % Find the indices of NaN values
    missingIdx = find(isnan(rowData));
    
    % If there are valid values, interpolate
    if ~isempty(validIdx)
        % Perform interpolation with extrapolation
        interpValues = interp1(validIdx, rowData(validIdx), (1:cols)', 'linear','extrap');
        
        % Identify the first and last known value
        firstValidIdx = validIdx(1);
        lastValidIdx = validIdx(end);
        
        % Find the number of NaNs at the end
        startNaNs = missingIdx(missingIdx < firstValidIdx);
        endNaNs = missingIdx(missingIdx > lastValidIdx);
        
        % Limit extrapolation to the first two NaNs
        if ~isempty(startNaNs)
            numStartNaNs = length(startNaNs);
            if numStartNaNs > 2
                interpValues(startNaNs(1:startNaNs(end)-2)) = NaN(numStartNaNs-2,1);
            end
        end
        % Limit extrapolation to the last two NaNs
        if ~isempty(endNaNs)
            numEndNaNs = length(endNaNs);
            if numEndNaNs > 2
                interpValues(endNaNs(3:end)) = NaN(numEndNaNs-2,1);
            end
        end      
    
        % Update the filledData
        filledData(row,:) = interpValues;
    end
end

%% Correct wrongly extrapolated trends after start of medication

week_vector = 0:2:104;

for i = 1:length(IDs_StartMedication)

    % trend_tremor_time_medicated_filled(i+length(IDs_BaselineMedicated),week_vector<start_week(i)) = NaN(1,length(find(week_vector<start_week(i))));
    % trend_median_tremor_power_medicated_filled(i+length(IDs_BaselineMedicated),week_vector<start_week(i)) = NaN(1,length(find(week_vector<start_week(i))));
    trend_modal_tremor_power_medicated_filled(i+length(IDs_BaselineMedicated),week_vector<start_week(i)) = NaN(1,length(find(week_vector<start_week(i))));
    % trend_perc90_tremor_power_medicated_filled(i+length(IDs_BaselineMedicated),week_vector<start_week(i)) = NaN(1,length(find(week_vector<start_week(i))));

    % trend_tremor_time_unmedicated_filled(i,week_vector>=start_week(i)) = NaN(1,length(find(week_vector>=start_week(i))));
    % trend_median_tremor_power_unmedicated_filled(i,week_vector>=start_week(i)) = NaN(1,length(find(week_vector>=start_week(i))));
    trend_modal_tremor_power_unmedicated_filled(i,week_vector>=start_week(i)) = NaN(1,length(find(week_vector>=start_week(i))));
    % trend_perc90_tremor_power_unmedicated_filled(i,week_vector>=start_week(i)) = NaN(1,length(find(week_vector>=start_week(i))));

end

%% Another correction (remove interpolated points if there was no data in the surrounding +- 4 datapoints)

measure = trend_modal_tremor_power_unmedicated; % Modify to run for the other measures

for i = 1:length(measure)
    isNan = isnan(measure(i,:));                   % Logical array: 1 where NaN, 0 elsewhere
    window = ones(1,5);                 % Look for 5 consecutive NaNs
    match = conv(double(isNan), window, 'valid') == 5;
    indices = find(match);
    remove_idx = indices+2; % Remove the middle point of the 5 consecutive NaNs
    trend_modal_tremor_power_unmedicated_filled(i,remove_idx) = NaN; % Modify
end

%% Correct negative values (set to NaN)
% 
% trend_tremor_time_medicated_filled(trend_tremor_time_medicated_filled<0) = NaN;
% trend_tremor_time_unmedicated_filled(trend_tremor_time_unmedicated_filled<0) = NaN;
% trend_median_tremor_power_medicated_filled(trend_median_tremor_power_medicated_filled<0) = NaN;
% trend_median_tremor_power_unmedicated_filled(trend_median_tremor_power_unmedicated_filled<0) = NaN;
trend_modal_tremor_power_medicated_filled(trend_modal_tremor_power_medicated_filled<0) = NaN;
trend_modal_tremor_power_unmedicated_filled(trend_modal_tremor_power_unmedicated_filled<0) = NaN;
% trend_perc90_tremor_power_medicated_filled(trend_perc90_tremor_power_medicated_filled<0) = NaN;
% trend_perc90_tremor_power_unmedicated_filled(trend_perc90_tremor_power_unmedicated_filled<0) = NaN;

%% Individual figures
week_vector = 0:2:104;
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_L1trend\trends2.mat')

close all;
C = colororder('glow12');
figure(); hold on;
scatter(week_vector,100*tremor_time_unmedicated(92,:),'k','filled')
plot(week_vector,100*trend_tremor_time_unmedicated_filled(92,:),'Color',C(2,:),'LineWidth',1.5)
xlabel('Weeks since baseline')
ylabel('Tremor time (% of inactive time)')
ylim([0 100])
xlim([0 104])

% Second figure
figure(); hold on;
% Create dummy scatter for 'Tremor time' legend entry
p1_dummy = scatter(nan, nan, 'k', 'filled');  % dummy for legend
p2 = scatter(week_vector, modal_tremor_power_unmedicated(92,:), 'kx');
p3 = scatter(week_vector, perc90_tremor_power_unmedicated(92,:), 'k^', 'filled');
plot(week_vector, trend_modal_tremor_power_unmedicated_filled(92,:), 'Color', C(2,:), 'LineWidth', 1.5)
p4 = plot(week_vector, trend_perc90_tremor_power_unmedicated_filled(92,:), 'Color', C(2,:), 'LineWidth', 1.5);

xlabel('Weeks since baseline')
ylabel('Tremor power (log values)')
legend([p1_dummy, p2, p3, p4], {'Tremor time', 'Modal tremor power', '90th percentile of tremor power', 'Piecewise linear trend'})
xlim([0 104])
ylim([0 4])

%%

close all;
C = colororder('glow12');
figure(); hold on;
% scatter(week_vector,100*tremor_time_unmedicated(73,:),'k','filled')
plot(week_vector,100*trend_tremor_time_unmedicated_filled(73,:),'Color',C(2,:),'LineWidth',1.5)
xlabel('Weeks since baseline')
ylabel('Tremor time (% of inactive time)')
ylim([0 100])
xlim([0 104])

% Second figure
figure(); hold on;
% Create dummy scatter for 'Tremor time' legend entry
p1_dummy = scatter(nan, nan, 'k', 'filled');  % dummy for legend
% p2 = scatter(week_vector, modal_tremor_power_unmedicated(73,:), 'kx');
% p3 = scatter(week_vector, perc90_tremor_power_unmedicated(73,:), 'k^', 'filled');
plot(week_vector, trend_modal_tremor_power_unmedicated_filled(73,:), 'Color', C(2,:), 'LineWidth', 1.5)
p4 = plot(week_vector, trend_perc90_tremor_power_unmedicated_filled(73,:), 'Color', C(2,:), 'LineWidth', 1.5);

xlabel('Weeks since baseline')
ylabel('Tremor power (log values)')
legend([p1_dummy, p2, p3, p4], {'Tremor time', 'Modal tremor power', '90th percentile of tremor power', 'Piecewise linear trend'})
xlim([0 104])
ylim([0 4])


%%
close all;
C = colororder('glow12');
figure(); hold on;
scatter(week_vector,100*tremor_time_medicated(116,:),'k','filled')
plot(week_vector,100*trend_tremor_time_medicated_filled(116,:),'Color',C(2,:),'LineWidth',1.5)
xlabel('Weeks since baseline')
ylabel('Tremor time (% of inactive time)')
ylim([0 100])
xlim([0 104])

% Second figure
figure(); hold on;
% Create dummy scatter for 'Tremor time' legend entry
p1_dummy = scatter(nan, nan, 'k', 'filled');  % dummy for legend
p2 = scatter(week_vector, modal_tremor_power_medicated(116,:), 'kx');
p3 = scatter(week_vector, perc90_tremor_power_medicated(116,:), 'k^', 'filled');
plot(week_vector, trend_modal_tremor_power_medicated_filled(116,:), 'Color', C(2,:), 'LineWidth', 1.5)
p4 = plot(week_vector, trend_perc90_tremor_power_medicated_filled(116,:), 'Color', C(2,:), 'LineWidth', 1.5);

xlabel('Weeks since baseline')
ylabel('Tremor power (log values)')
legend([p1_dummy, p2, p3, p4], {'Tremor time', 'Modal tremor power', '90th percentile of tremor power', 'Piecewise linear trend'})
xlim([0 104])
ylim([0 4])


%%
close all;
C = colororder('glow12');
figure(); hold on;
scatter(week_vector,100*tremor_time_medicated(180,:),'k','filled')
plot(week_vector,100*trend_tremor_time_medicated_filled(180,:),'Color',C(2,:),'LineWidth',1.5)
xlabel('Weeks since baseline')
ylabel('Tremor time (% of inactive time)')
ylim([0 100])
xlim([0 104])

% Second figure
figure(); hold on;
% Create dummy scatter for 'Tremor time' legend entry
p1_dummy = scatter(nan, nan, 'k', 'filled');  % dummy for legend
p2 = scatter(week_vector, modal_tremor_power_medicated(180,:), 'kx');
p3 = scatter(week_vector, perc90_tremor_power_medicated(180,:), 'k^', 'filled');
plot(week_vector, trend_modal_tremor_power_medicated_filled(180,:), 'Color', C(2,:), 'LineWidth', 1.5)
p4 = plot(week_vector, trend_perc90_tremor_power_medicated_filled(180,:), 'Color', C(2,:), 'LineWidth', 1.5);

xlabel('Weeks since baseline')
ylabel('Tremor power (log values)')
legend([p1_dummy, p2, p3, p4], {'Tremor time', 'Modal tremor power', '90th percentile of tremor power', 'Piecewise linear trend'})
xlim([0 104])
ylim([0 4])
