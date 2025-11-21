%% Step 1: Participant selection
% Apply the following exclusion criteria:
% - alternative diagnosis
% - switch in watch side during the study period
% - fewer than eight weeks of sensor data 
% - starting symptomatic medication within the first eight weeks of data
% - having an unknown start date of symptomatic medication 

clear all; close all;

%% Load clinical data

% De Novo data:
Visit1DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit1_DeNovo.csv");
Visit2DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit2_DeNovo.csv");
Visit3DeNovo = readtable("C:\Users\z835211\Documents\Data\DeNovo\csv_files\Visit3_DeNovo.csv");

% PPP data:
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

% Create struct to store data
Inclusion = [];
Inclusion.ID = [Visit1PPP.id; Visit1DeNovo.id];
Inclusion = struct2table(Inclusion);
Inclusion.Group = [repmat({'PPP'},length(Visit1PPP.id),1); repmat({'DeNovo'},length(Visit1DeNovo.id),1)];

%% Exclude alternative diagnosis

NonPDsubjects = {'POMU04E638F5EC3A95C0','POMU065F753B3E97FF42','POMU10F33F607BAF7942','POMU4861269DC8E779A6','POMU584ADDDA24B571C1',...
    'POMU785C84D79EA74864','POMU7C7AE0CD7DBDF0E2','POMU8AAE6ADB5A1BE3C6','POMUA20B1A38DFC93AEB','POMUA2EA522320DEB7B4','POMUC27BD4A6AF046175'...
    'POMUCEA19060522828C9','POMUDC96BAA36834CD4C','POMUECAA325B2EB33E7B','POMUF2A9725453BD4C5B'}';
Inclusion.PDDiagnosis = ~ismember(Inclusion.ID,NonPDsubjects);

%% Exclude switch in watch side

for i = 1:length(Inclusion.ID)
    id = Inclusion.ID(i);
    if contains(Inclusion.Group(i),{'PPP'})
        Watchside_1(i) = Visit1PPP.WatchSide(contains(Visit1PPP.id,id));
        if ~isempty(find(contains(Visit2PPP.id,id)))
            Watchside_2(i) = Visit2PPP.WatchSide(contains(Visit2PPP.id,id));
        else
            Watchside_2(i) = NaN;
        end
        if ~isempty(find(contains(Visit3PPP.id,id)))
            Watchside_3(i) = Visit3PPP.WatchSide(contains(Visit3PPP.id,id));
        else
            Watchside_3(i) = NaN;
        end
    else
        Watchside_1(i) = Visit1DeNovo.WatchSide(contains(Visit1DeNovo.id,id));
        if ~isempty(find(contains(Visit2DeNovo.id,id)))
            Watchside_2(i) = Visit2DeNovo.WatchSide(contains(Visit2DeNovo.id,id));
        else
            Watchside_2(i) = NaN;
        end
        if ~isempty(find(contains(Visit3DeNovo.id,id)))
            Watchside_3(i) = Visit3DeNovo.WatchSide(contains(Visit3DeNovo.id,id));
        else
            Watchside_3(i) = NaN;
        end
    end
end

watchside_switch_12 = find(Watchside_1 ~= Watchside_2 & ~isnan(Watchside_1) & ~isnan(Watchside_2)); % Check if switch occurred between visit 1 and 2 
watchside_switch_23 = find(Watchside_2 ~= Watchside_3 & ~isnan(Watchside_2) & ~isnan(Watchside_3)); % Check if switch occurred between visit 2 and 3
watchside_switch_13 = find(Watchside_1 ~= Watchside_3 & ~isnan(Watchside_1) & ~isnan(Watchside_3)); % Check if switch occurred between visit 1 and 3 (if there was no visit 2 data)
watchside_switch = unique([watchside_switch_12 watchside_switch_13]); % find all unique indices

Inclusion.SameWatchside = ones(length(Inclusion.ID),1);
Inclusion.SameWatchside(watchside_switch) = 0;

%% Exclude <8 weeks with enough sensor data

% First load sensor data
addpath(genpath('utils\jsonlab'))
filename = 'Tremor_aggregates.json';
weeks = 0:2:104;
tremor_time = [];

for i = 1:length(Inclusion.ID)
    id = Inclusion.ID{i};
    for k = 1:length(weeks)
        week = weeks(k);
        if contains(Inclusion.Group(i),{'PPP'})
            file = ['C:\Users\z835211\Documents\Data\PPP\aggregated_output_191125\ppp\' num2str(week) '\' id '\' filename];
        else
            file = ['C:\Users\z835211\Documents\Data\DeNovo\aggregated_output_201125\denovo\' num2str(week) '\' id '\' filename];
        end
        if isfile(file)
            Tremor_aggregates = loadjson(file);
            if Tremor_aggregates.metadata.nr_valid_days >= 3 % Week should contain at least 3 valid days (with 10 hours of data)
                tremor_time(i,k) = Tremor_aggregates.aggregated_tremor_measures.perc_windows_tremor; 
            else
                tremor_time(i,k) = NaN;
            end
        else
            tremor_time(i,k) = NaN;
        end
    end
    % Check if there are at least 8 valid weeks of sensor data
        if length(find(~isnan(tremor_time(i,:))))>=8
            Inclusion.SensorData(i) = 1;
        else
            Inclusion.SensorData(i) = 0;
        end
end

%% Extract start week of symptomatic drug treatment

% First determine which participants already started before the study, and which participants do not start with drug treatment during the study period

% For PPP
for i = 1:size(Visit1PPP,1)
    id = Inclusion.ID(i);
    if Visit1PPP.ParkinMedUser(i)==1 % Participant has already started
        Inclusion.StartWeek(i) = -1;
    elseif Visit1PPP.ParkinMedUser(i)==0
        if ~isempty(Visit2PPP.ParkinMedUser(ismember(Visit2PPP.id,id)))
            if  Visit2PPP.ParkinMedUser(ismember(Visit2PPP.id,id))==0 && Visit3PPP.ParkinMedUser(ismember(Visit3PPP.id,id))==0
                Inclusion.StartWeek(i) = Inf;
            elseif isnan(Visit2PPP.ParkinMedUser(ismember(Visit2PPP.id,id))) && Visit3PPP.ParkinMedUser(ismember(Visit3PPP.id,id))==0
                Inclusion.StartWeek(i) = Inf;
            else
                Inclusion.StartWeek(i) = NaN;
            end
        elseif ~isempty(Visit3PPP.ParkinMedUser(ismember(Visit3PPP.id,id)))
            if Visit3PPP.ParkinMedUser(ismember(Visit3PPP.id,id))==0
                Inclusion.StartWeek(i) = Inf;
            else
                Inclusion.StartWeek(i) = NaN;
            end
        else
            Inclusion.StartWeek(i) = NaN;
        end

    else
        Inclusion.StartWeek(i) = NaN;
    end
end

% For De Novo:
for i = size(Visit1PPP,1)+1:size(Inclusion,1)
    id = Inclusion.ID(i);
    if ~isempty(Visit2DeNovo.Up3OfParkMedic(ismember(Visit2DeNovo.id,id))) && ~isempty(Visit3DeNovo.Up3OfParkMedic(ismember(Visit3DeNovo.id,id)))
        if Visit2DeNovo.Up3OfParkMedic(ismember(Visit2DeNovo.id,id))==1 || Visit3DeNovo.Up3OfParkMedic(ismember(Visit3DeNovo.id,id))==1
            Inclusion.StartWeek(i) = NaN;
        elseif Visit2DeNovo.Up3OfParkMedic(ismember(Visit2DeNovo.id,id))==0 && Visit3DeNovo.Up3OfParkMedic(ismember(Visit3DeNovo.id,id))==0
            Inclusion.StartWeek(i) = Inf;
        elseif isnan(Visit2DeNovo.Up3OfParkMedic(ismember(Visit2DeNovo.id,id))) && Visit3DeNovo.Up3OfParkMedic(ismember(Visit3DeNovo.id,id))==0
            Inclusion.StartWeek(i) = Inf;
        else
            Inclusion.StartWeek(i) = NaN;
        end
    elseif ~isempty(Visit3DeNovo.Up3OfParkMedic(ismember(Visit3DeNovo.id,id)))
        if Visit3DeNovo.Up3OfParkMedic(ismember(Visit3DeNovo.id,id))==0
            Inclusion.StartWeek(i) = Inf;
        elseif Visit3DeNovo.Up3OfParkMedic(ismember(Visit3DeNovo.id,id))==1
            Inclusion.StartWeek(i) = NaN;
        else
            Inclusion.StartWeek(i) = NaN;
        end
    else
        Inclusion.StartWeek(i) = NaN;
    end
end

% Then add start week

% For PPP:
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU048E610B74CC7A47')) = 24;
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU30163E78B5A9CCAB')) = 121; % starts 121 weeks after visit 1 (and visit 3 is in that week)
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU382082A2F77477DE')) = 14;
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU3CB991DD7E324AE6')) = 18;
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU42C333778EDD7E88')) = 4;
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU5A02B606856FA12F')) = 2; % starts with propranolol before study, but propranolol is not an exclusion criterium for De Novo
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU6A4763DC58628F76')) = 34;
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU8A79BAF84EA7AAA9')) = 61;
Inclusion.StartWeek(ismember(Inclusion.ID,'POMU947ACDE0DB8BD887')) = 103;

start_week_denovo = readtable('C:\Users\z835211\Documents\Data\DeNovo\start_week_denovo.csv');

for i = size(Visit1PPP,1)+1:size(Inclusion,1)
    id = Inclusion.ID(i);
    if ~isempty(find(ismember(start_week_denovo.ID,id)))
        Inclusion.StartWeek(i) = start_week_denovo.StartWeek(ismember(start_week_denovo.ID,id));
    end
end

%% Check if there are at least eight weeks of sensor data when unmedicated (necessary for piecewise linear trend estimation)

for i = 1:length(Inclusion.ID)
    if Inclusion.StartWeek(i)>=0
        start_week = Inclusion.StartWeek(i);
        if length(find(~isnan(tremor_time(i,weeks<start_week))))<8
            Inclusion.SensorData(i) = 0;
        end
    end
end

%% Divide participants into baseline medicated and unmedicated 

IDs_BaselineMedicated = Inclusion.ID(Inclusion.PDDiagnosis==1 & Inclusion.SameWatchside==1 & Inclusion.SensorData==1 & Inclusion.StartWeek==-1);
IDs_BaselineUnmedicated = Inclusion.ID(Inclusion.PDDiagnosis==1 & Inclusion.SameWatchside==1 & Inclusion.SensorData==1 & Inclusion.StartWeek>=0);

%% Save IDs and inclusion table

save('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IDs_selected.mat',"IDs_BaselineUnmedicated", "IDs_BaselineMedicated")
save('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Inclusion.mat',"Inclusion")
