%%
clear all; close all;

%%
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250410_matching_info\matching_info.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250410_matching_info\IDs_BaselineMedicated_matched.mat')

BaselineMedicated = BaselineMedicated(ismember(BaselineMedicated.ID,IDs_BaselineMedicated_matched),:);
StartMedication = BaselineUnmedicated(ismember(BaselineUnmedicated.ID,IDs_StartMedication_L1trend),:);

%% Create survival curve
n = length(StartWeek);

% Sort weeks
sortedWeeks = sort(StartWeek);

% Compute number remaining at each unique week
[counts, uniqueWeeks] = groupcounts(sortedWeeks');

numRemaining = n - cumsum(counts);

% Plot
figure;
stairs([0; uniqueWeeks-1], [n; numRemaining],'k', 'LineWidth', 1.5);
xlabel('Weeks since baseline');
ylabel('Number of participants in unmedicated group');
ylim([0 80])
xlim([0 100])
%% Age
% median(BaselineMedicated.Age)
% prctile(BaselineMedicated.Age,25)
% prctile(BaselineMedicated.Age,75)

% median(BaselineUnmedicated.Age)
% prctile(BaselineUnmedicated.Age,25)
% prctile(BaselineUnmedicated.Age,75)

median(StartMedication.Age)
prctile(StartMedication.Age,25)
prctile(StartMedication.Age,75)

% [h,p] = ttest2(BaselineMedicated.Age,BaselineUnmedicated.Age)

%% Disease duration
% median(BaselineMedicated.Disease_duration)
% prctile(BaselineMedicated.Disease_duration,25)
% prctile(BaselineMedicated.Disease_duration,75)
% 
% median(BaselineUnmedicated.Disease_duration,'omitnan')
% prctile(BaselineUnmedicated.Disease_duration,25)
% prctile(BaselineUnmedicated.Disease_duration,75)

median(StartMedication.Disease_duration,'omitnan')
prctile(StartMedication.Disease_duration,25)
prctile(StartMedication.Disease_duration,75)
% 
% p = ranksum(BaselineMedicated.Disease_duration,BaselineUnmedicated.Disease_duration)

%% Most-affected
n1 = length(find(BaselineMedicated.Most_affected==1))
n2 = length(find(BaselineMedicated.Most_affected==0))
n3 = length(find(StartMedication.Most_affected==1))
n4 = length(find(StartMedication.Most_affected==0))

most_affected = [1 * ones(n1,1); 2 * ones(n2,1); 1 * ones(n3,1); 2 * ones(n4,1)]; % 1 = Male, 2 = Female
group  = [1 * ones(77,1); 2 * ones(76,1)]; % 1 = Group 1, 2 = Group 2
[table, chi2, p, labels] = crosstab(most_affected, group)

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
%%
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\IDs.mat')

%% Age
% BaselineMedicated_gender = Visit1PPP.Gender(ismember(Visit1PPP.id,IDs_BaselineMedicated_matched));
BaselineUnmedicated_gender = Visit1PPP.Gender(ismember(Visit1PPP.id,IDs_StartMedication_L1trend));
BaselineUnmedicated_gender = [BaselineUnmedicated_gender; Visit1DeNovo.Gender(ismember(Visit1DeNovo.id,IDs_StartMedication_L1trend))];

% n1 = length(find(BaselineMedicated_gender==1))
% n2 = length(find(BaselineMedicated_gender==2))
n3 = length(find(BaselineUnmedicated_gender==1))
n4 = length(find(BaselineUnmedicated_gender==2))

gender = [1 * ones(n1,1); 2 * ones(n2,1); 1 * ones(n3,1); 2 * ones(n4,1)]; % 1 = Male, 2 = Female
group  = [1 * ones(77,1); 2 * ones(77,1)]; % 1 = Group 1, 2 = Group 2
[table, chi2, p, labels] = crosstab(gender, group)

%% UPDRS 1
Visit1 = Visit1DeNovo;
UPDRS1a = [Visit1.Up1aCognit, Visit1.Up1aHalPsy, Visit1.Up1aDepres, Visit1.Up1aAnxious, Visit1.Up1aApathy, Visit1.Up1aDopDysSyn];
UPDRS1a_visit1 = sum(UPDRS1a,2);

UPDRS1b = [Visit1.Updrs2It07,Visit1.Updrs2It08,Visit1.Updrs2It09,Visit1.Updrs2It10,Visit1.Updrs2It11,Visit1.Updrs2It12,Visit1.Updrs2It13];
UPDRS1b_visit1 = sum(UPDRS1b,2);

UPDRS1_visit1 = UPDRS1a_visit1 + UPDRS1b_visit1;

UPDRS1a_PPP = [Visit1PPP.Up1aCognit, Visit1PPP.Up1aHalPsy, Visit1PPP.Up1aDepres, Visit1PPP.Up1aAnxious, Visit1PPP.Up1aApathy, Visit1PPP.Up1aDopDysSyn];
UPDRS1a_visit1_PPP = sum(UPDRS1a_PPP,2);

UPDRS1b_PPP = [Visit1PPP.Updrs2It07,Visit1PPP.Updrs2It08,Visit1PPP.Updrs2It09,Visit1PPP.Updrs2It10,Visit1PPP.Updrs2It11,Visit1PPP.Updrs2It12,Visit1PPP.Updrs2It13];
UPDRS1b_visit1_PPP = sum(UPDRS1b_PPP,2);

UPDRS1_visit1_PPP = UPDRS1a_visit1_PPP + UPDRS1b_visit1_PPP;

% BaselineMedicated_UPDRS1 = UPDRS1_visit1_PPP(ismember(Visit1PPP.id,IDs_BaselineMedicated_matched));
BaselineUnmedicated_UPDRS1 = UPDRS1_visit1_PPP(ismember(Visit1PPP.id,IDs_StartMedication_L1trend));
BaselineUnmedicated_UPDRS1 = [BaselineUnmedicated_UPDRS1; UPDRS1_visit1(ismember(Visit1.id,IDs_StartMedication_L1trend))];

% median(BaselineMedicated_UPDRS1,'omitnan')
% prctile(BaselineMedicated_UPDRS1,25)
% prctile(BaselineMedicated_UPDRS1,75)

median(BaselineUnmedicated_UPDRS1,'omitnan')
prctile(BaselineUnmedicated_UPDRS1,25)
prctile(BaselineUnmedicated_UPDRS1,75)

% p = ranksum(BaselineMedicated_UPDRS1,BaselineUnmedicated_UPDRS1)

% length(find(isnan(BaselineUnmedicated_UPDRS1)))

%% UPDRS 2
UPDRS2 = [Visit1.Updrs2It14,Visit1.Updrs2It15,Visit1.Updrs2It16,Visit1.Updrs2It17,Visit1.Updrs2It18,Visit1.Updrs2It19,Visit1.Updrs2It20,...
    Visit1.Updrs2It21,Visit1.Updrs2It22,Visit1.Updrs2It23,Visit1.Updrs2It24,Visit1.Updrs2It25,Visit1.Updrs2It26];
UPDRS2_visit1 = sum(UPDRS2,2);

UPDRS2_PPP = [Visit1PPP.Updrs2It14,Visit1PPP.Updrs2It15,Visit1PPP.Updrs2It16,Visit1PPP.Updrs2It17,Visit1PPP.Updrs2It18,Visit1PPP.Updrs2It19,Visit1PPP.Updrs2It20,...
    Visit1PPP.Updrs2It21,Visit1PPP.Updrs2It22,Visit1PPP.Updrs2It23,Visit1PPP.Updrs2It24,Visit1PPP.Updrs2It25,Visit1PPP.Updrs2It26];
UPDRS2_visit1_PPP = sum(UPDRS2_PPP,2);

% BaselineMedicated_UPDRS2 = UPDRS2_visit1_PPP(ismember(Visit1PPP.id,IDs_BaselineMedicated_matched));
BaselineUnmedicated_UPDRS2 = UPDRS2_visit1_PPP(ismember(Visit1PPP.id,IDs_StartMedication_L1trend));
BaselineUnmedicated_UPDRS2 = [BaselineUnmedicated_UPDRS2; UPDRS2_visit1(ismember(Visit1.id,IDs_StartMedication_L1trend))];
% 
% median(BaselineMedicated_UPDRS2,'omitnan')
% prctile(BaselineMedicated_UPDRS2,25)
% prctile(BaselineMedicated_UPDRS2,75)

median(BaselineUnmedicated_UPDRS2,'omitnan')
prctile(BaselineUnmedicated_UPDRS2,25)
prctile(BaselineUnmedicated_UPDRS2,75)

% p = ranksum(BaselineMedicated_UPDRS2,BaselineUnmedicated_UPDRS2)

length(find(isnan(BaselineUnmedicated_UPDRS2)))

%% UPDRS 3
clc;

UPDRS3OFF = [Visit1.Up3OfSpeech, Visit1.Up3OfFacial, Visit1.Up3OfRigNec, Visit1.Up3OfRigRue, Visit1.Up3OfRigLue, Visit1.Up3OfRigRle,...
    Visit1.Up3OfRigLle, Visit1.Up3OfFiTaNonDev, Visit1.Up3OfFiTaYesDev, Visit1.Up3OfHaMoNonDev, Visit1.Up3OfHaMoYesDev, Visit1.Up3OfProSNonDev,...
    Visit1.Up3OfProSYesDev, Visit1.Up3OfToTaNonDev, Visit1.Up3OfToTaYesDev, Visit1.Up3OfLAgiNonDev, Visit1.Up3OfLAgiYesDev, Visit1.Up3OfArise,...
    Visit1.Up3OfGait, Visit1.Up3OfFreez, Visit1.Up3OfStaPos, Visit1.Up3OfPostur, Visit1.Up3OfSpont, Visit1.Up3OfPosTNonDev, Visit1.Up3OfPosTYesDev,...
    Visit1.Up3OfKinTreNonDev, Visit1.Up3OfKinTreYesDev, Visit1.Up3OfRAmpArmNonDev, Visit1.Up3OfRAmpArmYesDev, Visit1.Up3OfRAmpLegNonDev,...
    Visit1.Up3OfRAmpLegYesDev, Visit1.Up3OfRAmpJaw, Visit1.Up3OfConstan];
UPDRS3OFF_visit1 = sum(UPDRS3OFF,2);

UPDRS3OFF_PPP = [Visit1PPP.Up3OfSpeech, Visit1PPP.Up3OfFacial, Visit1PPP.Up3OfRigNec, Visit1PPP.Up3OfRigRue, Visit1PPP.Up3OfRigLue, Visit1PPP.Up3OfRigRle,...
    Visit1PPP.Up3OfRigLle, Visit1PPP.Up3OfFiTaNonDev, Visit1PPP.Up3OfFiTaYesDev, Visit1PPP.Up3OfHaMoNonDev, Visit1PPP.Up3OfHaMoYesDev, Visit1PPP.Up3OfProSNonDev,...
    Visit1PPP.Up3OfProSYesDev, Visit1PPP.Up3OfToTaNonDev, Visit1PPP.Up3OfToTaYesDev, Visit1PPP.Up3OfLAgiNonDev, Visit1PPP.Up3OfLAgiYesDev, Visit1PPP.Up3OfArise,...
    Visit1PPP.Up3OfGait, Visit1PPP.Up3OfFreez, Visit1PPP.Up3OfStaPos, Visit1PPP.Up3OfPostur, Visit1PPP.Up3OfSpont, Visit1PPP.Up3OfPosTNonDev, Visit1PPP.Up3OfPosTYesDev,...
    Visit1PPP.Up3OfKinTreNonDev, Visit1PPP.Up3OfKinTreYesDev, Visit1PPP.Up3OfRAmpArmNonDev, Visit1PPP.Up3OfRAmpArmYesDev, Visit1PPP.Up3OfRAmpLegNonDev,...
    Visit1PPP.Up3OfRAmpLegYesDev, Visit1PPP.Up3OfRAmpJaw, Visit1PPP.Up3OfConstan];
UPDRS3OFF_visit1_PPP = sum(UPDRS3OFF_PPP,2);

% BaselineMedicated_UPDRS3 = UPDRS3OFF_visit1_PPP(ismember(Visit1PPP.id,IDs_BaselineMedicated_matched));
BaselineUnmedicated_UPDRS3 = UPDRS3OFF_visit1_PPP(ismember(Visit1PPP.id,IDs_StartMedication_L1trend));
BaselineUnmedicated_UPDRS3 = [BaselineUnmedicated_UPDRS3; UPDRS3OFF_visit1(ismember(Visit1.id,IDs_StartMedication_L1trend))];

% median(BaselineMedicated_UPDRS3,'omitnan')
% prctile(BaselineMedicated_UPDRS3,25)
% prctile(BaselineMedicated_UPDRS3,75)

median(BaselineUnmedicated_UPDRS3,'omitnan')
prctile(BaselineUnmedicated_UPDRS3,25)
prctile(BaselineUnmedicated_UPDRS3,75)

% p = ranksum(BaselineMedicated_UPDRS3,BaselineUnmedicated_UPDRS3)

length(find(isnan(BaselineUnmedicated_UPDRS3)))

%% UPDRS 4
clc;
UPDRS4_visit1 = Visit1.MotComTotal;
UPDRS4_visit1_PPP = Visit1PPP.MotComTotal;

% BaselineMedicated_UPDRS4 = UPDRS4_visit1_PPP(ismember(Visit1PPP.id,IDs_BaselineMedicated_matched));
BaselineUnmedicated_UPDRS4 = UPDRS4_visit1_PPP(ismember(Visit1PPP.id,IDs_StartMedication_L1trend));
BaselineUnmedicated_UPDRS4 = [BaselineUnmedicated_UPDRS4; UPDRS4_visit1(ismember(Visit1.id,IDs_StartMedication_L1trend))];

% median(BaselineMedicated_UPDRS4,'omitnan')
% prctile(BaselineMedicated_UPDRS4,25)
% prctile(BaselineMedicated_UPDRS4,75)

median(BaselineUnmedicated_UPDRS4,'omitnan')
prctile(BaselineUnmedicated_UPDRS4,25)
prctile(BaselineUnmedicated_UPDRS4,75)

% p = ranksum(BaselineMedicated_UPDRS4,BaselineUnmedicated_UPDRS4)

length(find(isnan(BaselineUnmedicated_UPDRS4)))

%% UPDRS 3 tremor subscore
clc;

UPDRS3OFF = [Visit1.Up3OfPosTNonDev, Visit1.Up3OfPosTYesDev, Visit1.Up3OfKinTreNonDev, Visit1.Up3OfKinTreYesDev,...
    Visit1.Up3OfRAmpArmNonDev, Visit1.Up3OfRAmpArmYesDev, Visit1.Up3OfRAmpLegNonDev,...
    Visit1.Up3OfRAmpLegYesDev, Visit1.Up3OfRAmpJaw, Visit1.Up3OfConstan];
UPDRS3OFF_visit1 = sum(UPDRS3OFF,2);

UPDRS3OFF_PPP = [Visit1PPP.Up3OfPosTNonDev, Visit1PPP.Up3OfPosTYesDev, Visit1PPP.Up3OfKinTreNonDev, Visit1PPP.Up3OfKinTreYesDev, ...
    Visit1PPP.Up3OfRAmpArmNonDev, Visit1PPP.Up3OfRAmpArmYesDev, Visit1PPP.Up3OfRAmpLegNonDev,...
    Visit1PPP.Up3OfRAmpLegYesDev, Visit1PPP.Up3OfRAmpJaw, Visit1PPP.Up3OfConstan];
UPDRS3OFF_visit1_PPP = sum(UPDRS3OFF_PPP,2);

% BaselineMedicated_UPDRS3 = UPDRS3OFF_visit1_PPP(ismember(Visit1PPP.id,IDs_BaselineMedicated_matched));
BaselineUnmedicated_UPDRS3 = UPDRS3OFF_visit1_PPP(ismember(Visit1PPP.id,IDs_StartMedication_L1trend));
BaselineUnmedicated_UPDRS3 = [BaselineUnmedicated_UPDRS3; UPDRS3OFF_visit1(ismember(Visit1.id,IDs_StartMedication_L1trend))];

% median(BaselineMedicated_UPDRS3,'omitnan')
% prctile(BaselineMedicated_UPDRS3,25)
% prctile(BaselineMedicated_UPDRS3,75)

median(BaselineUnmedicated_UPDRS3,'omitnan')
prctile(BaselineUnmedicated_UPDRS3,25)
prctile(BaselineUnmedicated_UPDRS3,75)

% p = ranksum(BaselineMedicated_UPDRS3,BaselineUnmedicated_UPDRS3)

length(find(isnan(BaselineUnmedicated_UPDRS3)))

%% Rest tremor device sided arm
clc;
Rest_tremor_DeNovo = Visit1.Up3OfRAmpArmYesDev;
Rest_tremor_PPP = Visit1PPP.Up3OfRAmpArmYesDev;

% BaselineMedicated_Rest_tremor = Rest_tremor_PPP(ismember(Visit1PPP.id,IDs_BaselineMedicated_matched));
BaselineUnmedicated_Rest_tremor = Rest_tremor_PPP(ismember(Visit1PPP.id,IDs_StartMedication_L1trend));
BaselineUnmedicated_Rest_tremor = [BaselineUnmedicated_Rest_tremor; Rest_tremor_DeNovo(ismember(Visit1.id,IDs_StartMedication_L1trend))];
% 
% median(BaselineMedicated_Rest_tremor,'omitnan')
% prctile(BaselineMedicated_Rest_tremor,25)
% prctile(BaselineMedicated_Rest_tremor,75)

median(BaselineUnmedicated_Rest_tremor,'omitnan')
prctile(BaselineUnmedicated_Rest_tremor,25)
prctile(BaselineUnmedicated_Rest_tremor,75)

% p = ranksum(BaselineMedicated_Rest_tremor,BaselineUnmedicated_Rest_tremor)

length(find(isnan(BaselineUnmedicated_Rest_tremor)))

%% Check visit date
% Determine visit closest to week 100
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\visit_week_numbers.mat')

k = size(visit_week_numbers,1);
visit2_week_diff = NaN(1,k);
visit3_week_diff = NaN(1,k);
for i = 1:k
    if ~isempty(visit_week_numbers.visit2{i})
        visit2_week_diff(i) = visit_week_numbers.visit2{i}-104;
    end
    if ~isempty(visit_week_numbers.visit3{i})
        visit3_week_diff(i) = visit_week_numbers.visit3{i}-104;
    end
end 
visit2_week_diff(find(-visit2_week_diff<visit3_week_diff))
visit3_week_diff(find(-visit2_week_diff<visit3_week_diff))
visit_week_numbers.ID(find(-visit2_week_diff<visit3_week_diff))
