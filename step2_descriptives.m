%% Step 2: Extract clinical and demographic characteristics of the baseline medicated and unmedicated groups
clear all; close all;

%% Load IDs
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IDs_selected.mat'); % load selected IDs

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

%% Extract demographics and clinical characteristics for medicated group

Age = [];
Gender = [];
Disease_duration = [];
Final_most_affected = [];
UPDRS1 = [];
UPDRS2 = [];
UPDRS3_OFF = [];
UPDRS4 = [];
UPDRS_tremor_OFF = [];
UPDRS_devside_tremor_OFF = [];

Visit1 = Visit1PPP;

for i = 1:length(IDs_BaselineMedicated)

    ppp_pep_number = IDs_BaselineMedicated{i};
    
    Age(i) = Visit1.Age(contains(Visit1.id,ppp_pep_number));
    Gender(i) = Visit1.Gender(contains(Visit1.id,ppp_pep_number));
    Disease_duration(i) = Visit1.MonthSinceDiag(contains(Visit1.id,ppp_pep_number));
   
    Tremor_score_dev_side(i) = Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number))+...
        Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number))+...
        Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number))+...
        Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number));
    Tremor_score_nondev_side(i) = Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number))+...
        Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number))+...
        Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number))+...
        Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number));
    
    Watchside(i) = Visit1.WatchSide(contains(Visit1.id,ppp_pep_number)); % side on which the watch is worn
    Most_affected_side(i) = Visit1.MostAffSide(contains(Visit1.id,ppp_pep_number));

    if Watchside(i) == 1 % right side
        
        Watchside_UPOf3 = [Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
            Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), ...
            Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number))];
        
        NonWatchside_UPOf3 = [Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)),...
            Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)), ...
            Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number))];

        if Most_affected_side(i) == 1 || Most_affected_side(i) == 3  % re>li affected
            Most_affected_subj(i)=1; % watch is worn on patient reported most-affected side
        elseif Most_affected_side(i) == 2 || Most_affected_side(i) == 4 % li>re affected
            Most_affected_subj(i)=0;  % watch is not worn on patient reported least-affected side
        elseif Most_affected_side(i) == 5 || Most_affected_side(i) == 6 % both sides equally affected
            Most_affected_subj(i)=2; % watch is worn on one of equally affected sides
        end

    elseif Watchside(i) == 2 % left side

        NonWatchside_UPOf3 = [Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
            Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), ...
            Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number))];
        
        Watchside_UPOf3 = [Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)),...
            Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)), ...
            Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number))];

        if Most_affected_side(i) == 1 || Most_affected_side(i) == 3 % re>li affected
            Most_affected_subj(i)=0; % watch is not worn on patient reported most-affected side
        elseif Most_affected_side(i) == 2 || Most_affected_side(i) == 4 % li>re affected
            Most_affected_subj(i)=1; % watch is worn on patient reported least-affected side
        elseif Most_affected_side(i) == 5 || Most_affected_side(i) == 6 % both side equally affected
            Most_affected_subj(i)=2; % watch is worn on one of equally affected sides
        end
    end

    Watchside_UPOf3_total = sum(Watchside_UPOf3); % total UPDRS 3 OFF score on watchside
    NonWatchside_UPOf3_total = sum(NonWatchside_UPOf3); % total UPDRS 3 OFF score on non-watchside

    if Watchside_UPOf3_total>NonWatchside_UPOf3_total
        Most_affected_obj(i)=1; % Watch on most_affected side
    elseif Watchside_UPOf3_total<NonWatchside_UPOf3_total
        Most_affected_obj(i)=0; % Watch on least affected side
    elseif Watchside_UPOf3_total==NonWatchside_UPOf3_total
        Most_affected_obj(i)=2; % Watch on one of equally affected sides
    end
    
    if Tremor_score_dev_side(i) < Tremor_score_nondev_side(i)
        Final_most_affected(i) = 0;
    elseif Tremor_score_dev_side(i) > Tremor_score_nondev_side(i)
        Final_most_affected(i) = 1;
    elseif Most_affected_subj(i) == 0 || Most_affected_subj(i) == 1
        Final_most_affected(i) = Most_affected_subj(i);
    elseif Most_affected_obj(i) == 0 || Most_affected_obj(i) == 1
        Final_most_affected(i) = Most_affected_obj(i);
    else
        Final_most_affected(i) = NaN;
    end

    UPDRS1a = [Visit1.Up1aCognit(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aHalPsy(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up1aDepres(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aAnxious(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up1aApathy(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aDopDysSyn(contains(Visit1.id,ppp_pep_number))];
    UPDRS1a_visit1 = sum(UPDRS1a,2);
    UPDRS1b = [Visit1.Updrs2It07(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It08(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It09(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It10(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It11(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It12(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It13(contains(Visit1.id,ppp_pep_number))];
    UPDRS1b_visit1 = sum(UPDRS1b,2);
    UPDRS1(i) = UPDRS1a_visit1 + UPDRS1b_visit1;

    UPDRS2_visit1 = [Visit1.Updrs2It14(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It15(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It16(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It17(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It18(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It19(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It20(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It21(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It22(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It23(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It24(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It25(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It26(contains(Visit1.id,ppp_pep_number))];
    UPDRS2(i) = sum(UPDRS2_visit1,2);

    UPDRS3OFF_visit1 = [Visit1.Up3OfSpeech(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFacial(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigNec(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfArise(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfGait(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFreez(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfStaPos(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPostur(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfSpont(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number)),Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpJaw(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfConstan(contains(Visit1.id,ppp_pep_number))];
    UPDRS3_OFF(i) = sum(UPDRS3OFF_visit1,2);

    UPDRS4(i) = Visit1.MotComTotal(contains(Visit1.id,ppp_pep_number));
    
    UPDRS_devside_tremor_OFF(i) = Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number));

    UPDRS3OFF_tremor = [Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number)),Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpJaw(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfConstan(contains(Visit1.id,ppp_pep_number))];
    UPDRS_tremor_OFF(i) = sum(UPDRS3OFF_tremor,2);

end

Disease_duration(Disease_duration<0) = NaN;

BaselineMedicated = [];
BaselineMedicated.ID = IDs_BaselineMedicated;
BaselineMedicated.Age = Age';
BaselineMedicated.Gender = Gender';
BaselineMedicated.Disease_duration = Disease_duration';
BaselineMedicated.Most_affected = Final_most_affected';
BaselineMedicated.UPDRS1 = UPDRS1';
BaselineMedicated.UPDRS2 = UPDRS2';
BaselineMedicated.UPDRS3_OFF = UPDRS3_OFF';
BaselineMedicated.UPDRS4 = UPDRS4';
BaselineMedicated.UPDRS_tremor_OFF = UPDRS_tremor_OFF';
BaselineMedicated.UPDRS_devside_tremor_OFF = UPDRS_devside_tremor_OFF';

BaselineMedicated = struct2table(BaselineMedicated);

%% Extract demographics and clinical characteristics for the unmedicated group

Age = [];
Gender = [];
Disease_duration = [];
Final_most_affected = [];
UPDRS1 = [];
UPDRS2 = [];
UPDRS3_OFF = [];
UPDRS4 = [];
UPDRS_tremor_OFF = [];
UPDRS_devside_tremor_OFF = [];
for i = 1:length(IDs_BaselineUnmedicated)

    ppp_pep_number = IDs_BaselineUnmedicated{i};

    if ~isempty(find(contains(Visit1PPP.id,ppp_pep_number)))

        Visit1 = Visit1PPP;
        
        Age(i) = Visit1.Age(contains(Visit1.id,ppp_pep_number));
        Gender(i) = Visit1.Gender(contains(Visit1.id,ppp_pep_number));
        Disease_duration(i) = Visit1.MonthSinceDiag(contains(Visit1.id,ppp_pep_number));

        Tremor_score_dev_side(i) = Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number));
        Tremor_score_nondev_side(i) = Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number));

        Watchside(i) = Visit1.WatchSide(contains(Visit1.id,ppp_pep_number)); % side on which the watch is worn
        Most_affected_side(i) = Visit1.MostAffSide(contains(Visit1.id,ppp_pep_number));

        if Watchside(i) == 1 % right side

            Watchside_UPOf3 = [Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
                Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), ...
                Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number))];

            NonWatchside_UPOf3 = [Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)),...
                Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)), ...
                Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number))];

            if Most_affected_side(i) == 1 || Most_affected_side(i) == 3  % re>li affected
                Most_affected_subj(i)=1; % watch is worn on patient reported most-affected side
            elseif Most_affected_side(i) == 2 || Most_affected_side(i) == 4 % li>re affected
                Most_affected_subj(i)=0;  % watch is not worn on patient reported least-affected side
            elseif Most_affected_side(i) == 5 || Most_affected_side(i) == 6 % both sides equally affected
                Most_affected_subj(i)=2; % watch is worn on one of equally affected sides
            end

        elseif Watchside(i) == 2 % left side

            NonWatchside_UPOf3 = [Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
                Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), ...
                Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number))];

            Watchside_UPOf3 = [Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)),...
                Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)), ...
                Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number))];

            if Most_affected_side(i) == 1 || Most_affected_side(i) == 3 % re>li affected
                Most_affected_subj(i)=0; % watch is not worn on patient reported most-affected side
            elseif Most_affected_side(i) == 2 || Most_affected_side(i) == 4 % li>re affected
                Most_affected_subj(i)=1; % watch is worn on patient reported least-affected side
            elseif Most_affected_side(i) == 5 || Most_affected_side(i) == 6 % both side equally affected
                Most_affected_subj(i)=2; % watch is worn on one of equally affected sides
            end
        end

        Watchside_UPOf3_total = sum(Watchside_UPOf3); % total UPDRS 3 OFF score on watchside
        NonWatchside_UPOf3_total = sum(NonWatchside_UPOf3); % total UPDRS 3 OFF score on non-watchside

        if Watchside_UPOf3_total>NonWatchside_UPOf3_total
            Most_affected_obj(i)=1; % Watch on most_affected side
        elseif Watchside_UPOf3_total<NonWatchside_UPOf3_total
            Most_affected_obj(i)=0; % Watch on least affected side
        elseif Watchside_UPOf3_total==NonWatchside_UPOf3_total
            Most_affected_obj(i)=2; % Watch on one of equally affected sides
        end

        if Tremor_score_dev_side(i) < Tremor_score_nondev_side(i)
            Final_most_affected(i) = 0;
        elseif Tremor_score_dev_side(i) > Tremor_score_nondev_side(i)
            Final_most_affected(i) = 1;
        elseif Most_affected_subj(i) == 0 || Most_affected_subj(i) == 1
            Final_most_affected(i) = Most_affected_subj(i);
        elseif Most_affected_obj(i) == 0 || Most_affected_obj(i) == 1
            Final_most_affected(i) = Most_affected_obj(i);
        else
            Final_most_affected(i) = NaN;
        end

            UPDRS1a = [Visit1.Up1aCognit(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aHalPsy(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up1aDepres(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aAnxious(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up1aApathy(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aDopDysSyn(contains(Visit1.id,ppp_pep_number))];
    UPDRS1a_visit1 = sum(UPDRS1a,2);
    UPDRS1b = [Visit1.Updrs2It07(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It08(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It09(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It10(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It11(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It12(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It13(contains(Visit1.id,ppp_pep_number))];
    UPDRS1b_visit1 = sum(UPDRS1b,2);
    UPDRS1(i) = UPDRS1a_visit1 + UPDRS1b_visit1;

    UPDRS2_visit1 = [Visit1.Updrs2It14(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It15(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It16(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It17(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It18(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It19(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It20(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It21(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It22(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It23(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It24(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It25(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It26(contains(Visit1.id,ppp_pep_number))];
    UPDRS2(i) = sum(UPDRS2_visit1,2);

    UPDRS3OFF_visit1 = [Visit1.Up3OfSpeech(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFacial(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigNec(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfArise(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfGait(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFreez(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfStaPos(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPostur(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfSpont(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number)),Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpJaw(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfConstan(contains(Visit1.id,ppp_pep_number))];
    UPDRS3_OFF(i) = sum(UPDRS3OFF_visit1,2);

    UPDRS4(i) = Visit1.MotComTotal(contains(Visit1.id,ppp_pep_number));
    
    UPDRS_devside_tremor_OFF(i) = Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number));

    UPDRS3OFF_tremor = [Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number)),Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpJaw(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfConstan(contains(Visit1.id,ppp_pep_number))];
    UPDRS_tremor_OFF(i) = sum(UPDRS3OFF_tremor,2);

    else

        Visit1 = Visit1DeNovo;

        Age(i) = Visit1.Age(contains(Visit1.id,ppp_pep_number));
        Gender(i) = Visit1.Gender(contains(Visit1.id,ppp_pep_number));

        Disease_duration_DeNovo = between(datetime(Visit1.DiagParkYear(contains(Visit1.id,ppp_pep_number)),Visit1.DiagParkMonth(contains(Visit1.id,ppp_pep_number)),1),...
            datetime(Visit1.AssessYear(contains(Visit1.id,ppp_pep_number)),Visit1.AssessMonth(contains(Visit1.id,ppp_pep_number)),1),'months');
        Disease_duration(i) = split(Disease_duration_DeNovo,{'months'});

        Tremor_score_dev_side(i) = Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number));
        Tremor_score_nondev_side(i) = Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number))+...
            Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number));

        Watchside(i) = Visit1.WatchSide(contains(Visit1.id,ppp_pep_number)); % side on which the watch is worn
        Most_affected_side(i) = Visit1.MostAffSide(contains(Visit1.id,ppp_pep_number));

        if Watchside(i) == 1 % right side

            Watchside_UPOf3 = [Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
                Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), ...
                Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number))];

            NonWatchside_UPOf3 = [Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)),...
                Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)), ...
                Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number))];

            if Most_affected_side(i) == 1 || Most_affected_side(i) == 3  % re>li affected
                Most_affected_subj(i)=1; % watch is worn on patient reported most-affected side
            elseif Most_affected_side(i) == 2 || Most_affected_side(i) == 4 % li>re affected
                Most_affected_subj(i)=0;  % watch is not worn on patient reported least-affected side
            elseif Most_affected_side(i) == 5 || Most_affected_side(i) == 6 % both sides equally affected
                Most_affected_subj(i)=2; % watch is worn on one of equally affected sides
            end

        elseif Watchside(i) == 2 % left side

            NonWatchside_UPOf3 = [Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
                Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), ...
                Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number))];

            Watchside_UPOf3 = [Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)),...
                Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)), ...
                Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number))];

            if Most_affected_side(i) == 1 || Most_affected_side(i) == 3 % re>li affected
                Most_affected_subj(i)=0; % watch is not worn on patient reported most-affected side
            elseif Most_affected_side(i) == 2 || Most_affected_side(i) == 4 % li>re affected
                Most_affected_subj(i)=1; % watch is worn on patient reported least-affected side
            elseif Most_affected_side(i) == 5 || Most_affected_side(i) == 6 % both side equally affected
                Most_affected_subj(i)=2; % watch is worn on one of equally affected sides
            end
        end

        Watchside_UPOf3_total = sum(Watchside_UPOf3); % total UPDRS 3 OFF score on watchside
        NonWatchside_UPOf3_total = sum(NonWatchside_UPOf3); % total UPDRS 3 OFF score on non-watchside

        if Watchside_UPOf3_total>NonWatchside_UPOf3_total
            Most_affected_obj(i)=1; % Watch on most_affected side
        elseif Watchside_UPOf3_total<NonWatchside_UPOf3_total
            Most_affected_obj(i)=0; % Watch on least affected side
        elseif Watchside_UPOf3_total==NonWatchside_UPOf3_total
            Most_affected_obj(i)=2; % Watch on one of equally affected sides
        end

        if Tremor_score_dev_side(i) < Tremor_score_nondev_side(i)
            Final_most_affected(i) = 0;
        elseif Tremor_score_dev_side(i) > Tremor_score_nondev_side(i)
            Final_most_affected(i) = 1;
        elseif Most_affected_subj(i) == 0 || Most_affected_subj(i) == 1
            Final_most_affected(i) = Most_affected_subj(i);
        elseif Most_affected_obj(i) == 0 || Most_affected_obj(i) == 1
            Final_most_affected(i) = Most_affected_obj(i);
        else
            Final_most_affected(i) = NaN;
        end

            UPDRS1a = [Visit1.Up1aCognit(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aHalPsy(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up1aDepres(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aAnxious(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up1aApathy(contains(Visit1.id,ppp_pep_number)), Visit1.Up1aDopDysSyn(contains(Visit1.id,ppp_pep_number))];
    UPDRS1a_visit1 = sum(UPDRS1a,2);
    UPDRS1b = [Visit1.Updrs2It07(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It08(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It09(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It10(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It11(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It12(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It13(contains(Visit1.id,ppp_pep_number))];
    UPDRS1b_visit1 = sum(UPDRS1b,2);
    UPDRS1(i) = UPDRS1a_visit1 + UPDRS1b_visit1;

    UPDRS2_visit1 = [Visit1.Updrs2It14(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It15(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It16(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It17(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It18(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It19(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It20(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It21(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It22(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It23(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Updrs2It24(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It25(contains(Visit1.id,ppp_pep_number)),Visit1.Updrs2It26(contains(Visit1.id,ppp_pep_number))];
    UPDRS2(i) = sum(UPDRS2_visit1,2);

    UPDRS3OFF_visit1 = [Visit1.Up3OfSpeech(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFacial(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigNec(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRue(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigLue(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRigRle(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRigLle(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFiTaNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfFiTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfHaMoNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfHaMoYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfProSNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfProSYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfToTaNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfToTaYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfLAgiNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfLAgiYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfArise(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfGait(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfFreez(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfStaPos(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPostur(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfSpont(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number)),Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpJaw(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfConstan(contains(Visit1.id,ppp_pep_number))];
    UPDRS3_OFF(i) = sum(UPDRS3OFF_visit1,2);

    UPDRS4(i) = Visit1.MotComTotal(contains(Visit1.id,ppp_pep_number));
    
    UPDRS_devside_tremor_OFF(i) = Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number));

    UPDRS3OFF_tremor = [Visit1.Up3OfPosTNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfPosTYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfKinTreNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfKinTreYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number)),Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number)),...
        Visit1.Up3OfRAmpJaw(contains(Visit1.id,ppp_pep_number)), Visit1.Up3OfConstan(contains(Visit1.id,ppp_pep_number))];
    UPDRS_tremor_OFF(i) = sum(UPDRS3OFF_tremor,2);
    end
end

Disease_duration(Disease_duration<0) = NaN;

BaselineUnmedicated = [];
BaselineUnmedicated.ID = IDs_BaselineUnmedicated;
BaselineUnmedicated.Age = Age';
BaselineUnmedicated.Gender = Gender';
BaselineUnmedicated.Disease_duration = Disease_duration';
BaselineUnmedicated.Most_affected = Final_most_affected';
BaselineUnmedicated.UPDRS1 = UPDRS1';
BaselineUnmedicated.UPDRS2 = UPDRS2';
BaselineUnmedicated.UPDRS3_OFF = UPDRS3_OFF';
BaselineUnmedicated.UPDRS4 = UPDRS4';
BaselineUnmedicated.UPDRS_tremor_OFF = UPDRS_tremor_OFF';
BaselineUnmedicated.UPDRS_devside_tremor_OFF = UPDRS_devside_tremor_OFF';

BaselineUnmedicated = struct2table(BaselineUnmedicated);

%% Save descriptives
save('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Descriptives.mat',"BaselineUnmedicated","BaselineMedicated")

%% Age
median(BaselineMedicated.Age)
prctile(BaselineMedicated.Age,25)
prctile(BaselineMedicated.Age,75)

median(BaselineUnmedicated.Age)
prctile(BaselineUnmedicated.Age,25)
prctile(BaselineUnmedicated.Age,75)

%% Gender
length(find(BaselineMedicated.Gender==1))
length(find(BaselineUnmedicated.Gender==1))

%% Disease duration
median(BaselineMedicated.Disease_duration,'omitnan')
prctile(BaselineMedicated.Disease_duration,25)
prctile(BaselineMedicated.Disease_duration,75)

median(BaselineUnmedicated.Disease_duration,'omitnan')
prctile(BaselineUnmedicated.Disease_duration,25)
prctile(BaselineUnmedicated.Disease_duration,75)

%% Most-affected
n1 = length(find(BaselineMedicated.Most_affected==1))
n2 = length(find(BaselineMedicated.Most_affected==0))
n3 = length(find(BaselineUnmedicated.Most_affected==1))
n4 = length(find(BaselineUnmedicated.Most_affected==0))

% most_affected = [1 * ones(n1,1); 2 * ones(n2,1); 1 * ones(n3,1); 2 * ones(n4,1)]; % 1 = Male, 2 = Female
% group  = [1 * ones(77,1); 2 * ones(76,1)]; % 1 = Group 1, 2 = Group 2
% [table, chi2, p, labels] = crosstab(most_affected, group)

%% UPDRS 1

median(BaselineMedicated.UPDRS1,'omitnan')
prctile(BaselineMedicated.UPDRS1,25)
prctile(BaselineMedicated.UPDRS1,75)

median(BaselineUnmedicated.UPDRS1,'omitnan')
prctile(BaselineUnmedicated.UPDRS1,25)
prctile(BaselineUnmedicated.UPDRS1,75)

%% UPDRS 2
median(BaselineMedicated.UPDRS2,'omitnan')
prctile(BaselineMedicated.UPDRS2,25)
prctile(BaselineMedicated.UPDRS2,75)

median(BaselineUnmedicated.UPDRS2,'omitnan')
prctile(BaselineUnmedicated.UPDRS2,25)
prctile(BaselineUnmedicated.UPDRS2,75)

%% UPDRS 3
median(BaselineMedicated.UPDRS3_OFF,'omitnan')
prctile(BaselineMedicated.UPDRS3_OFF,25)
prctile(BaselineMedicated.UPDRS3_OFF,75)

median(BaselineUnmedicated.UPDRS3_OFF,'omitnan')
prctile(BaselineUnmedicated.UPDRS3_OFF,25)
prctile(BaselineUnmedicated.UPDRS3_OFF,75)

%% UPDRS 4
median(BaselineMedicated.UPDRS4,'omitnan')
prctile(BaselineMedicated.UPDRS4,25)
prctile(BaselineMedicated.UPDRS4,75)

median(BaselineUnmedicated.UPDRS4,'omitnan')
prctile(BaselineUnmedicated.UPDRS4,25)
prctile(BaselineUnmedicated.UPDRS4,75)

%% UPDRS 3 tremor subscore
median(BaselineMedicated.UPDRS_tremor_OFF,'omitnan')
prctile(BaselineMedicated.UPDRS_tremor_OFF,25)
prctile(BaselineMedicated.UPDRS_tremor_OFF,75)

median(BaselineUnmedicated.UPDRS_tremor_OFF,'omitnan')
prctile(BaselineUnmedicated.UPDRS_tremor_OFF,25)
prctile(BaselineUnmedicated.UPDRS_tremor_OFF,75)

%% Rest tremor device sided arm
median(BaselineMedicated.UPDRS_devside_tremor_OFF,'omitnan')
prctile(BaselineMedicated.UPDRS_devside_tremor_OFF,25)
prctile(BaselineMedicated.UPDRS_devside_tremor_OFF,75)

median(BaselineUnmedicated.UPDRS_devside_tremor_OFF,'omitnan')
prctile(BaselineUnmedicated.UPDRS_devside_tremor_OFF,25)
prctile(BaselineUnmedicated.UPDRS_devside_tremor_OFF,75)

%% Create survival curve (until treatment initiation)
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Inclusion.mat');
start_week = Inclusion.StartWeek(ismember(Inclusion.ID,IDs_BaselineUnmedicated));
n = length(start_week);

% Sort weeks
sortedWeeks = sort(start_week);

% Compute number remaining at each unique week
[counts, uniqueWeeks] = groupcounts(sortedWeeks);

numRemaining = n - cumsum(counts);

% Plot
figure;
stairs([0; uniqueWeeks-1], [n; numRemaining],'k', 'LineWidth', 1.5);
xlabel('Weeks since baseline');
ylabel('Number of participants in unmedicated group');
ylim([0 80])
xlim([0 100])