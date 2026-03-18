%% Step 2: Extract clinical and demographic characteristics for excluded participants
clear all; close all;

%% Load IDs
load('\\umcn.nl\nas\RBS\NEURO_AI4P\Users\Nienke Timmermans\Tremor progression\Derived_data\IDs_exluded_dropout.mat'); % load selected IDs
IDs_excluded = [IDs_ExcludedPD; IDs_MedicatedDropout; IDs_UnmedicatedDropout];

%% Load clinical data
% De Novo data:
Visit1DeNovo = readtable("\\umcn.nl\nas\RBS\NEURO_AI4P\Datasets\PPP_tremor\Tremor progression paper\clinical_data\Visit1_DeNovo.csv");
Visit2DeNovo = readtable("\\umcn.nl\nas\RBS\NEURO_AI4P\Datasets\PPP_tremor\Tremor progression paper\clinical_data\Visit2_DeNovo.csv");
Visit3DeNovo = readtable("\\umcn.nl\nas\RBS\NEURO_AI4P\Datasets\PPP_tremor\Tremor progression paper\clinical_data\Visit3_DeNovo.csv");

% PPP data:
Visit1PPP = readtable("\\umcn.nl\nas\RBS\NEURO_AI4P\Datasets\PPP_tremor\Tremor progression paper\clinical_data\General_visit1.csv");
Visit2PPP = readtable("\\umcn.nl\nas\RBS\NEURO_AI4P\Datasets\PPP_tremor\Tremor progression paper\clinical_data\General_visit2.csv");
Visit3PPP = readtable("\\umcn.nl\nas\RBS\NEURO_AI4P\Datasets\PPP_tremor\Tremor progression paper\clinical_data\General_visit3.csv");

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


%% Extract demographics and clinical characteristics

Age = [];
Gender = [];
Disease_duration = [];
UPDRS1 = [];
UPDRS2 = [];
UPDRS3_OFF = [];
UPDRS4 = [];
UPDRS_tremor_OFF = [];
UPDRS_devside_tremor_OFF = [];

for i = 1:length(IDs_excluded)

    ppp_pep_number = IDs_excluded{i};

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

Excluded = [];
Excluded.ID = IDs_excluded;
Excluded.Age = Age';
Excluded.Gender = Gender';
Excluded.Disease_duration = Disease_duration';
Excluded.UPDRS1 = UPDRS1';
Excluded.UPDRS2 = UPDRS2';
Excluded.UPDRS3_OFF = UPDRS3_OFF';
Excluded.UPDRS4 = UPDRS4';
Excluded.UPDRS_tremor_OFF = UPDRS_tremor_OFF';
Excluded.UPDRS_devside_tremor_OFF = UPDRS_devside_tremor_OFF';

Excluded = struct2table(Excluded);

%% Age
median(Excluded.Age,'omitnan')
prctile(Excluded.Age,25)
prctile(Excluded.Age,75)

%% Gender
length(find(Excluded.Gender==1))
length(find(Excluded.Gender==2))

%% Disease duration
median(Excluded.Disease_duration,'omitnan')
prctile(Excluded.Disease_duration,25)
prctile(Excluded.Disease_duration,75)


%% UPDRS 1
median(Excluded.UPDRS1,'omitnan')
prctile(Excluded.UPDRS1,25)
prctile(Excluded.UPDRS1,75)

%% UPDRS 2
median(Excluded.UPDRS2,'omitnan')
prctile(Excluded.UPDRS2,25)
prctile(Excluded.UPDRS2,75)

%% UPDRS 3
median(Excluded.UPDRS3_OFF,'omitnan')
prctile(Excluded.UPDRS3_OFF,25)
prctile(Excluded.UPDRS3_OFF,75)

%% UPDRS 4
median(Excluded.UPDRS4,'omitnan')
prctile(Excluded.UPDRS4,25)
prctile(Excluded.UPDRS4,75)

%% UPDRS 3 tremor subscore
median(Excluded.UPDRS_tremor_OFF,'omitnan')
prctile(Excluded.UPDRS_tremor_OFF,25)
prctile(Excluded.UPDRS_tremor_OFF,75)

%% Rest tremor device sided arm
median(Excluded.UPDRS_devside_tremor_OFF,'omitnan')
prctile(Excluded.UPDRS_devside_tremor_OFF,25)
prctile(Excluded.UPDRS_devside_tremor_OFF,75)
