%% 
clear all; close all;

%% Load IDs
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\IDs.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250402_l1trend\IDs_BaselineUnmedicated_L1trend.mat')

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

%% Extract clinical info baseline medicated group

Disease_duration = [];
Age = [];
Final_most_affected = [];
Visit1 = Visit1PPP;

for i = 1:length(IDs_BaselineMedicated)

    ppp_pep_number = IDs_BaselineMedicated{i};

    Disease_duration(i) = Visit1.MonthSinceDiag(contains(Visit1.id,ppp_pep_number));
    Age(i) = Visit1.Age(contains(Visit1.id,ppp_pep_number));
   
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

end

Disease_duration(Disease_duration<0) = NaN;

BaselineMedicated = [];
BaselineMedicated.ID = IDs_BaselineMedicated;
BaselineMedicated.Disease_duration = Disease_duration';
BaselineMedicated.Age = Age';
BaselineMedicated.Most_affected = Final_most_affected';

BaselineMedicated = struct2table(BaselineMedicated);

%% Extract clinical info baseline unmedicated group

Disease_duration = [];
Age = [];
Final_most_affected = [];

for i = 1:length(IDs_BaselineUnmedicated_L1trend)

    ppp_pep_number = IDs_BaselineUnmedicated_L1trend{i};

    if ~isempty(find(contains(Visit1PPP.id,ppp_pep_number)))

        Visit1 = Visit1PPP;

        Disease_duration(i) = Visit1.MonthSinceDiag(contains(Visit1.id,ppp_pep_number));
        Age(i) = Visit1.Age(contains(Visit1.id,ppp_pep_number));

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

    else

        Visit1 = Visit1DeNovo;

        Age(i) = Visit1.Age(contains(Visit1.id,ppp_pep_number));

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
    end
end

Disease_duration(Disease_duration<0) = NaN;

BaselineUnmedicated = [];
BaselineUnmedicated.ID = IDs_BaselineUnmedicated_L1trend;
BaselineUnmedicated.Disease_duration = Disease_duration';
BaselineUnmedicated.Age = Age';
BaselineUnmedicated.Most_affected = Final_most_affected';

BaselineUnmedicated = struct2table(BaselineUnmedicated);

%% Propensity score matching

treatment = [zeros(1,length(IDs_BaselineMedicated)) ones(1,length(IDs_BaselineUnmedicated_L1trend))]';
X  = [BaselineMedicated(:,2:4).Variables; BaselineUnmedicated(:,2:4).Variables];

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