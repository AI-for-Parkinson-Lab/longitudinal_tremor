%% Step 8: Perform regression analysis to assess the effect of dopaminergic medication
% Covariates include change in LEDD, disease duration and watch side (more- or less- affected side)
clear all; close all;

%% Load sensor data and IDs
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Trends_filled.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IDs_selected.mat'); 
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\visit_week_numbers.mat') % load visit week numbers

%% Load clinical data (only data from PPP as only the baseline medicated participants are used)
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

%% Extract disease duration and watchside

Visit1 = Visit1PPP;
Final_most_affected = [];
MonthSinceDiag = [];

for i = 1:length(IDs_BaselineMedicated)

    ppp_pep_number = IDs_BaselineMedicated(i);
    
    MonthSinceDiag(i) = Visit1.MonthSinceDiag(contains(Visit1.id,ppp_pep_number));
   
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

% Apply corrections (negative disease duration)
MonthSinceDiag(MonthSinceDiag<0) = NaN;

%% Load LEDD values
LEDD_visit1 = [];
LEDD_visit2 = [];
LEDD_visit3 = [];
LEDDPPP = readtable("C:\Users\z835211\Documents\Data\PPP\Medication\LEDD_PPP.csv");

for i = 1:length(IDs_BaselineMedicated)

    id = IDs_BaselineMedicated{i};

    if ~isempty(find(ismember(LEDDPPP.id,id)))
        
        if LEDDPPP.Visit1_missing_med(ismember(LEDDPPP.id,id))>0
            LEDD_visit1(i) = 0;
        else
            LEDD_visit1(i) = LEDDPPP.Visit1(ismember(LEDDPPP.id,id));
        end
        visit2_week_number = cell2mat(visit_week_numbers.visit2(contains(visit_week_numbers.ID,id)));
        if visit2_week_number <= 60
            if LEDDPPP.Visit2_missing_med(ismember(LEDDPPP.id,id))>0
                LEDD_visit2(i) = 0;
            else
                LEDD_visit2(i) = LEDDPPP.Visit2(ismember(LEDDPPP.id,id));
            end
        else
            LEDD_visit2(i) = NaN;
        end

        visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));

        if visit3_week_number <= 120
            if LEDDPPP.Visit3_missing_med(ismember(LEDDPPP.id,id))>0
                LEDD_visit3(i) = 0;
            else
                LEDD_visit3(i) = LEDDPPP.Visit3(ismember(LEDDPPP.id,id));
            end
        else
            LEDD_visit3(i) = NaN;
        end
    else
        LEDD_visit1(i) = NaN;
        LEDD_visit2(i) = NaN;
        LEDD_visit3(i) = NaN;
    end
end

IDs_incorrect_baseline_LEDD = {'POMU066326B8F70E150E','POMU0A109E0D97672361','POMU428FEF5AA8B909DC','POMU6080FFB910C10DC3','POMUA249715D7C6FDE8F'};
LEDD_visit1(ismember(IDs_BaselineMedicated,IDs_incorrect_baseline_LEDD)) = NaN;

%% Compute 2-year changes

start_idx = 2;
end_idx = 51;

idx = contains([IDs_BaselineMedicated; IDs_BaselineUnmedicated],IDs_BaselineMedicated);

Delta = [];
Delta.tremor_time = logit(trend_tremor_time_medicated_filled(idx,end_idx)) - logit(trend_tremor_time_medicated_filled(idx,start_idx));
Delta.tremor_power_mode = trend_modal_tremor_power_medicated_filled(idx,end_idx) - trend_modal_tremor_power_medicated_filled(idx,start_idx);
Delta.tremor_power_90th = trend_perc90_tremor_power_medicated_filled(idx,end_idx) - trend_perc90_tremor_power_medicated_filled(idx,start_idx);
Delta.LEDD = LEDD_visit3(idx)' - LEDD_visit1(idx)';
Delta.Disease_duration = MonthSinceDiag(idx)';
Delta.Least_affected = 1 - Final_most_affected(idx)';
Delta_table = struct2table(Delta);


%% Multivariable linear regression
close all;

Delta_standardized = (Delta_table.Variables - mean(Delta_table.Variables,'omitnan'))./std(Delta_table.Variables,'omitnan'); % Standardize variables
X = [ones(size(Delta_standardized,1),1) Delta_standardized(:,[4 5 6])];

delta_tremor_time = Delta_standardized(:,1);
delta_modal_tremor_power = Delta_standardized(:,2);
delta_tremor_power_90th = Delta_standardized(:,3);

[b,bint,r,rint,stats] = regress(delta_tremor_time,X);
b_standardized_tremor_time = b(2:end)'
bint_standardized_tremor_time = bint(2:end,:)'

% Check residuals
% figure(); histogram(r)
% figure(); scatter(delta_tremor_time,r)
% [h,p] = kstest(r)

[b,bint,r,rint,stats] = regress(delta_modal_tremor_power,X);
b_standardized_modal_tremor_power = b(2:end)'
bint_standardized_modal_tremor_power = bint(2:end,:)'

[b,bint,r,rint,stats] = regress(delta_tremor_power_90th,X);
b_standardized_tremor_power_90th = b(2:end)'
bint_standardized_tremor_power_90th = bint(2:end,:)'

% Check residuals
% figure(); histogram(r)
% figure(); scatter(delta_tremor_power_90th,r)
% [h,p] = kstest(r)


%% Plot coefficients 

cmap = colororder();
n=3;
m = 11;
figure(); hold on;
for i = 1:n
    line([bint_standardized_tremor_time(1,i), bint_standardized_tremor_time(2,i)], [m, m], 'Color',cmap(i,:), 'LineWidth', 1.5); % CI line
    plot(b_standardized_tremor_time(i), m, 'ko', 'MarkerFaceColor', cmap(i,:), 'MarkerSize', 6); % SRM point

    line([bint_standardized_modal_tremor_power(1,i), bint_standardized_modal_tremor_power(2,i)], [m-4, m-4], 'Color', cmap(i,:), 'LineWidth', 1.5); % CI line
    plot(b_standardized_modal_tremor_power(i), m-4, 'ko', 'MarkerFaceColor', cmap(i,:), 'MarkerSize', 6); % SRM point

    line([bint_standardized_tremor_power_90th(1,i), bint_standardized_tremor_power_90th(2,i)], [m-8, m-8], 'Color', cmap(i,:), 'LineWidth', 1.5); % CI line
    plot(b_standardized_tremor_power_90th(i), m-8, 'ko', 'MarkerFaceColor', cmap(i,:), 'MarkerSize', 6); % SRM point

    m = m-1;
end

% Add labels
yticks([2 6 10])
yticklabels({'\Delta90th percentile of tremor power (n=122)','\Deltamodal tremor power (n=122)','\Deltatremor time (n=359)'})
xline(0,'--')
ylim([0 12])
xlim([-0.5 0.5])
xlabel('Standardized \beta coefficients with 95%-CI')
legend({'\DeltaLEDD','','','','','','Disease duration','','','','','','Least-affected side'}')

%% Scatterplot

Delta = Delta_table.Variables;
X = [ones(size(Delta,1),1) Delta(:,[4 5 6])];

delta_tremor_time = Delta(:,1);

[b_tremor_time] = regress(delta_tremor_time,X);

corrected_delta_tremor_time = delta_tremor_time - b_tremor_time(2)*X(:,2) - b_tremor_time(4)*X(:,4);

figure(); hold on;
scatter(X(:,3),corrected_delta_tremor_time,'filled','k','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.4)
plot(0:90,b_tremor_time(1)+[0:90]*b_tremor_time(3),'k','LineWidth',2)
yline(0,'k--')
text(70,0.3,'Increase in tremor time')
text(70,-0.25,'Decrease in tremor time')
ylabel('\Delta tremor time (log odds-ratio)')
xlabel('Disease duration (months)')


