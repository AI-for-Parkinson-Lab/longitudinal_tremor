%% 
clear all; close all;
%%
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\IDs.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250523_L1trend\trends2.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\MATLAB\Tremor progression\Paper\20250327_inclusion_exclusion\visit_week_numbers.mat')

%% Load clinical data
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

%% Extract clinical information
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

    % Rest_tremor_OFF(i) = Visit1.Up3OfRAmpLegNonDev(contains(Visit1.id,ppp_pep_number)) + Visit1.Up3OfRAmpLegYesDev(contains(Visit1.id,ppp_pep_number)) + Visit1.Up3OfRAmpArmYesDev(contains(Visit1.id,ppp_pep_number)) + Visit1.Up3OfRAmpArmNonDev(contains(Visit1.id,ppp_pep_number)) + Visit1.Up3OfConstan(contains(Visit1.id,ppp_pep_number));
    % Rest_tremor_ON(i) = Visit1.Up3OnRAmpLegNonDev(contains(Visit1.id,ppp_pep_number)) + Visit1.Up3OnRAmpLegYesDev(contains(Visit1.id,ppp_pep_number)) + Visit1.Up3OnRAmpArmYesDev(contains(Visit1.id,ppp_pep_number)) + Visit1.Up3OnRAmpArmNonDev(contains(Visit1.id,ppp_pep_number)) + Visit1.Up3OnConstan(contains(Visit1.id,ppp_pep_number));
    Rest_tremor_OFF(i) = Visit1.Up3OfConstan(contains(Visit1.id,ppp_pep_number));
    Rest_tremor_ON(i) = Visit1.Up3OnConstan(contains(Visit1.id,ppp_pep_number));
    Proportional_change(i) = (Rest_tremor_OFF(i)  - Rest_tremor_ON(i))/Rest_tremor_OFF(i);

    if isnan(Proportional_change(i))
        Levadopa_responsiveness(i) = NaN;
    elseif Proportional_change(i) >= 0.25
        Levadopa_responsiveness(i) = 1;
    else
        Levadopa_responsiveness(i) = 0;
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

IDs_incorrect_baseline_LEDD = ['POMU066326B8F70E150E','POMU0A109E0D97672361','POMU428FEF5AA8B909DC','POMU6080FFB910C10DC3','POMUA249715D7C6FDE8F'];
LEDD_visit1(ismember(IDs_BaselineMedicated,IDs_incorrect_baseline_LEDD)) = NaN;

%% Create delta table (2 year)
start_idx = 2;
end_idx = 51;

% IDs_BaselineMedicated_tremor = IDs_BaselineMedicated(all(~isnan(trend_median_tremor_power_medicated_filled(1:length(IDs_BaselineMedicated),[2 51])),2));

idx = contains([IDs_BaselineMedicated; IDs_StartMedication],IDs_BaselineMedicated);

Delta = [];
Delta.tremor_time = logit(trend_tremor_time_medicated_filled(idx,end_idx)) - logit(trend_tremor_time_medicated_filled(idx,start_idx));
Delta.tremor_power_median = trend_median_tremor_power_medicated_filled(idx,end_idx) - trend_median_tremor_power_medicated_filled(idx,start_idx);
Delta.tremor_power_mode = trend_modal_tremor_power_medicated_filled(idx,end_idx) - trend_modal_tremor_power_medicated_filled(idx,start_idx);
Delta.tremor_power_90th = trend_perc90_tremor_power_medicated_filled(idx,end_idx) - trend_perc90_tremor_power_medicated_filled(idx,start_idx);
Delta.LEDD = LEDD_visit3(idx)' - LEDD_visit1(idx)';
Delta.Disease_duration = MonthSinceDiag(idx)';
Delta.Least_affected = 1 - Final_most_affected(idx)';
Delta.Levadopa_responsiveness = Levadopa_responsiveness(idx)';
Delta_table = struct2table(Delta);
% Delta_table.Interaction_LEDD_Disease_duration = Delta_table.LEDD.*Delta_table.Disease_duration;
Delta_table.LEDD_responsiveness_interaction = Delta_table.LEDD.*Delta_table.Levadopa_responsiveness;


%% Multivariable linear regression
close all;

Delta_standardized = (Delta_table.Variables - mean(Delta_table.Variables,'omitnan'))./std(Delta_table.Variables,'omitnan');
X = [ones(size(Delta_standardized,1),1) Delta_standardized(:,[5 6 7 8 9])];

delta_tremor_time = Delta_standardized(:,1);
delta_median_tremor_power = Delta_standardized(:,2);
delta_modal_tremor_power = Delta_standardized(:,3);
delta_tremor_power_90th = Delta_standardized(:,4);

% use same group for tremor time
% delta_tremor_time(isnan(delta_tremor_power_90th)) = NaN; 

[b,bint,r,rint,stats] = regress(delta_tremor_time,X);
b_standardized_tremor_time = b(2:end)'
bint_standardized_tremor_time = bint(2:end,:)'

mdl = fitlm(X, delta_tremor_time);
disp(mdl.Coefficients)

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

%%
figure(); hold on;
scatter(Delta_table.LEDD(Delta_table.Levadopa_responsiveness==1),Delta_table.tremor_time(Delta_table.Levadopa_responsiveness==1))
scatter(Delta_table.LEDD(Delta_table.Levadopa_responsiveness==0),Delta_table.tremor_time(Delta_table.Levadopa_responsiveness==0))




%% Forest plot
% Names = { '\DeltaLEDD','Disease duration','Least-affected side'};
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

%% Scatterplots

Delta = Delta_table.Variables;
X = [ones(size(Delta,1),1) Delta(:,[5 6 7])];

delta_tremor_time = Delta(:,1);

% delta_tremor_time(isnan(Delta(:,2)))=NaN;

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

%%
Delta = Delta_table.Variables;
X = [ones(size(Delta,1),1) Delta(:,[5 6 7])];

delta_tremor_power_90th = Delta(:,4);

[b_tremor_power] = regress(delta_tremor_power_90th,X);

corrected_delta_tremor_power = delta_tremor_power_90th - b_tremor_power(2)*X(:,2) - b_tremor_power(4)*X(:,4) ;

figure(); hold on;
scatter(X(:,3),corrected_delta_tremor_power,'filled','k','MarkerFaceAlpha',0.4,['MarkerEdgeAlpha'],0.4)
plot(0:90,b_tremor_power(1)+[0:90]*b_tremor_power(3),'k','LineWidth',2)
yline(0,'k--')
text(70,0.3,'Increase in tremor power')
text(70,-0.3,'Decrease in tremor power')
ylabel('\Delta 90th percentile of tremor power')
xlabel('Disease duration (months)')

