%% Step 7: Assess correlation to changes in patient-reported tremor score
clear all; close all;

%% Load sensor data and IDs
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\Trends_filled.mat')
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\IDs_selected.mat'); 

%% Merge medicated and unmedicated data
tremor_time = trend_tremor_time_medicated_filled;
tremor_time_unmedicated_total = trend_tremor_time_medicated_filled(length(IDs_BaselineMedicated)+1:end,:);
tremor_time_unmedicated_total(isnan(tremor_time_unmedicated_total)) = trend_tremor_time_unmedicated_filled(isnan(tremor_time_unmedicated_total));
tremor_time(length(IDs_BaselineMedicated)+1:end,:) = tremor_time_unmedicated_total;

modal_tremor_power = trend_modal_tremor_power_medicated_filled;
modal_tremor_power_unmedicated_total = trend_modal_tremor_power_medicated_filled(length(IDs_BaselineMedicated)+1:end,:);
modal_tremor_power_unmedicated_total(isnan(modal_tremor_power_unmedicated_total)) = trend_modal_tremor_power_unmedicated_filled(isnan(modal_tremor_power_unmedicated_total));
modal_tremor_power(length(IDs_BaselineMedicated)+1:end,:) = modal_tremor_power_unmedicated_total;

perc90_tremor_power = trend_perc90_tremor_power_medicated_filled;
perc90_tremor_power_unmedicated_total = trend_perc90_tremor_power_medicated_filled(length(IDs_BaselineMedicated)+1:end,:);
perc90_tremor_power_unmedicated_total(isnan(perc90_tremor_power_unmedicated_total)) = trend_perc90_tremor_power_unmedicated_filled(isnan(perc90_tremor_power_unmedicated_total));
perc90_tremor_power(length(IDs_BaselineMedicated)+1:end,:) = perc90_tremor_power_unmedicated_total;

%% Calculate 2-year delta of whole group

delta_tremor_time = logit(tremor_time(:,51)) - logit(tremor_time(:,2));
delta_modal_tremor_power = modal_tremor_power(:,51) - modal_tremor_power(:,2);
delta_perc90_tremor_power = perc90_tremor_power(:,51) - perc90_tremor_power(:,2);

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

%% Extract UPDRS scores
load('C:\Users\z835211\OneDrive - Radboudumc\Documents\Tremor progression paper\Matlab_results\visit_week_numbers.mat') % load visit week numbers

UPDRS_317OFF_1 = [];
UPDRS_318OFF_1 = [];
UPDRS_210_1 = [];

UPDRS_317OFF_3 = [];
UPDRS_318OFF_3 = [];
UPDRS_210_3 = [];

UPDRS_317ON_1 = [];
UPDRS_318ON_1 = [];
UPDRS_317ON_3 = [];
UPDRS_318ON_3 = [];

IDs_included = [IDs_BaselineMedicated; IDs_BaselineUnmedicated];

for i = 1:length(IDs_included)
    id = IDs_included{i};

    if ~isempty(find(contains(Visit1PPP.id,id)))

        UPDRS_317OFF_1(i) = Visit1PPP.Up3OfRAmpArmYesDev(contains(Visit1PPP.id,id));
        UPDRS_318OFF_1(i) = Visit1PPP.Up3OfConstan(contains(Visit1PPP.id,id));
        UPDRS_210_1(i) = Visit1PPP.Updrs2It23(contains(Visit1PPP.id,id));

        UPDRS_317ON_1(i) = Visit1PPP.Up3OnRAmpArmYesDev(contains(Visit1PPP.id,id));
        UPDRS_318ON_1(i) = Visit1PPP.Up3OnConstan(contains(Visit1PPP.id,id));

        if ~isempty(find(contains(Visit3PPP.id,id)))
            visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));
            if visit3_week_number <= 120
                UPDRS_317OFF_3(i) = Visit3PPP.Up3OfRAmpArmYesDev(contains(Visit3PPP.id,id));
                UPDRS_318OFF_3(i) = Visit3PPP.Up3OfConstan(contains(Visit3PPP.id,id));
                UPDRS_210_3(i) = Visit3PPP.Updrs2It23(contains(Visit3PPP.id,id));

                UPDRS_317ON_3(i) = Visit3PPP.Up3OnRAmpArmYesDev(contains(Visit3PPP.id,id));
                UPDRS_318ON_3(i) = Visit3PPP.Up3OnConstan(contains(Visit3PPP.id,id));
            else
                UPDRS_317OFF_3(i) = NaN;
                UPDRS_318OFF_3(i) = NaN;
                UPDRS_210_3(i) = NaN;
            end
        else
            UPDRS_317OFF_3(i) = NaN;
            UPDRS_318OFF_3(i) = NaN;
            UPDRS_210_3(i) = NaN;
        end

    else
        UPDRS_317OFF_1(i) = Visit1DeNovo.Up3OfRAmpArmYesDev(contains(Visit1DeNovo.id,id));
        UPDRS_318OFF_1(i) = Visit1DeNovo.Up3OfConstan(contains(Visit1DeNovo.id,id));
        UPDRS_210_1(i) = Visit1DeNovo.Updrs2It23(contains(Visit1DeNovo.id,id));
        UPDRS_317ON_1(i) =  NaN;
        UPDRS_318ON_1(i) =  NaN;
        UPDRS_317ON_3(i) =  NaN;
        UPDRS_318ON_3(i) =  NaN;

        if contains(id,'POMU600C11F136E6FB4D')

            UPDRS_317OFF_3(i) = Visit2DeNovo.Up3OfRAmpArmYesDev(contains(Visit2DeNovo.id,id));
            UPDRS_318OFF_3(i) = Visit2DeNovo.Up3OfConstan(contains(Visit2DeNovo.id,id));
            UPDRS_210_3(i) = Visit2DeNovo.Updrs2It23(contains(Visit2DeNovo.id,id));


        else
            if ~isempty(find(contains(Visit3DeNovo.id,id)))
                visit3_week_number = cell2mat(visit_week_numbers.visit3(contains(visit_week_numbers.ID,id)));
                if visit3_week_number <= 120
                    UPDRS_317OFF_3(i) = Visit3DeNovo.Up3OfRAmpArmYesDev(contains(Visit3DeNovo.id,id));
                    UPDRS_318OFF_3(i) = Visit3DeNovo.Up3OfConstan(contains(Visit3DeNovo.id,id));
                    UPDRS_210_3(i) = Visit3DeNovo.Updrs2It23(contains(Visit3DeNovo.id,id));
                else
                    UPDRS_317OFF_3(i) = NaN;
                    UPDRS_318OFF_3(i) = NaN;
                    UPDRS_210_3(i) = NaN;

                end
            else
                UPDRS_317OFF_3(i) = NaN;
                UPDRS_318OFF_3(i) = NaN;
                UPDRS_210_3(i) = NaN;
            end
        end
    end
end

%% Calculate 2-year delta of UPDRS scores

Delta_UPDRS317OFF = UPDRS_317OFF_3 - UPDRS_317OFF_1;
Delta_UPDRS318OFF = UPDRS_318OFF_3 - UPDRS_318OFF_1;
Delta_UPDRS317ON = UPDRS_317ON_3 - UPDRS_317ON_1;
Delta_UPDRS318ON = UPDRS_318ON_3 - UPDRS_318ON_1;
Delta_UPDRS210 = UPDRS_210_3 - UPDRS_210_1;

%% Spearman rank correlation
[Rho,pval] = corr([delta_tremor_time delta_modal_tremor_power delta_perc90_tremor_power Delta_UPDRS317OFF' Delta_UPDRS318OFF' Delta_UPDRS317ON' Delta_UPDRS318ON' Delta_UPDRS210' ],'Type','Spearman','rows','pairwise');
rho_values = round(Rho,2);

% Mask upper triangle
mask = triu(true(size(rho_values)), 1);
rho_values(mask) = NaN;
pval(mask) = NaN;

figure();
colormap('sky')
h = imagesc(rho_values);               % Create image object
h.AlphaData = ~isnan(rho_values);      % Only show non-NaN values
colorbar;
set(gca, 'XTick',1:8,'XTickLabel',{'\Delta tremor time','\Delta modal tremor power','\Delta 90th perc of tremor power',...
    '\Delta rest tremor severity in OFF','\Delta rest tremor constancy in OFF','\Delta rest tremor severity in ON','\Delta rest tremor constancy in ON',...
    '\Delta patient-reported tremor'});
set(gca, 'YTick',1:8,'YTickLabel',{'\Delta tremor time','\Delta modal tremor power','\Delta 90th perc of tremor power',...
    '\Delta rest tremor severity in OFF','\Delta rest tremor constancy in OFF','\Delta rest tremor severity in ON','\Delta rest tremor constancy in ON',...
    '\Delta patient-reported tremor'});
ax = gca;
ax.FontSize = 16;
ax.TickLength = [0,0];

% Add significance level annotations
ax = gca;
for i = 1:size(rho_values, 1)
    for j = 1:i-1 % Only loop through lower triangle
        if isnan(rho_values(i,j))
            continue;
        end
        % Determine significance level
        if pval(i, j) < 0.001
            sig_level = '***'; % Highly significant
        elseif pval(i, j) < 0.01
            sig_level = '**';  % Significant
        elseif pval(i, j) < 0.05
            sig_level = '*';   % Marginally significant
        else
            sig_level = '';    % Not significant
        end
        
        % Add text annotation
        text_str = sprintf('%.2f\n%s', rho_values(i, j), sig_level);
        text(j, i, text_str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color', 'black');
    end
end
clim([0 1])
box off

%% Find the number of participants used for each correlation 

length(find(~isnan(Delta_UPDRS317OFF) & ~isnan(Delta_UPDRS210)))