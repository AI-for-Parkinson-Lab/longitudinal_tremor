

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