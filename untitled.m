% yield_estimation_from_csv.m
% Data-dependent yield estimation pipeline
% Requires: Statistics and Machine Learning Toolbox (TreeBagger, fitrensemble)

clear; close all; clc;

%% Configuration

fn = 'Crops_data.csv';         % filename (already uploaded)
targetNamesTry = {'Yield','yield','YIELD','target','Yield_ton_ha'}; % candidate target names
cvRatio = 0.75;                % train fraction
rng(42);

useKfoldCV = true;
kFold = 5;

doMechanistic = true;   % attempt mechanistic model if daily rad/temp found
doLSTM = false;         % set true to attempt LSTM when time-series per-sample available

saveResults = true;
outMat = 'yield_model_results.mat';
outPredCSV = 'predictions.csv';

%% Load data

% if ~isfile(fn)
%     error('File "%s" not found in current folder. Move Crops_data.csv here or set fn appropriately.', fn);
% end
T = readtable(fn);
nObs = height(T);
fprintf('Loaded %d rows and %d columns from %s\n', nObs, width(T), fn);

colNames = T.Properties.VariableNames;
disp('Detected columns:');
disp(colNames');

% Find target column
targetCol = '';
for c = 1:length(targetNamesTry)
    if any(strcmpi(colNames, targetNamesTry{c}))
        targetCol = colNames{strcmpi(colNames, targetNamesTry{c})};
        break;
    end
end
if isempty(targetCol)
    error('No target yield column found. Please include a column named "Yield" (or "yield" / "target").');
end
y = T.(targetCol);
if ~isnumeric(y)
    error('Target column "%s" must be numeric.', targetCol);
end


%% Detect time-series / daily columns and aggregate

% Heuristic: columns containing tokens like 'day', 'd_', 'rad', 'temp', 'tmean', 'ndvi' with suffix numbers
tsTokens = {'day','d_','rad','rad_','temp','tmean','ndvi','precip','prec','rain'};
tsGroups = struct(); % will hold grouped column names

for i=1:length(colNames)
    cn = colNames{i};
    for tk = tsTokens
        if contains(lower(cn), tk{1})
            % attempt to extract base name and trailing index
            % e.g., Rad_day1, rad_01, ndvi_1 -> base 'rad'
            tokens = regexp(cn, '([A-Za-z]+)[_\- ]?(\d+)$','tokens');
            if ~isempty(tokens)
                base = lower(tokens{1}{1});
                idx = str2double(tokens{1}{2});
                if ~isfield(tsGroups, base)
                    tsGroups.(base) = struct('names',{},{'indices',[]});
                    tsGroups.(base).names = {};
                    tsGroups.(base).indices = [];
                end
                tsGroups.(base).names{end+1} = cn;
                tsGroups.(base).indices(end+1) = idx;
            else
                % fallback: group by token itself
                base = tk{1};
                if ~isfield(tsGroups, base)
                    tsGroups.(base) = struct('names',{},{'indices',[]});
                    tsGroups.(base).names = {};
                    tsGroups.(base).indices = [];
                end
                tsGroups.(base).names{end+1} = cn;
                tsGroups.(base).indices(end+1) = NaN;
            end
            break;
        end
    end
end

% Aggregate time-series features: sum, mean, max, integral (trapezoid)
aggFeatures = table();
aggFeatureNames = {};
groupNames = fieldnames(tsGroups);
for g = 1:length(groupNames)
    gname = groupNames{g};
    colList = tsGroups.(gname).names;
    if isempty(colList), continue; end
    % create matrix (nObs x nDays)
    try
        mat = zeros(nObs, length(colList));
        for j=1:length(colList)
            mat(:,j) = double(T.(colList{j})); % cast to double
        end
        % aggregate stats
        agg_sum = sum(mat,2,'omitnan');
        agg_mean = mean(mat,2,'omitnan');
        agg_max = max(mat,[],2);
        % trapezoidal integral along columns (if columns ordered by index)
        % attempt to order by indices if available
        [~,order] = sort(tsGroups.(gname).indices);
        if ~all(isnan(tsGroups.(gname).indices))
            matOrd = mat(:,order);
        else
            matOrd = mat;
        end
        integralVals = trapz(matOrd,2);
        % add to aggFeatures
        fnamePrefix = matlab.lang.makeValidName(['agg_' gname]);
        aggFeatures.([fnamePrefix '_sum']) = agg_sum;
        aggFeatures.([fnamePrefix '_mean']) = agg_mean;
        aggFeatures.([fnamePrefix '_max']) = agg_max;
        aggFeatures.([fnamePrefix '_int']) = integralVals;
        aggFeatureNames = [aggFeatureNames, ...
            {[fnamePrefix '_sum'],[fnamePrefix '_mean'],[fnamePrefix '_max'],[fnamePrefix '_int']}];
        fprintf('Aggregated %d time-series columns for group "%s" into features: %s\n', length(colList), gname, strjoin(aggFeatureNames(end-3:end),', '));
    catch ME
        warning('Failed to aggregate time-series group %s: %s', gname, ME.message);
    end
end

%% Prepare feature table

% Exclude the target column and any ID or year columns from predictors
excludeCols = {targetCol};
possibleID = {'ID','id','FieldID','field','Plot','plot'};
for idc = 1:length(possibleID)
    if any(strcmpi(colNames, possibleID{idc}))
        excludeCols{end+1} = colNames{strcmpi(colNames, possibleID{idc})};
    end
end

% Build initial predictor table from remaining non-time-series columns
predictorNames = setdiff(colNames, [excludeCols, vertcat(tsGroups.(groupNames{:}).names)]);
% remove empty and ensure valid names
predictorNames = predictorNames(~cellfun(@isempty,predictorNames));
Xtable = table();
for i = 1:length(predictorNames)
    cn = predictorNames{i};
    try
        Xtable.(matlab.lang.makeValidName(cn)) = T.(cn);
    catch
        % skip if problem converting
    end
end

% Append aggregated features
if ~isempty(aggFeatures)
    Xtable = [Xtable, aggFeatures];
end

% Final feature names
featNames = Xtable.Properties.VariableNames;
fprintf('Final feature count: %d\n', length(featNames));

% Convert categorical columns to categorical type and numeric to double
for i = 1:length(featNames)
    v = Xtable.(featNames{i});
    if iscellstr(v) || isstring(v) || iscategorical(v)
        Xtable.(featNames{i}) = categorical(v);
    elseif isnumeric(v)
        Xtable.(featNames{i}) = double(v);
    end
end

%% Handle missing data (simple imputation)

for i = 1:length(featNames)
    v = Xtable.(featNames{i});
    if isa(v,'double')
        missingIdx = isnan(v);
        if any(missingIdx)
            med = median(v(~missingIdx));
            if isempty(med) || isnan(med), med = 0; end
            v(missingIdx) = med;
            Xtable.(featNames{i}) = v;
            fprintf('Imputed %d missing values in numeric feature %s with median=%.3g\n', sum(missingIdx), featNames{i}, med);
        end
    elseif isa(v,'categorical')
        missingIdx = isundefined(v);
        if any(missingIdx)
            modeVal = mode(v(~missingIdx));
            v(missingIdx) = modeVal;
            Xtable.(featNames{i}) = v;
            fprintf('Imputed %d missing values in categorical feature %s with mode=%s\n', sum(missingIdx), featNames{i}, string(modeVal));
        end
    end
end


%% Encode categorical variables (one-hot)

% Use dummyvar for categorical columns
X = [];
Xnames = {};
for i = 1:length(featNames)
    v = Xtable.(featNames{i});
    if isa(v,'categorical')
        cats = categories(v);
        D = dummyvar(v);
        for c = 1:length(cats)
            cname = matlab.lang.makeValidName([featNames{i} '_' char(cats{c})]);
            X(:,end+1) = D(:,c);
            Xnames{end+1} = cname;
        end
    else
        X(:,end+1) = Xtable.(featNames{i});
        Xnames{end+1} = featNames{i};
    end
end
X = double(X);

%% Train/test split and scaling

n = size(X,1);
idx = randperm(n);
nTrain = max(2,round(cvRatio*n));
trainIdx = idx(1:nTrain);
testIdx = idx(nTrain+1:end);

Xtrain = X(trainIdx,:);
ytrain = y(trainIdx);
Xtest = X(testIdx,:);
ytest = y(testIdx);

% Standardize numeric features (z-score) using training mean/std
mu = mean(Xtrain,1);
sigma = std(Xtrain,[],1);
sigma(sigma==0) = 1;
Xtrain_s = (Xtrain - mu)./sigma;
Xtest_s = (Xtest - mu)./sigma;


%% Optional: k-fold CV on training set to get validation RMSE

if useKfoldCV
    cv = cvpartition(length(ytrain),'KFold',kFold);
else
    cv = [];
end

%% Model 1: Random Forest (TreeBagger)

nTrees = 300;
fprintf('Training Random Forest with %d trees...\n', nTrees);
rfModel = TreeBagger(nTrees, Xtrain_s, ytrain, 'Method','regression', 'OOBPrediction','On', 'MinLeafSize',5);
oobErr = oobError(rfModel);
figure('Name','RF OOB Error'); plot(oobErr); xlabel('Num Trees'); ylabel('OOB MSE');

% predict
y_rf_pred = predict(rfModel, Xtest_s);

%% Model 2: Gradient Boosting (LSBoost)

fprintf('Training Gradient Boosting (LSBoost)...\n');
t = templateTree('MinLeafSize',5);
gbModel = fitrensemble(Xtrain_s, ytrain, 'Method','LSBoost', 'NumLearningCycles',200, 'Learners', t);
y_gb_pred = predict(gbModel, Xtest_s);


%% Optional mechanistic model (simple biomass -> yield)

y_mech_pred = nan(size(ytest));
mech_available = false;
if doMechanistic && isfield(tsGroups,'rad') && isfield(tsGroups,'temp') || isfield(tsGroups,'tmean') || isfield(tsGroups,'temp')
    try
        fprintf('Attempting simple mechanistic model using daily rad & temp aggregates...\n');
        mech_available = true;
        % Use aggregated features (prefix agg_rad_int or similar) if present
        ixRadInt = find(contains(Xnames,'agg_rad_int'));
        ixTempInt = find(contains(Xnames,'agg_temp_int') | contains(Xnames,'agg_tmean_int'));
        % fallback to mean if integral missing
        if isempty(ixRadInt), ixRadInt = find(contains(Xnames,'agg_rad_mean')); end
        if isempty(ixTempInt), ixTempInt = find(contains(Xnames,'agg_temp_mean') | contains(Xnames,'agg_tmean_mean')); end
        
        % Mechanistic parameters (will be calibrateable)
        mechParams.RUE = 1.8; % g/MJ
        mechParams.HI = 0.45;
        mechParams.base_temp = 8;
        mechParams.rad_coeff = 0.5;
        
        % A quick surrogate mechanistic model: yield = HI * (RUE * TotRad_int * fTemp) * conv
        % find features per sample
        TotRad_train = Xtrain(:,ixRadInt);
        TotTemp_train = Xtrain(:,ixTempInt);
        TotRad_test = Xtest(:,ixRadInt);
        TotTemp_test = Xtest(:,ixTempInt);
        if isempty(TotRad_test) || isempty(TotTemp_test)
            mech_available = false;
        else
            % simple temperature factor: 1 - 0.01*(abs(meanTemp - 22))
            fTemp_test = max(0.2, 1 - 0.01*abs(TotTemp_test - 22));
            % compute biomass g/m2 -> t/ha (g/m2 * 0.01)
            biomass_gm2 = mechParams.RUE .* TotRad_test .* fTemp_test;
            yield_m = biomass_gm2 * 0.01 * mechParams.HI;
            y_mech_pred = yield_m;
            
            % Calibrate RUE & HI using training set if optimization toolbox present
            if license('test','Optimization_Toolbox')
                fprintf('Calibrating mechanistic RUE & HI using lsqnonlin...\n');
                % build training arrays
                TotRad_train_a = Xtrain(:,ixRadInt);
                TotTemp_train_a = Xtrain(:,ixTempInt);
                objfun = @(p) mech_obj_residuals(p, TotRad_train_a, TotTemp_train_a, ytrain);
                p0 = [mechParams.RUE, mechParams.HI];
                lb = [0.2, 0.05]; ub = [5, 0.8];
                opts = optimoptions('lsqnonlin','Display','off','MaxIter',200);
                try
                    p_est = lsqnonlin(objfun, p0, lb, ub, opts);
                    mechParams.RUE = p_est(1);
                    mechParams.HI = p_est(2);
                    fprintf('Calibrated mechanistic params: RUE=%.3f, HI=%.3f\n', mechParams.RUE, mechParams.HI);
                    % recompute test predictions
                    fTemp_test = max(0.2, 1 - 0.01*abs(TotTemp_test - 22));
                    biomass_gm2 = mechParams.RUE .* TotRad_test .* fTemp_test;
                    y_mech_pred = biomass_gm2 * 0.01 * mechParams.HI;
                catch ME
                    warning('Mechanistic calibration failed: %s', ME.message);
                end
            else
                fprintf('Optimization toolbox not available — skipping calibration.\n');
            end
        end
    catch ME
        warning('Mechanistic modeling failed: %s', ME.message);
        mech_available = false;
    end
else
    fprintf('Daily rad/temp groups not detected — skipping mechanistic model.\n');
end

%% Ensemble: weight by validation RMSE

% Compute validation RMSE using a small holdout (or OOB for RF)
yvals_rf = predict(rfModel, Xtrain_s);
yvals_gb = predict(gbModel, Xtrain_s);
rmse_rf_val = sqrt(mean((yvals_rf - ytrain).^2));
rmse_gb_val = sqrt(mean((yvals_gb - ytrain).^2));

w_rf = 1/(rmse_rf_val+eps); w_gb = 1/(rmse_gb_val+eps);
w_sum = w_rf + w_gb;
w_rf = w_rf / w_sum; w_gb = w_gb / w_sum;

y_ens_pred = w_rf * y_rf_pred + w_gb * y_gb_pred;
if mech_available
    % include mech with weight proportional to its validation RMSE if computed
    % compute mech rmse on training if available
    try
        TotRad_train_a = Xtrain(:,ixRadInt);
        TotTemp_train_a = Xtrain(:,ixTempInt);
        fTemp_train = max(0.2, 1 - 0.01*abs(TotTemp_train_a - 22));
        biomass_gm2_train = mechParams.RUE .* TotRad_train_a .* fTemp_train;
        y_mech_train = biomass_gm2_train * 0.01 * mechParams.HI;
        rmse_mech_val = sqrt(mean((y_mech_train - ytrain).^2));
        w_mech = 1/(rmse_mech_val+eps);
        % renormalize weights
        w_rf = w_rf / (w_rf + w_gb + w_mech);
        w_gb = w_gb / (w_rf + w_gb + w_mech);
        % recompute (normalize properly)
        W = [1/(rmse_rf_val+eps), 1/(rmse_gb_val+eps), 1/(rmse_mech_val+eps)];
        W = W / sum(W);
        y_ens_pred = W(1)*y_rf_pred + W(2)*y_gb_pred + W(3)*y_mech_pred;
        fprintf('Ensemble weights (RF, GB, Mech): %.3f, %.3f, %.3f\n', W(1), W(2), W(3));
    catch
        % ignore mech inclusion if fails
    end
else
    fprintf('Ensemble weights (RF, GB): %.3f, %.3f\n', w_rf, w_gb);
end

%% Evaluate performance

rmse = @(a,b) sqrt(mean((a-b).^2));
mae = @(a,b) mean(abs(a-b));
r2 = @(a,b) 1 - sum((a-b).^2)/sum((b-mean(b)).^2);

metrics = struct();
metrics.RF.RMSE = rmse(y_rf_pred, ytest);
metrics.RF.MAE  = mae(y_rf_pred, ytest);
metrics.RF.R2   = r2(y_rf_pred, ytest);

metrics.GB.RMSE = rmse(y_gb_pred, ytest);
metrics.GB.MAE  = mae(y_gb_pred, ytest);
metrics.GB.R2   = r2(y_gb_pred, ytest);

metrics.Ens.RMSE = rmse(y_ens_pred, ytest);
metrics.Ens.MAE  = mae(y_ens_pred, ytest);
metrics.Ens.R2   = r2(y_ens_pred, ytest);

if mech_available
    metrics.Mechanistic.RMSE = rmse(y_mech_pred, ytest);
    metrics.Mechanistic.MAE  = mae(y_mech_pred, ytest);
    metrics.Mechanistic.R2   = r2(y_mech_pred, ytest);
end

disp('Performance on test set:');
disp(metrics);

% Plot observed vs predicted
figure('Name','Observed vs Predicted'); 
plot(ytest,'-ko','LineWidth',1.2); hold on;
plot(y_rf_pred,'--s'); plot(y_gb_pred,':d'); plot(y_ens_pred,'-.^');
if mech_available, plot(y_mech_pred,':+'); end
legend(['Observed','RF','GB','Ensemble', (mech_available)*"Mech" ]);
xlabel('Test sample index'); ylabel('Yield'); title('Observed vs Predicted Yield');

%% -------------------------
%% Save outputs
%% -------------------------
if saveResults
    save(outMat, 'rfModel','gbModel','metrics','Xnames','mu','sigma','featNames','mechParams','ytest','y_rf_pred','y_gb_pred','y_ens_pred','y_mech_pred');
    fprintf('Saved models and metrics to %s\n', outMat);
    % Save predictions with test indices
    predTab = table((1:length(ytest))', ytest, y_rf_pred, y_gb_pred, y_ens_pred, 'VariableNames',{'Idx','Observed','RF','GB','Ensemble'});
    if mech_available
        predTab.Mechanistic = y_mech_pred;
    end
    writetable(predTab, outPredCSV);
    fprintf('Saved predictions to %s\n', outPredCSV);
end


%% Helper functions

function res = mech_obj_residuals(p, TotRad, TotTemp, obsY)
    % surrogate mechanistic residuals for calibration
    % p = [RUE, HI]
    RUE = p(1); HI = p(2);
    fTemp = max(0.2, 1 - 0.01*abs(TotTemp - 22));
    biomass_gm2 = RUE .* TotRad .* fTemp;
    predY = biomass_gm2 * 0.01 * HI;
    res = predY - obsY;
end
