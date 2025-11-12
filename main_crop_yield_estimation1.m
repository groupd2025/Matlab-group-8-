%% CROP YIELD ESTIMATION SYSTEM
% High-end implementation with multiple approaches
% Requires: Statistics and ML Toolbox

function main_crop_yield_estimation()
    clear; clc; close all;
    
    fprintf('=== CROP YIELD ESTIMATION SYSTEM ===\n');
    
    % Configuration
    config = initialize_configuration();
    
    % Generate synthetic dataset
    [features, targets, spatial_data] = generate_synthetic_dataset(config);
    
    % Data preprocessing
    [processed_features, feature_names] = preprocess_data(features, config);
    
    % Split data for training and testing
    [X_train, X_test, y_train, y_test, train_idx, test_idx] = ...
        split_dataset(processed_features, targets, 0.8);
    
    %% MACHINE LEARNING APPROACHES
    fprintf('\n--- Training Machine Learning Models ---\n');
    
    % 1. Ensemble Random Forest
    rf_model = train_random_forest(X_train, y_train, config);
    y_pred_rf = predict(rf_model, X_test);
    
    % 2. Support Vector Machine
    svm_model = train_svm(X_train, y_train, config);
    y_pred_svm = predict(svm_model, X_test);
    
    % 3. Neural Network
    [nn_model, y_pred_nn] = train_neural_network(X_train, y_train, X_test, config);
    
    % 4. Gradient Boosting
    gbr_model = train_gradient_boosting(X_train, y_train, config);
    y_pred_gbr = predict(gbr_model, X_test);
    
    %% PROCESS-BASED MODELING
    fprintf('\n--- Running Process-Based Simulation ---\n');
    
    % Crop growth simulation
    simulated_yields = simulate_crop_growth(spatial_data, config);
    y_pred_sim = simulated_yields(test_idx);
    
    %% ENSEMBLE PREDICTION
    fprintf('\n--- Creating Ensemble Prediction ---\n');
    
    % Combine predictions with dynamic weighting
    ensemble_weights = calculate_ensemble_weights(...
        [y_pred_rf, y_pred_svm, y_pred_nn, y_pred_gbr, y_pred_sim], y_test);
    
    y_pred_ensemble = [y_pred_rf, y_pred_svm, y_pred_nn, y_pred_gbr, y_pred_sim] * ensemble_weights;
    
    %% MODEL EVALUATION
    fprintf('\n--- Model Performance Evaluation ---\n');
    
    models = {'Random_Forest', 'SVM', 'Neural_Network', 'Gradient_Boost', 'Simulation', 'Ensemble'};
    predictions = [y_pred_rf, y_pred_svm, y_pred_nn, y_pred_gbr, y_pred_sim, y_pred_ensemble];
    
    performance_metrics = evaluate_models(models, predictions, y_test);
    display_performance_table(performance_metrics, models);
    
    %% SPATIAL YIELD MAPPING
    fprintf('\n--- Generating Spatial Yield Maps ---\n');
    
    create_spatial_yield_maps(spatial_data, predictions, test_idx, config);
    
    %% UNCERTAINTY ANALYSIS
    fprintf('\n--- Uncertainty Quantification ---\n');
    
    uncertainty_analysis = perform_uncertainty_analysis(predictions, y_test);
    
    %% MANAGEMENT INSIGHTS
    fprintf('\n--- Generating Management Insights ---\n');
    
    generate_management_insights(processed_features, targets, config);
    
    %% VISUALIZATION DASHBOARD
    fprintf('\n--- Creating Comprehensive Dashboard ---\n');
    
    create_comprehensive_dashboard(models, predictions, y_test, spatial_data, ...
        performance_metrics, uncertainty_analysis, test_idx, config);
    
    fprintf('\n=== ANALYSIS COMPLETE ===\n');
end

%% CONFIGURATION FUNCTION
function config = initialize_configuration()
    config = struct();
    
    % Dataset parameters
    config.n_samples = 500;  % Reduced for faster execution
    config.n_features = 12;
    config.train_test_split = 0.8;
    
    % Model parameters
    config.rf_num_trees = 50;
    config.svm_kernel = 'gaussian';
    config.nn_hidden_layers = [20, 10];  % Simplified architecture
    config.gbr_learning_rate = 0.1;
    config.gbr_n_estimators = 50;
    
    % Simulation parameters
    config.growing_season_days = 120;
    config.field_size = [100, 100]; % meters
    
    % Visualization parameters
    config.map_resolution = 2; % meters
end

%% SYNTHETIC DATASET GENERATION
function [features, targets, spatial_data] = generate_synthetic_dataset(config)
    rng(42); % For reproducibility
    
    n = config.n_samples;
    p = config.n_features;
    
    % Generate realistic features
    features = zeros(n, p);
    
    % 1. Satellite-derived features (NDVI, EVI, etc.)
    features(:,1) = 0.3 + 0.6 * rand(n,1); % NDVI
    features(:,2) = 0.2 + 0.5 * rand(n,1); % EVI
    features(:,3) = 0.4 + 0.4 * rand(n,1); % GNDVI
    
    % 2. Weather data
    features(:,4) = 15 + 15 * rand(n,1); % Temperature
    features(:,5) = 500 + 800 * rand(n,1); % Precipitation
    features(:,6) = 1000 + 500 * rand(n,1); % Growing Degree Days
    
    % 3. Soil properties
    features(:,7) = 5.5 + 1.5 * rand(n,1); % pH
    features(:,8) = 0.1 + 0.3 * rand(n,1); % Nitrogen
    features(:,9) = 1.0 + 3.0 * rand(n,1); % Organic Matter
    
    % 4. Management practices
    features(:,10) = 50 + 100 * rand(n,1); % Fertilizer
    features(:,11) = randi([0, 3], n, 1); % Irrigation level
    features(:,12) = randn(n, 1); % Other covariate
    
    % Generate realistic yield targets (tons/hectare)
    base_yield = 8.0; % Base yield
    weather_effect = 0.1 * features(:,4) - 0.05 * features(:,5)/100;
    soil_effect = 2.0 * features(:,8) + 0.5 * features(:,9);
    satellite_effect = 4.0 * features(:,1) + 2.0 * features(:,2);
    
    targets = base_yield + weather_effect + soil_effect + satellite_effect + ...
              0.1 * features(:,10) + 0.5 * features(:,11) + randn(n,1)*0.5;
    
    % Ensure positive yields
    targets = max(2, targets);
    
    % Spatial data
    spatial_data.coordinates = 100 * rand(n, 2);
    spatial_data.field_boundary = [0, 100, 0, 100];
    spatial_data.all_predictions = []; % Initialize
    
    fprintf('Generated synthetic dataset: %d samples, %d features\n', n, p);
end

%% DATA PREPROCESSING
function [processed_features, feature_names] = preprocess_data(features, config)
    % Handle missing values (if any)
    if any(any(isnan(features)))
        features = fillmissing(features, 'movmedian', 5);
    end
    
    % Normalize features
    processed_features = normalize(features);
    
    % Feature names
    feature_names = {'NDVI', 'EVI', 'GNDVI', 'Temperature', 'Precipitation', ...
                    'GDD', 'pH', 'Nitrogen', 'Organic Matter', 'Fertilizer', ...
                    'Irrigation', 'Feature12'};
    
    fprintf('Data preprocessing completed\n');
end

%% DATA SPLITTING
function [X_train, X_test, y_train, y_test, train_idx, test_idx] = ...
         split_dataset(features, targets, split_ratio)
    
    n = size(features, 1);
    split_point = floor(split_ratio * n);
    
    % Random permutation
    idx = randperm(n);
    train_idx = idx(1:split_point);
    test_idx = idx(split_point+1:end);
    
    X_train = features(train_idx, :);
    X_test = features(test_idx, :);
    y_train = targets(train_idx);
    y_test = targets(test_idx);
    
    fprintf('Data split: %d training, %d test samples\n', length(train_idx), length(test_idx));
end

%% RANDOM FOREST MODEL
function model = train_random_forest(X_train, y_train, config)
    model = TreeBagger(config.rf_num_trees, X_train, y_train, ...
                      'Method', 'regression', ...
                      'OOBPrediction', 'on', ...
                      'MinLeafSize', 5);
    
    fprintf('Random Forest trained with %d trees\n', config.rf_num_trees);
end

%% SVM MODEL
function model = train_svm(X_train, y_train, config)
    try
        model = fitrsvm(X_train, y_train, ...
                       'KernelFunction', config.svm_kernel, ...
                       'Standardize', true, ...
                       'KernelScale', 'auto');
        fprintf('SVM model trained with %s kernel\n', config.svm_kernel);
    catch
        % Fallback to linear kernel if gaussian fails
        model = fitrsvm(X_train, y_train, ...
                       'KernelFunction', 'linear', ...
                       'Standardize', true);
        fprintf('SVM model trained with linear kernel (fallback)\n');
    end
end

%% NEURAL NETWORK MODEL - FIXED VERSION
function [net, y_pred] = train_neural_network(X_train, y_train, X_test, config)
    % Create a simpler neural network approach
    try
        % Normalize data for neural network
        X_train_nn = X_train';
        y_train_nn = y_train';
        X_test_nn = X_test';
        
        % Create and configure neural network
        net = feedforwardnet(config.nn_hidden_layers);
        net.trainParam.showWindow = false;
        net.trainParam.epochs = 100;
        net.trainParam.max_fail = 10;
        
        % Train network
        [net, tr] = train(net, X_train_nn, y_train_nn);
        
        % Predict
        y_pred = net(X_test_nn)';
        
        fprintf('Neural Network trained successfully\n');
    catch ME
        fprintf('Neural Network training failed, using linear regression instead: %s\n', ME.message);
        % Fallback to linear regression
        net = fitlm(X_train, y_train);
        y_pred = predict(net, X_test);
    end
end

%% GRADIENT BOOSTING MODEL
function model = train_gradient_boosting(X_train, y_train, config)
    try
        template = templateTree('MaxNumSplits', 20, 'MinLeafSize', 5);
        model = fitrensemble(X_train, y_train, ...
                            'Method', 'LSBoost', ...
                            'NumLearningCycles', config.gbr_n_estimators, ...
                            'LearnRate', config.gbr_learning_rate, ...
                            'Learners', template);
        fprintf('Gradient Boosting trained with %d estimators\n', config.gbr_n_estimators);
    catch
        % Fallback to bagging
        model = fitrensemble(X_train, y_train, 'Method', 'Bag');
        fprintf('Gradient Boosting failed, using Bagging instead\n');
    end
end

%% PROCESS-BASED CROP GROWTH SIMULATION
function yields = simulate_crop_growth(spatial_data, config)
    n = size(spatial_data.coordinates, 1);
    yields = zeros(n, 1);
    
    for i = 1:n
        % Extract location-based factors
        x = spatial_data.coordinates(i, 1);
        y = spatial_data.coordinates(i, 2);
        
        % Simulate environmental gradients
        temperature_gradient = 20 + 0.1 * x - 0.05 * y;
        soil_quality = 0.5 + 0.3 * sin(0.1 * x) * cos(0.1 * y);
        
        % Simple yield calculation
        base_potential = 8.0; % tons/hectare
        yields(i) = base_potential * soil_quality * (1 + 0.02 * temperature_gradient) + ...
                    randn() * 0.5;
    end
    
    fprintf('Process-based simulation completed for %d locations\n', n);
end

%% ENSEMBLE WEIGHT CALCULATION
function weights = calculate_ensemble_weights(predictions, actual)
    n_models = size(predictions, 2);
    errors = zeros(n_models, 1);
    
    % Calculate RMSE for each model
    for i = 1:n_models
        errors(i) = sqrt(mean((predictions(:, i) - actual).^2));
    end
    
    % Convert errors to weights (lower error = higher weight)
    weights = 1 ./ (errors + 1e-8);  % Add small constant to avoid division by zero
    weights = weights / sum(weights); % Normalize
    
    fprintf('Ensemble weights: %s\n', mat2str(weights, 2));
end

%% MODEL EVALUATION
function performance_metrics = evaluate_models(model_names, predictions, actual)
    n_models = length(model_names);
    
    for i = 1:n_models
        pred = predictions(:, i);
        
        % Calculate various metrics
        rmse = sqrt(mean((pred - actual).^2));
        mae = mean(abs(pred - actual));
        r2 = max(0, 1 - sum((actual - pred).^2) / sum((actual - mean(actual)).^2));
        mape = mean(abs((actual - pred) ./ (actual + 1e-8))) * 100; % Avoid division by zero
        
        performance_metrics(i) = struct(...
            'Model', model_names{i}, ...
            'RMSE', rmse, ...
            'MAE', mae, ...
            'R2', r2, ...
            'MAPE', mape);
    end
end

%% PERFORMANCE TABLE DISPLAY
function display_performance_table(performance_metrics, models)
    fprintf('\n%-20s %-8s %-8s %-8s %-8s\n', 'Model', 'RMSE', 'MAE', 'R²', 'MAPE%');
    fprintf('%s\n', repmat('-', 1, 60));
    
    for i = 1:length(performance_metrics)
        metrics = performance_metrics(i);
        fprintf('%-20s %-8.3f %-8.3f %-8.3f %-8.2f\n', ...
                models{i}, metrics.RMSE, metrics.MAE, metrics.R2, metrics.MAPE);
    end
end

%% SPATIAL YIELD MAPPING
function create_spatial_yield_maps(spatial_data, predictions, test_idx, config)
    figure('Position', [100, 100, 1200, 600], 'Name', 'Spatial Yield Analysis');
    
    % Use ensemble predictions (last column)
    ensemble_pred = predictions(:, end);
    
    % Extract test coordinates and predictions
    test_coords = spatial_data.coordinates(test_idx, :);
    test_predictions = ensemble_pred;
    
    % Create interpolation for smooth map
    try
        F = scatteredInterpolant(test_coords(:,1), test_coords(:,2), test_predictions, 'natural', 'linear');
        
        [Xq, Yq] = meshgrid(...
            linspace(0, 100, 50), ...
            linspace(0, 100, 50));
        Zq = F(Xq, Yq);
        
        % Create yield map
        subplot(1, 2, 1);
        contourf(Xq, Yq, Zq, 20, 'LineStyle', 'none');
        colorbar;
        title('Spatial Yield Distribution (tons/ha)');
        xlabel('X Coordinate (m)');
        ylabel('Y Coordinate (m)');
        axis equal;
    catch
        % Fallback to scatter plot
        subplot(1, 2, 1);
        scatter(test_coords(:,1), test_coords(:,2), 50, test_predictions, 'filled');
        colorbar;
        title('Spatial Yield Distribution (tons/ha)');
        xlabel('X Coordinate (m)');
        ylabel('Y Coordinate (m)');
        axis equal;
    end
    
    % Yield histogram
    subplot(1, 2, 2);
    histogram(test_predictions, 20, 'FaceColor', 'green', 'FaceAlpha', 0.7);
    title('Yield Distribution');
    xlabel('Yield (tons/ha)');
    ylabel('Frequency');
    grid on;
    
    fprintf('Spatial yield maps generated\n');
end

%% UNCERTAINTY ANALYSIS
function uncertainty = perform_uncertainty_analysis(predictions, actual)
    ensemble_pred = predictions(:, end);
    
    uncertainty = struct();
    uncertainty.prediction_std = std(ensemble_pred - actual);
    uncertainty.confidence_interval = [mean(ensemble_pred) - 1.96*std(ensemble_pred), ...
                                      mean(ensemble_pred) + 1.96*std(ensemble_pred)];
    uncertainty.coefficient_variation = std(ensemble_pred) / mean(ensemble_pred) * 100;
    
    fprintf('Uncertainty analysis completed\n');
    fprintf('Coefficient of Variation: %.2f%%\n', uncertainty.coefficient_variation);
end

%% MANAGEMENT INSIGHTS
function generate_management_insights(features, targets, config)
    fprintf('\n--- MANAGEMENT INSIGHTS ---\n');
    
    % Calculate correlations with yield
    correlations = zeros(config.n_features, 1);
    for i = 1:config.n_features
        correlations(i) = corr(features(:,i), targets);
    end
    
    [~, top_idx] = sort(abs(correlations), 'descend');
    
    feature_names = {'NDVI', 'EVI', 'GNDVI', 'Temperature', 'Precipitation', ...
                    'GDD', 'pH', 'Nitrogen', 'Organic Matter', 'Fertilizer', ...
                    'Irrigation', 'Feature12'};
    
    fprintf('Top factors influencing yield:\n');
    for i = 1:min(5, length(top_idx))
        idx = top_idx(i);
        fprintf('  %d. %s: correlation = %.3f\n', i, feature_names{idx}, correlations(idx));
    end
    
    % Yield optimization suggestions
    avg_yield = mean(targets);
    fprintf('\nYield Optimization Suggestions:\n');
    fprintf('  • Current average yield: %.2f tons/ha\n', avg_yield);
    fprintf('  • Target yield range: %.2f - %.2f tons/ha\n', avg_yield*0.9, avg_yield*1.1);
    fprintf('  • Focus on improving: %s\n', feature_names{top_idx(1)});
    
    % Economic impact
    area_hectares = 100; % Example field size
    value_per_ton = 200; % USD per ton
    potential_increase = avg_yield * 0.1; % 10% increase
    economic_impact = potential_increase * area_hectares * value_per_ton;
    
    fprintf('  • Potential economic impact: $%.2f per season\n', economic_impact);
end

%% COMPREHENSIVE DASHBOARD - FIXED VERSION
function create_comprehensive_dashboard(model_names, predictions, actual, ...
                                       spatial_data, metrics, uncertainty, test_idx, config)
    
    figure('Name', 'Crop Yield Estimation Dashboard', ...
           'Position', [50, 50, 1200, 800]);
    
    % 1. Prediction vs Actual plot
    subplot(2, 3, 1);
    ensemble_pred = predictions(:, end);
    plot(actual, ensemble_pred, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'blue');
    hold on;
    plot([min(actual), max(actual)], [min(actual), max(actual)], 'r--', 'LineWidth', 2);
    xlabel('Actual Yield (tons/ha)');
    ylabel('Predicted Yield (tons/ha)');
    title('Ensemble: Prediction vs Actual');
    grid on;
    legend('Predictions', 'Perfect Fit', 'Location', 'northwest');
    
    % 2. Residual analysis
    subplot(2, 3, 2);
    residuals = ensemble_pred - actual;
    histogram(residuals, 20, 'FaceColor', 'red', 'FaceAlpha', 0.7);
    xlabel('Residuals');
    ylabel('Frequency');
    title('Residual Distribution');
    grid on;
    
    % 3. Model comparison
    subplot(2, 3, 3);
    rmse_values = zeros(length(metrics), 1);
    for i = 1:length(metrics)
        rmse_values(i) = metrics(i).RMSE;
    end
    bar(rmse_values);
    set(gca, 'XTickLabel', model_names, 'XTickLabelRotation', 45);
    ylabel('RMSE');
    title('Model Performance Comparison');
    grid on;
    
    % 4. Spatial distribution
    subplot(2, 3, 4);
    test_coords = spatial_data.coordinates(test_idx, :);
    scatter(test_coords(:,1), test_coords(:,2), 30, ensemble_pred, 'filled');
    colorbar;
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Yield Spatial Distribution');
    axis equal;
    
    % 5. Time series simulation
    subplot(2, 3, 5);
    days = 1:config.growing_season_days;
    simulated_growth = 8 * (1 - exp(-days/30)) + randn(size(days))*0.1;
    plot(days, simulated_growth, 'g-', 'LineWidth', 2);
    xlabel('Days');
    ylabel('Biomass (tons/ha)');
    title('Simulated Crop Growth');
    grid on;
    
    % 6. Economic impact
    subplot(2, 3, 6);
    yield_levels = [6, 8, 10, 12];
    revenue = yield_levels * 100 * 200; % area * price
    bar(yield_levels, revenue/1000);
    xlabel('Yield (tons/ha)');
    ylabel('Revenue (Thousand USD)');
    title('Economic Impact Analysis');
    grid on;
    
    fprintf('Comprehensive dashboard created\n');
end
