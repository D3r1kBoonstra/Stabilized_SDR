addpath(genpath(['.' filesep]));
data = readmatrix('ends_sim/data.csv', 'Delimiter', ',');
u = readmatrix('ends_sim/u.csv', 'Delimiter', ',');
ncores = readmatrix('ends_sim/ncores.csv', 'Delimiter', ',');
nsims = readmatrix('ends_sim/nsims.csv', 'Delimiter', ',');
rng(1);
lam = ENDS_CV(data(:,2:end),data(:,1), u, 0:0.01:1);
err = zeros(nsims, 1);

parpool(ncores, 'IdleTimeout', 180)
parfor sim = 1:nsims

    rng(sim);
    X = data(:, 2:end); 
    y = data(:, 1); 
    
    % Create 10-fold cross-validation partition
    cv = cvpartition(y, 'KFold', 10, 'Stratify', true);
    
    % Initialize an array to store validation errors
    validationErrors = zeros(cv.NumTestSets, 1);
    
        for i = 1:cv.NumTestSets
            % Get training and test indices for the current fold
            trainIdx = training(cv, i);
            testIdx = test(cv, i);
            
            % Extract training and test data
            XTrain = X(trainIdx, :);
            yTrain = y(trainIdx, :);
            XTest = X(testIdx, :);
            yTest = y(testIdx, :);
            
            % Train your model (for example, a simple linear model)
            Ghat = ENDS(XTrain, yTrain, u, lam.qda);
            qda = fitcdiscr(XTrain*Ghat, yTrain,'DiscrimType','PseudoQuadratic');
            yhat = predict(qda,XTest*Ghat);
            
            % Compute validation error for the current fold
            validationErrors(i) = mean(yhat ~= yTest);
        end
    
    % Compute average validation error
    err(sim) = mean(validationErrors);
end
delete(gcp('nocreate'));
csvwrite("ends_sim/ENDS_error.csv", err);

