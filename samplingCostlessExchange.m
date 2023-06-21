function [outputStruct] = samplingCostlessExchange(model1, model2,EnsembleThreshold, nSamples)
%% Load Models
% modelNames = fieldnames(models);
% nS = numSpecies;
%speciesPairCombos = repelem(nchoosek(modelNames,nS),1,1);
iterations = nSamples;
%% Sample Models alone: Using polytope sampler
% initSampler
% P = struct;
% P.lb = model1.lb;
% P.ub = model1.ub;
% P.beq = model1.b;
% P.Aeq = model1.S;
% Model1_samplingResult = sample(P, iterations);
% Model1_samples = Model1_samplingResult.samples;
% 
% P = struct;
% P.lb = model2.lb;
% P.ub = model2.ub;
% P.beq = model2.b;
% P.Aeq = model2.S;
% Model2_samplingResult = sample(P, iterations);
% Model2_samples = Model2_samplingResult.samples;
% 
% % Save Samples
% outputStruct.Alone.fluxesModel1 = Model1_samples;
% outputStruct.Alone.fluxesModel2 = Model2_samples;

%% Sample Models alone: Using RHMC
Samplingoptions.nStepsPerPoint = 100; %sampling density
Samplingoptions.nPointsReturned = iterations; %number of points returned
Samplingoptions.toRound = 0; %whether or not the polytope is rounded
Samplingoptions.optPercentage = 0;
[~, Model1_samples] =  sampleCbModel(model1, [], 'RHMC', Samplingoptions);
[~, Model2_samples] =  sampleCbModel(model2, [], 'RHMC', Samplingoptions);

% % Save Samples
outputStruct.Alone.fluxesModel1 = Model1_samples;
outputStruct.Alone.fluxesModel2 = Model2_samples;

%% Get secreted metabolites alone
SecByModel1 ={};
SecByModel2 ={};
minSamples = min(size(Model1_samples, 2), size(Model2_samples, 2));
for i = 1:minSamples
    a = Model1_samples(:,i);
    SecByModel1{i} = model1.rxns(intersect(find(a > 1e-6), find(strncmp('EX_',model1.rxns,3))));
    b = Model2_samples(:,i);
    SecByModel2{i} = model2.rxns(intersect(find(b > 1e-6), find(strncmp('EX_',model2.rxns,3))));
    
end
% Exchange Reactions that are turned on (secreting metabolites)
outputStruct.Alone.SecByModel1 = SecByModel1;
outputStruct.Alone.SecByModel2 = SecByModel2;

%% Ensemble Thresholds
if EnsembleThreshold == 'any'
    SecRxnsModel1 = union_several(SecByModel1{1,:});
    SecRxnsModel2 = union_several(SecByModel2{1,:});
    SecRxns = union(SecRxnsModel1, SecRxnsModel2);
elseif EnsembleThreshold == 'all'
    SecRxnsModel1 = mintersect(SecByM1{1,:});
    SecRxnsModel2 = mintersect(SecByM2{1,:});
    SecRxns = union(SecRxnsModel1, SecRxnsModel2);
else
    fprintf('Please select correct ensemble threshold')
end
% Consensus Exchange Reactions that are turned on (secreting metabolites)

outputStruct.Alone.EnsembleSecByModel1 = SecRxnsModel1;
outputStruct.Alone.EnsembleSecByModel2 = SecRxnsModel2;
outputStruct.Alone.EnsembleSecMets = SecRxns;
outputStruct.Alone.Model1 = model1;
outputStruct.Alone.Model2 = model2;

expansions = 1;
newSecretions = SecRxns;
length(newSecretions)
everSecreted = newSecretions;
%% Update Media Conditions
while ~isempty(newSecretions)
    
    % Update Models: set new extracellular environment
    C    = cell(1,size(SecRxns, 1));
    D    = cell(1,size(SecRxns, 1));
    C(:) = {'l'};
    D(:) = {'u'};
    model1 = changeRxnBounds(model1, SecRxns, -1000, C);
    model1 = changeRxnBounds(model1, SecRxns, 1000, D);
    model2 = changeRxnBounds(model2, SecRxns, -1000, C);
    model2 = changeRxnBounds(model2, SecRxns, 1000, D);
    
    expansions = expansions + 1;
    round = strcat('round', num2str(expansions));
    
    [ round ]
    
    % Save Models
    name = 'model1';
    outputStruct.(round).(name) = model1;
    name = 'model2';
    outputStruct.(round).(name) = model2;
    
    % Sample: polytope sampler 
%     initSampler
%     fprintf('sampling model 1')
%     P = struct;
%     P.lb = model1.lb;
%     P.ub = model1.ub;
%     P.beq = model1.b;
%     P.Aeq = model1.S;
%     Model1_samplingResult = sample(P, iterations);
%     Model1_samples = Model1_samplingResult.samples;
%     
%     P = struct;
%     fprintf('sampling model 2')
%     P.lb = model2.lb;
%     P.ub = model2.ub;
%     P.beq = model2.b;
%     P.Aeq = model2.S;
%     Model2_samplingResult = sample(P, iterations);
%     Model2_samples = Model2_samplingResult.samples;
    
% Sample: RHMC
    [~, Model1_samples] =  sampleCbModel(model1, [], 'RHMC', Samplingoptions);
    [~, Model2_samples] =  sampleCbModel(model2, [], 'RHMC', Samplingoptions);
    
    % Save Samples
    name = 'fluxesModel1';
    outputStruct.(round).(name) = Model1_samples;
    name = 'fluxesModel2';
    outputStruct.(round).(name) = Model2_samples;
    
    SecByModel1 ={};
    SecByModel2 ={};
    minSamples = min(size(Model1_samples, 2), size(Model2_samples, 2));
    for i = 1:minSamples
        a = Model1_samples(:,i);
        SecByModel1{i} = model1.rxns(intersect(find(a > 1e-6), find(strncmp('EX_',model1.rxns,3))));
        b = Model2_samples(:,i);
        SecByModel2{i} = model2.rxns(intersect(find(b > 1e-6), find(strncmp('EX_',model2.rxns,3))));
        
    end
    % Exchange Reactions that are turned on (secreting metabolites)
    name = 'SecByModel1';
    outputStruct.(round).(name) = SecByModel1;
    name = 'SecByModel2';
    outputStruct.(round).(name) = SecByModel2;
    
    % Ensemble Thresholds
    if EnsembleThreshold == 'any'
        SecRxnsModel1 = union_several(SecByModel1{1,:});
        SecRxnsModel2 = union_several(SecByModel2{1,:});
        SecRxns = union(SecRxnsModel1, SecRxnsModel2);
    elseif EnsembleThreshold == 'all'
        SecRxnsModel1 = mintersect(SecByM1{1,:});
        SecRxnsModel2 = mintersect(SecByM2{1,:});
        SecRxns = union(SecRxnsModel1, SecRxnsModel2);
    else
        fprintf('Please select correct ensemble threshold')
    end
    % Consensus Exchange Reactions that are turned on (secreting metabolites)
    name = 'SecRxnsModel1';
    outputStruct.(round).(name) = SecRxnsModel1;
    name = 'SecRxnsModel2';
    outputStruct.(round).(name) = SecRxnsModel2;
    name = 'SecRxns';
    outputStruct.(round).(name) = SecRxns;
    
    newSecretions = setdiff(everSecreted, SecRxns);
    length(newSecretions)
    everSecreted = unique(union(everSecreted, newSecretions));
    
    save('outputStruct.mat', 'outputStruct');
end



end