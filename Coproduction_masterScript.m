%to study co-production
initCobraToolbox
solver = 'gurobi';
solverOK = changeCobraSolver(solver,'all');

%%choose organism and conditions

%%ecoli as organism
%model = readCbModel('iML1515.xml'); %for GSMM
model = readCbModel('e_coli_core.xml'); %core model

%%yeast as organism
% model = readCbModel('iMM904.xml');
% %Anaerobic condition %anaerobic growth of yeast requires fatty acid supplementation in medium
% model = changeRxnBounds(model,{'EX_ergst_e';'EX_zymst_e';'EX_hdcea_e';'EX_ocdca_e';'EX_ocdcea_e'; 'EX_ocdcya_e' },-1000,'l');
% model = changeRxnBounds(model,{'EX_ergst_e';'EX_zymst_e';'EX_hdcea_e';'EX_ocdca_e';'EX_ocdcea_e'; 'EX_ocdcya_e' },1000,'u');

%%Anaerobic condition
% model = changeRxnBounds(model,{'EX_o2_e'},0,'b');

%% to study co-production for all pairs of metabolites in a model
minBM = 0.25; %minimum biomass of mutant - given in percentage of wild-type biomass
modelSol = optimizeCbModel(model);
coFseofScoreTable = coFSEOF(model,minBM,solver);

%%o study co-production for a single set of metabolites

%TargetProducts = {'EX_ibutoh_e';'EX_succ_e'}; %example - for iMM904 model
TargetProducts = {'EX_etoh_e';'EX_succ_e'}; %example - for e_coli_core model or iML1515 model
TargetsScoreTable= coFSEOFTargets(model,minBM,solver,TargetProducts);
