%to study co-production
initCobraToolbox
solverOK = changeCobraSolver('ibm_cplex','all');

%%choose organism and conditions

% %ecoli as organism
%model = readCbModel('iML1515.xml'); %for GSMM
model = readCbModel('e_coli_core.xml'); %core model

%%Anaerobic condition
% model = changeRxnBounds(model,{'EX_o2_e'},0,'b');

%%yeast as organism
% model = readCbModel('iMM904.xml');
% %Anaerobic condition %anaerobic growth of yeast requires fatty acid supplementation in medium
% model = changeRxnBounds(model,{'EX_ergst_e';'EX_zymst_e';'EX_hdcea_e';'EX_ocdca_e';'EX_ocdcea_e'; 'EX_ocdcya_e' },-1000,'l');
% model = changeRxnBounds(model,{'EX_ergst_e';'EX_zymst_e';'EX_hdcea_e';'EX_ocdca_e';'EX_ocdcea_e'; 'EX_ocdcya_e' },1000,'u');
% model = changeRxnBounds(model,{'EX_o2_e'},0,'b');

%% to study co-production for all pairs of metabolites in a model
modelSol = optimizeCbModel(model);
[FseofAll,multiFseof,Amp,KO] = multiFSEOF(model);
multiFseofScoreTable = multiFseof_Score(model,multiFseof,Amp,KO);

%% to study co-production for a single set of metabolites

%TargetProducts = {'EX_ibutoh_e';'EX_succ_e'}; %example - for iML1515 model
TargetProducts = {'EX_etoh_e';'EX_succ_e'}; %example - for e_coli_core model

[FseofAll,multiFseof]= FSEOFTargets(model,TargetProducts);
TargetsScoreTable = Targets_HO_score(model,multiFseof);
