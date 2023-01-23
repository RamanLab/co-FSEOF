# Co-production
Algorithm to study co-production of multiple metabolites

Analysis 1:
Choose a metabolic model and study all pairs of metabolites that can be co-produced/co-optimized with single intervention (amplification/deletion) in the model.
The intervention strategies along with appropriate scores can be found in the output multiFseofScoreTable

Analysis 2:
Choose a specific set of metabolites and identify higher order interventions of size upto three.
The intervention strategies along with appropriate scores can be found in the output TargetsScoreTable

Coproduction_masterScript.m is the main code for studying co-production. It has files for both all metabolites and higher-order interventions for single set of metabolites.
Organisms studied in master script: Ecoli (iML15151.xml) and yeast (iMM904.xml)

Scripts common to both analysis: FSEOF.m

Scripts for co-production of all pairs of metabolites in a model: 
FSEOFall.m; coFSEOF.m; testresultsFVA.m

Scripts for co-production of specific set of metabolites in a model: 
coFSEOFTargets.m; testresultsFVATargets.m; testresultsFVAMixedTargets.m

All the analysis were performed on MATLAB R2018a and COBRAToolbox v.3.0 and IBM Cplex 12.8 solver.
