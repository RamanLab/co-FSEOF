% to test and score higher-order interventions
function [TargetsScoreTable] = Target_HO_score(model,multiFseof)

%model: the GSMM with appropriate medium bounds applied
%multiFseof: A matrix of reaction for product A,reaction for product B, amplification targets, deletion targets
%TargetScoreTable: A table with all interventions and scores

%to test higher-order interventions
[Amp, KO, AmpKO_2, AmpKO_3] = deal({});
Amp = testresultsFVATargets(model,multiFseof{:,1},multiFseof{:,2},{'amp'}); %multiple amps (A)
KO = testresultsFVATargets(model,multiFseof{:,1},multiFseof{:,3},{'ko'}); %multiple KOs (K)
AmpKO_2 = testresultsFVAMixedTargets(model,multiFseof{:,1},multiFseof{:,2},multiFseof{:,3},2); %Amp+KO (AK)
AmpKO_3 = testresultsFVAMixedTargets(model,multiFseof{:,1},multiFseof{:,2},multiFseof{:,3},3); %Amp+Amp+KO, Amp+KO+KO (AAK,AKK)

[~,maxFlux] = fluxVariability(model,100,'max',multiFseof{1,1});
Targets = vertcat(Amp,KO,AmpKO_2,AmpKO_3);
modelSol = optimizeCbModel(model);

for i=1: size(Targets)
    Targets{i,6} = (Targets{i,2}-maxFlux(1))/(modelSol.f-Targets{i,4});
    Targets{i,7} = (Targets{i,3}-maxFlux(2))/(modelSol.f-Targets{i,4});
    Targets{i,8} = Targets{i,6}+Targets{i,7};
end

TargetsScoreTable = sortrows(Targets,8,{'descend'}); %sorted in descending order of Score A+B

header = {'Intervention set','Mutant flux of prd A','Mutant flux of prd B','Mutant biomass flux','Type of intervention','Score A','Score B','Score A+B'};
TargetsScoreTable = [header;TargetsScoreTable];

end
