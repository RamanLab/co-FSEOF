%to perform FSEOF for given set of metabolites in a model
function[TargetsScoreTable]= coFSEOFTargets(model,minBM,solver,TargetProducts)

%%% input and output parameters and other data
%model: the GSMM with appropriate medium bounds applied
%TargetScoreTable: A table with all interventions and scores
%TargetProduct: list of products for which intervention strategies have to be found
%FSEOFAllTargets: A matrix of reactions corr to each product of interest, amplification targets, deletion targets
%coFseofTargets: A matrix of reactions corr to sets of products of interest, amplification targets, deletion targets
%% set up elimination lists
    %finding exchange and transport rxns
    excRxns = model.rxns(findExcRxns(model)==1);
    transRxns = union(excRxns,findTransRxns(model));

    %elimination list for amplification targets
    eliListInc = union(transRxns,'ATPM');
    eliListInc = union(eliListInc, model.rxns(model.c==1));

    %elimination list for knockout targets
    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model);
    EssRxns = model.rxns(grRatio < 1e-3);
    eliListDec = union(eliListInc,EssRxns);
 %% running FSEOF for all exchange metabolites
    for i = 1:length(TargetProducts)
        [FluxMatrix, IncFlux, DecFlux] = FSEOF (model,TargetProducts(i),eliListInc, eliListDec);
        if ~isempty(IncFlux) 
            IncAll{i,1} = IncFlux(:,4); 
        else
            IncAll{i,1} = [];
        end
        if ~isempty(DecFlux)
            DecAll{i,1} = DecFlux(:,4); 
        else 
            DecAll{i,1} = [];
        end
    end
    Fseof_temp = horzcat(TargetProducts,IncAll,DecAll);

    %removing empty entries
    temp=0;
    for j=1:length(TargetProducts)    
        if ~isempty(Fseof_temp{j,2}) || ~isempty(Fseof_temp{j,3})
            temp=temp+1;
            FseofAllTargets(temp,:) = Fseof_temp(j,:);
        end
    end

%%     
coFseofTargets{:,1} = (horzcat(FseofAllTargets(1:end,1)));
coFseofTargets{:,2} = unique(vertcat(FseofAllTargets{1:end,2}));
coFseofTargets{:,3} = unique(vertcat(FseofAllTargets{1:end,3}));

%%scoring the intervention targets
if length(coFseofTargets{1,1})>1
    %to test higher-order interventions
    [Amp, KO, AmpKO_2, AmpKO_3] = deal({});
    Amp = testresultsFVATargets(model,minBM,solver,coFseofTargets{:,1},coFseofTargets{:,2},{'amp'}); %multiple amps (A)
    KO = testresultsFVATargets(model,minBM,solver,coFseofTargets{:,1},coFseofTargets{:,3},{'ko'}); %multiple KOs (K)
    AmpKO_2 = testresultsFVAMixedTargets(model,minBM,solver,coFseofTargets{:,1},coFseofTargets{:,2},coFseofTargets{:,3},2); %Amp+KO (AK)
    AmpKO_3 = testresultsFVAMixedTargets(model,minBM,solver,coFseofTargets{:,1},coFseofTargets{:,2},coFseofTargets{:,3},3); %Amp+Amp+KO, Amp+KO+KO (AAK,AKK)

    [~,maxFlux] = fluxVariability(model,100,'max',coFseofTargets{1,1});
    Targets = vertcat(Amp,KO,AmpKO_2,AmpKO_3);
    modelSol = optimizeCbModel(model);

    for i=1: size(Targets)
        Targets{i,6} = (Targets{i,2}-maxFlux(1))/(modelSol.f-Targets{i,4});
        Targets{i,7} = (Targets{i,3}-maxFlux(2))/(modelSol.f-Targets{i,4});
        Targets{i,8} = Targets{i,6}+Targets{i,7};
    end

    TargetsSorted = sortrows(Targets,8,{'descend'}); %sorted in descending order of Score A+B
    TargetsScoreTable(:,4:10) = TargetsSorted(:,2:8);
    for i= 1:length(TargetsSorted)
        if length(TargetsSorted{i,1}) == 1
            TargetsScoreTable{i,1} = TargetsSorted{i,1};
            TargetsScoreTable{i,2} = {};
            TargetsScoreTable{i,3} = {};
        end
        if length(TargetsSorted{i,1}) == 2
            TargetsScoreTable{i,1} = TargetsSorted{i,1}{1,1};
            TargetsScoreTable{i,2} = TargetsSorted{i,1}{1,2};
            TargetsScoreTable{i,3} = {};
        end
        if length(TargetsSorted{i,1}) == 3
            TargetsScoreTable{i,1} = TargetsSorted{i,1}{1,1};
            TargetsScoreTable{i,2} = TargetsSorted{i,1}{1,2};
            TargetsScoreTable{i,3} = TargetsSorted{i,1}{1,3};
        end
    end
    
    header = {'Intervention1','Intervention2','Intervention3','Mutant flux of prd A','Mutant flux of prd B','Mutant biomass flux','Type of intervention','Score A','Score B','Score A+B'};
    TargetsScoreTable = [header;TargetsScoreTable];
else
    disp ('The given combination of products cannot be co-produced');
    TargetsScoreTable = [];
end

end


