%to perform FSEOF for given set of metabolites in a model
function[FseofAll,multiFseof]= FSEOFTargets(model,TargetProducts)

%%% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%FSEOFAll: A matrix of reactions corr to product of interest, amplification targets, deletion targets
%FSEOFAll: A matrix of reactions corr to sets of products of interest, amplification targets, deletion targets
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
            FseofAll(temp,:) = Fseof_temp(j,:);
        end
    end

%%     
multiFseof{:,1} = (horzcat(FseofAll(1:end,1)));
multiFseof{:,2} = unique(vertcat(FseofAll{1:end,2}));
multiFseof{:,3} = unique(vertcat(FseofAll{1:end,3}));

end


