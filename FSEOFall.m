%to perform FSEOF for all secretory or exchange metabolites in a model
function[FseofAll]= FSEOFall(model)
%%% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%FSEOFAll: A matrix of reaction, amplification targets, deletion targets
    %% set up elimination lists
    %finding exchange and transport rxns
    excRxns = model.rxns(findExcRxns(model)==1); % to remove exchange reactions
    transRxns = union(excRxns,findTransRxns(model)); % to remove transport reactions

    %elimination list for amplification targets
    eliListInc = union(transRxns,'ATPM');
    eliListInc = union(eliListInc, model.rxns(model.c==1));

    %elimination list for knockout targets
    fastSL(model,0.01,1,model.rxns(findExcRxns(model)==1)); 
    load(horzcat(model.description,'_Rxn_lethals.mat'));
    eliListDec = union(eliListInc,Jsl); % to remove single lethals

    %% running FSEOF for all exchange metabolites
    for i = 1:length(excRxns)
        [FluxMatrix, IncFlux, DecFlux] = FSEOF (model,excRxns(i),eliListInc, eliListDec);
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
    Fseof_temp = horzcat(excRxns,IncAll,DecAll);

    %removing empty entries
    temp=0;
    for j=1:length(excRxns)    
        if ~isempty(Fseof_temp{j,2}) || ~isempty(Fseof_temp{j,3})
            temp=temp+1;
            FseofAll(temp,:) = Fseof_temp(j,:);
        end
    end


end
