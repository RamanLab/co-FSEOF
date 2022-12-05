%Uses fastFVA to test the effect of multiple amplification/deletion of any target on
%the flux of specific products and biomass
function [Result]=testresultsFVATargets(model,ProductRxnsList,TargetRxns,var)

%%%% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%ProductRxnsList: list of products for which the analysis is done
%TargetRxns: the potential amp/ko target to be tested
%var: nature of intervention - amp/ko
%Result: array of intervention, min and max mutant flux of products,mutant
%biomass flux, wild-type product fluxes


[InterventionMin,fluxMinMet1,fluxMinMet2,mutantBiomassMin,InterventionMax,fluxMaxMet1,fluxMaxMet2,mutantBiomassMax]= deal({});
[fvaMaxMet1,fvaMaxMet2,biomass,fvaMinMet1,fvaMinMet2] = deal([]);

%maximum flux for secreted metabolites in wild type
[minFlux,maxFlux] = fastFVA(model,100,'max','ibm_cplex',ProductRxnsList);
minMetWT = {minFlux(1),minFlux(2)};
maxMetWT = {maxFlux(1),maxFlux(2)};

%set minimum biomass for mutant %default=25%
fbaWT = optimizeCbModel(model);
model = changeRxnBounds(model,model.rxns(model.c==1),0.25*fbaWT.f,'l');

%% amplification of target reactions
if strcmp(var,'amp')
    inv = 0;
    for intrvSize = 1:3
        testRxns = nchoosek(TargetRxns,intrvSize);

        for i=1:length(testRxns)
            %for amplification targets set value to maxflux
            for j = 1: size(testRxns,2)
                index(j) = find(strcmp(model.rxns,testRxns(i,j)));
                modeltest = changeObjective(model,model.rxns(index(j)),1);
                soltest = optimizeCbModel(modeltest);
                soltestFlux(j) = soltest.x(index(j));
                modelMut = changeRxnBounds(model,testRxns(i,:),soltestFlux(j)-1e-3,'l'); %%%%%%%%%%%%%
                modelMut = changeRxnBounds(modelMut,testRxns(i,:),soltestFlux(j)+1e-3,'u');%%%%%%%%%%%%%%%
                solMut = optimizeCbModel(modelMut);
                if solMut.stat==1
                    inv = inv+1;
                    biomass(inv) = solMut.f;
                    [minFluxMut,maxFluxMut] = fastFVA(modelMut,100,'max','ibm_cplex',ProductRxnsList);
                    fvaMinMet1(inv) = minFluxMut(1); fvaMinMet2(inv) = minFluxMut(2);
                    fvaMaxMet1(inv) = maxFluxMut(1); fvaMaxMet2(inv) = maxFluxMut(2);
                    Intervention{inv} = testRxns(i,:);
                    clearAllMemoizedCaches
                end
            end
        end
    end

%% deletion of target reactions
elseif strcmp(var,'ko')
    inv = 0;
    for intrvSize = 1:3
        testRxns = nchoosek(TargetRxns,intrvSize);
        % for knockout targets set value to zero
     for i=1:length(testRxns)
        modelMut = changeRxnBounds(model,testRxns(i,:),0,'b');
        solMut = optimizeCbModel(modelMut);
        if solMut.stat==1
            inv=inv+1;
            biomass(inv) = solMut.f;
            [minFluxMut,maxFluxMut] = fastFVA(modelMut,100,'max','ibm_cplex',ProductRxnsList);            
            fvaMinMet1(inv) = minFluxMut(1); fvaMinMet2(inv) = minFluxMut(2);
            fvaMaxMet1(inv) = maxFluxMut(1); fvaMaxMet2(inv) = maxFluxMut(2);
            Intervention{inv} = testRxns(i,:);
            clearAllMemoizedCaches
        end
    end
    end
end

%% choose targets with non-zero flux improvement
temp=0;
for i=1:inv
    if fvaMinMet1(i)>0.01 && fvaMinMet2(i)> 0.01 && fvaMinMet1(i)>1.05*minFlux(1) && fvaMinMet2(i)>1.05*minFlux(2) && biomass(i)>0
        temp=temp+1;
        InterventionMin(temp) = Intervention(i);
        fluxMinMet1(temp)=num2cell(fvaMinMet1(i));
        fluxMinMet2(temp)=num2cell(fvaMinMet2(i));
        mutantBiomassMin(temp) = num2cell(biomass(i));  
    end
end
temp=0;
for i=1:inv
    if fvaMaxMet1(i)>0.01 && fvaMaxMet2(i)> 0.01 && fvaMaxMet1(i)>1.05*maxFlux(1) && fvaMaxMet2(i)>1.05*maxFlux(2) && biomass(i)>0
        temp=temp+1;
        InterventionMax(temp) = Intervention(i);
        fluxMaxMet1(temp)=num2cell(fvaMaxMet1(i));
        fluxMaxMet2(temp)=num2cell(fvaMaxMet2(i));
        mutantBiomassMax(temp) = num2cell(biomass(i));  
    end
end
Result= vertcat(InterventionMax,fluxMaxMet1,fluxMaxMet2,mutantBiomassMax,InvTypeMax(1:length(InterventionMax)))';

end
