%Uses fastFVA to test the effect of amplification/deletion of any target on
%the flux of products and biomass
function [Result]=testresultsFVA(model,minBM,testRxns,met1,met2,var)

%%%% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%testRxns: the potential amp/ko target to be tested
%met1: exchange reaction of product A
%met2: exchange reaction of product B
%var: nature of intervention - amp/ko
%Result: array of intervention, min and max mutant flux of products,mutant
%biomass flux, wild-type product fluxes

[InterventionMin,fluxMinMet1,fluxMinMet2,mutantBiomassMin,InterventionMax,fluxMaxMet1,fluxMaxMet2,mutantBiomassMax]= deal({});
[fbaMet1,fbaMet2,fvaMaxMet1,fvaMaxMet2,biomass,fvaMinMet1,fvaMinMet2]= deal(zeros(length(testRxns),1));
rxnsList = horzcat(met1,met2);

%maximum flux for secreted metabolites in wild type
[minFlux,maxFlux] = fluxVariability(model,100,'max',rxnsList);
minMetWT = {minFlux(1),minFlux(2)};
maxMetWT = {maxFlux(1),maxFlux(2)};

%set minimum biomass for mutant %default=25%
fbaWT = optimizeCbModel(model);
model = changeRxnBounds(model,model.rxns(model.c==1),minBM*fbaWT.f,'l');


for i=1:length(testRxns)
    %for amplification targets set value to maxflux
    if strcmp(var,'amp')
        index = find(strcmp(model.rxns,testRxns(i)));
        modeltest = changeObjective(model,model.rxns(index),1);
        soltest = optimizeCbModel(modeltest);
        modelMut = changeRxnBounds(model,testRxns(i),soltest.x(index)-1e-3,'l');
        modelMut = changeRxnBounds(modelMut,testRxns(i),soltest.x(index)+1e-3,'u');
        solMut = optimizeCbModel(modelMut);
        if solMut.stat==1
            fbaMet1(i) = solMut.x(strcmp(model.rxns,met1));
            fbaMet2(i) = solMut.x(strcmp(model.rxns,met2));
            biomass(i) = solMut.f;
            rxnsList= horzcat(met1,met2);
            [minFluxMut,maxFluxMut] = fluxVariability(modelMut,100,'max',rxnsList);
            fvaMinMet1(i) = minFluxMut(1); fvaMinMet2(i) = minFluxMut(2); 
            fvaMaxMet1(i) = maxFluxMut(1); fvaMaxMet2(i) = maxFluxMut(2); 
            clearAllMemoizedCaches
        end
    %for knockout targets set value to zero
    elseif strcmp(var,'ko')
        modelMut = changeRxnBounds(model,testRxns(i),0,'b');
        solMut = optimizeCbModel(modelMut);
        if solMut.stat==1
            fbaMet1(i) = solMut.x(strcmp(model.rxns,met1));
            fbaMet2(i) = solMut.x(strcmp(model.rxns,met2));
            biomass(i) = solMut.f;
            rxnsList = horzcat(met1,met2);
            [minFluxMut,maxFluxMut] = fluxVariability(modelMut,100,'max',rxnsList);            
            fvaMinMet1(i) = minFluxMut(1); fvaMinMet2(i) = minFluxMut(2);
            fvaMaxMet1(i) = maxFluxMut(1); fvaMaxMet2(i) = maxFluxMut(2); 
            clearAllMemoizedCaches
        end
    end
end

%choose targets with non-zero flux improvement
temp=0;
for i=1:length(testRxns)
    if fvaMinMet1(i)>0.01 && fvaMinMet2(i)> 0.01 && fvaMinMet1(i)>1.05*minFlux(1) && fvaMinMet2(i)>1.05*minFlux(2) && biomass(i)>0
        temp=temp+1;
        InterventionMin(temp) = testRxns(i);
        fluxMinMet1(temp)=num2cell(fvaMinMet1(i));
        fluxMinMet2(temp)=num2cell(fvaMinMet2(i));
        mutantBiomassMin(temp) = num2cell(biomass(i));  
    end
end
temp=0;
for i=1:length(testRxns)
    if fvaMaxMet1(i)>0.01 && fvaMaxMet2(i)> 0.01 && fvaMaxMet1(i)>1.05*maxFlux(1) && fvaMaxMet2(i)>1.05*maxFlux(2) && biomass(i)>0
        temp=temp+1;
        InterventionMax(temp) = testRxns(i);
        fluxMaxMet1(temp)=num2cell(fvaMaxMet1(i));
        fluxMaxMet2(temp)=num2cell(fvaMaxMet2(i));
        mutantBiomassMax(temp) = num2cell(biomass(i));  
    end
end
Result={InterventionMin,fluxMinMet1,fluxMinMet2,mutantBiomassMin,minMetWT,InterventionMax,fluxMaxMet1,fluxMaxMet2,mutantBiomassMax,maxMetWT};
end
