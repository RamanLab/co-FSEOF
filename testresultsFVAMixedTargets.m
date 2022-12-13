%Uses fastFVA to test the effect of mixture of amplification/deletion of any target on
%the flux of specific products and biomass

function [Result]=testresultsFVAMixedTargets(model,minBM,solver,ProductRxnsList,AmpTargets,KoTargets,intrvSize)

%%%% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%ProductRxnsList: list of products for which the analysis is done
%AmpTargets: the potential amp target to be tested
%KOTargets: the potential KO target to be tested
%intrvSize: size of intervention
%Result: array of intervention, min and max mutant flux of products,mutant
%biomass flux, wild-type product fluxes


[InterventionMin,fluxMinMet1,fluxMinMet2,mutantBiomassMin,InterventionMax,fluxMaxMet1,fluxMaxMet2,mutantBiomassMax]= deal({});
[Intervention,invType,invTypeMin,invTypeMax]= deal({});
[fvaMinMet1,fvaMinMet2,fvaMaxMet1,fvaMaxMet2,biomass] = deal([]);

%maximum flux for secreted metabolites in wild type
if strcmp(solver,'ibm_cplex')
    [minFlux,maxFlux] = fastFVA(model,100,'max','ibm_cplex',ProductRxnsList);
else
    [minFlux,maxFlux] = fluxVariability(model,100,'max',ProductRxnsList);
end
minMetWT = {minFlux(1),minFlux(2)};
maxMetWT = {maxFlux(1),maxFlux(2)};

%set minimum biomass for mutant %default=25%
fbaWT = optimizeCbModel(model);
model = changeRxnBounds(model,model.rxns(model.c==1),minBM*fbaWT.f,'l');
%% Amp + KO
if intrvSize == 2
    inv=0;
    for i=1:length(AmpTargets)
        for j= 1:length(KoTargets)
            indexAmp1 = find(strcmp(model.rxns,AmpTargets(i)));
            modeltest = changeObjective(model,model.rxns(indexAmp1),1);
             soltest = optimizeCbModel(modeltest);
             soltestFlux = soltest.x(indexAmp1);
             modelMut = changeRxnBounds(model,AmpTargets(i),0.5*soltestFlux-1e-3,'l'); 
             modelMut = changeRxnBounds(modelMut,AmpTargets(i),0.5*soltestFlux+1e-3,'u'); 
             modelMut = changeRxnBounds(modelMut,KoTargets(j),0,'b');
             solMut = optimizeCbModel(modelMut);
           
             if solMut.stat==1
                inv=inv+1;
                biomass(inv) = solMut.f;
                if strcmp(solver,'ibm_cplex')
                    [minFluxMut,maxFluxMut] = fastFVA(modelMut,100,'max','ibm_cplex',ProductRxnsList);
                else
                    [minFluxMut,maxFluxMut] = fluxVariability(modelMut,100,'max',ProductRxnsList);
                end         
                fvaMinMet1(inv) = minFluxMut(1); fvaMinMet2(inv) = minFluxMut(2);
                fvaMaxMet1(inv) = maxFluxMut(1); fvaMaxMet2(inv) = maxFluxMut(2);
                Intervention{inv} = [AmpTargets(i),KoTargets(j)];
                invType(inv) = {'AK'};
                clearAllMemoizedCaches
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
            invTypeMin(temp) = invType(i);
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
            invTypeMax(temp) = invType(i);
        end
    end
    Result= vertcat(InterventionMax,fluxMaxMet1,fluxMaxMet2,mutantBiomassMax,invTypeMax)';
    

elseif intrvSize == 3
    %% AMP + KO + KO
 inv=0;
 for i=1:length(AmpTargets)
     for j= 1:length(KoTargets)-1  
         indexAmp1 = find(strcmp(model.rxns,AmpTargets(i)));
            modeltest = changeObjective(model,model.rxns(indexAmp1),1);
             soltest = optimizeCbModel(modeltest);
             soltestFlux = soltest.x(indexAmp1);
             modelMut = changeRxnBounds(model,AmpTargets(i),0.5*soltestFlux -1e-3,'l'); 
             modelMut = changeRxnBounds(modelMut,AmpTargets(i),0.5*soltestFlux +1e-3,'u'); 
             modelMut = changeRxnBounds(modelMut,KoTargets(j:j+1),0,'b');      
             solMut = optimizeCbModel(modelMut);
           
             if solMut.stat==1
                inv=inv+1;
                biomass(inv) = solMut.f;
               if strcmp(solver,'ibm_cplex')
                    [minFluxMut,maxFluxMut] = fastFVA(modelMut,100,'max','ibm_cplex',ProductRxnsList);
                else
                    [minFluxMut,maxFluxMut] = fluxVariability(modelMut,100,'max',ProductRxnsList);
                end              
                fvaMaxMet1(inv) = maxFluxMut(1); fvaMaxMet2(inv) = maxFluxMut(2);
                Intervention{inv} = [AmpTargets(i),KoTargets(j),KoTargets(j+1)];
                invType(inv) = {'AKK'};
                clearAllMemoizedCaches
             end
     end
 end
    
    
    %% AMP + AMP + KO
    for i=1:length(AmpTargets)-1
     for j= 1:length(KoTargets)  
         indexAmp1 = find(strcmp(model.rxns,AmpTargets(i)));
         modeltest = changeObjective(model,model.rxns(indexAmp1),1);
         soltest = optimizeCbModel(modeltest);
         soltestFlux = soltest.x(indexAmp1);
         modelMut = changeRxnBounds(model,AmpTargets(i),0.5*soltestFlux-1e-3,'l'); 
         modelMut = changeRxnBounds(modelMut,AmpTargets(i),0.5*soltestFlux+1e-3,'u'); 
         
         indexAmp2 = find(strcmp(model.rxns,AmpTargets(i+1)));
         modeltest = changeObjective(model,model.rxns(indexAmp2),1);
         soltest = optimizeCbModel(modeltest);
         soltestFlux = soltest.x(indexAmp2);
         modelMut = changeRxnBounds(modelMut,AmpTargets(i+1),0.5*soltestFlux-1e-3,'l'); 
         modelMut = changeRxnBounds(modelMut,AmpTargets(i+1),0.5*soltestFlux+1e-3,'u');     
         modelMut = changeRxnBounds(modelMut,KoTargets(j),0,'b');
         solMut = optimizeCbModel(modelMut);
         
         if solMut.stat==1
                inv=inv+1;
                biomass(inv) = solMut.f;
                if strcmp(solver,'ibm_cplex')
                    [minFluxMut,maxFluxMut] = fastFVA(modelMut,100,'max','ibm_cplex',ProductRxnsList);
                else
                    [minFluxMut,maxFluxMut] = fluxVariability(modelMut,100,'max',ProductRxnsList);
                end        
                fvaMinMet1(inv) = minFluxMut(1); fvaMinMet2(inv) = minFluxMut(2);
                fvaMaxMet1(inv) = maxFluxMut(1); fvaMaxMet2(inv) = maxFluxMut(2);
                Intervention{inv} = [AmpTargets(i),AmpTargets(i+1),KoTargets(j)];
                invType(inv) = {'AAK'};
                clearAllMemoizedCaches
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
            invTypeMin(temp) = invType(i);
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
            invTypeMax(temp) = invType(i);
        end
    end
    Result= vertcat(InterventionMax,fluxMaxMet1,fluxMaxMet2,mutantBiomassMax,invTypeMax)';
end
end     
