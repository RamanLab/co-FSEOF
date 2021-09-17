%to find co-optimization targets for all pairs of metabolites in a model
function[FseofAll,multiFseof,Amp,KO] = multiFSEOF(model)

%%%% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%FSEOFAll: A matrix of reaction, amplification targets, deletion targets
%multiFseof: A matrix of reaction for product A,reaction for product B, amplification targets, deletion targets
%Amp: data matrix with intervention, min and max flux improvement in both products,mutant biomass and wild-type fluxes of products for Amplification targets
%KO: data matrix with intervention, min and max flux improvement in both products,mutant biomass and wild-type fluxes of products for deletion targets

multiFseof = {};
temp =0;    
%running FSEOF for all exchange metabolites
FseofAll= FSEOFall(model);

    %finding common targets
    for i=1:length(FseofAll)
        for j=i+1:length(FseofAll)
            temp=temp+1;
            multiFseof{temp,1}=FseofAll(i,1);
            multiFseof{temp,2}=FseofAll(j,1);
            if ~isempty(FseofAll{i,2}) && ~isempty(FseofAll{j,2})
                a=FseofAll{i,2}; b=FseofAll{j,2};
                multiFseof{temp,3}=intersect(a,b);
            end
            if ~isempty(FseofAll{i,3}) && ~isempty(FseofAll{j,3})
                multiFseof{temp,4}=intersect(FseofAll{i,3},FseofAll{j,3});
            end
        end
    end
    %testing the targets using FVA
    [Amp,KO] = deal(cell(length(multiFseof),10));
    for i=1:length(multiFseof)
        Amp(i,:) = testresultsFVA(model,multiFseof{i,3},multiFseof{i,1},multiFseof{i,2},{'amp'});
        KO(i,:) = testresultsFVA(model,multiFseof{i,4},multiFseof{i,1},multiFseof{i,2},{'ko'});
    end
    
end
 