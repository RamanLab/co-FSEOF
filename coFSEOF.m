%to find co-optimization targets for all pairs of metabolites in a model
function[coFseofScoreTable] = coFSEOF(model,minBM,solver)

%%%% input and output parameters and other important data
%model: the GSMM with appropriate medium bounds applied
%minBM: minimum biomass of mutant - given in percentage of wild-type biomass
%solver: solver name
%coFseofScoreTable: output table with all successful intervention strategies along with their scores
%FSEOFAll: A matrix of reaction, amplification targets, deletion targets
%coFseof: A matrix of reaction for product A,reaction for product B, amplification targets, deletion targets
%Amp: data matrix with intervention, min and max flux improvement in both products,mutant biomass and wild-type fluxes of products for Amplification targets
%KO: data matrix with intervention, min and max flux improvement in both products,mutant biomass and wild-type fluxes of products for deletion targets

coFseof = {};
temp =0;    
%running FSEOF for all exchange metabolites
FseofAll= FSEOFall(model);

    %finding common targets
    for i=1:length(FseofAll)
        for j=i+1:length(FseofAll)
            temp=temp+1;
            coFseof{temp,1}=FseofAll(i,1);
            coFseof{temp,2}=FseofAll(j,1);
            if ~isempty(FseofAll{i,2}) && ~isempty(FseofAll{j,2})
                a=FseofAll{i,2}; b=FseofAll{j,2};
                coFseof{temp,3}=intersect(a,b);
            end
            if ~isempty(FseofAll{i,3}) && ~isempty(FseofAll{j,3})
                coFseof{temp,4}=intersect(FseofAll{i,3},FseofAll{j,3});
            end
        end
    end
    %testing the targets using FVA
    [Amp,KO] = deal(cell(length(coFseof),10));
    for i=1:length(coFseof)
        Amp(i,:) = testresultsFVA(model,minBM,solver,coFseof{i,3},coFseof{i,1},coFseof{i,2},{'amp'});
        KO(i,:) = testresultsFVA(model,minBM,solver,coFseof{i,4},coFseof{i,1},coFseof{i,2},{'ko'});
    end
 %% scoring the different targets obtained   
    coFseofScore = coFseof;
    %rearranging the table
    for i=1:size(coFseof,1)
       coFseofScore{i,5}=vertcat(Amp{i,1:4})';
       coFseofScore{i,6}=vertcat(Amp{i,6:9})';
       coFseofScore{i,7}=vertcat(KO{i,1:4})';
       coFseofScore{i,8}=vertcat(KO{i,6:9})';
       coFseofScore{i,9}=Amp{i,5};
       coFseofScore{i,10}=Amp{i,10};
    end
    biomassWT = optimizeCbModel(model);
    %Scoring flux improvement
    for k = 5:8
        for i=1:size(coFseofScore)
            %Min Flux improvement through amplification
            if ~isempty(coFseofScore{i,k})
                for j=1:size(coFseofScore{i,k},1)      
                        coFseofScore{i,k}{j,5} = (coFseofScore{i,k}{j,2}-coFseofScore{i,9}{1,1})/(biomassWT.f-coFseofScore{i,k}{j,4}); % Score A     
                        coFseofScore{i,k}{j,6} = (coFseofScore{i,k}{j,3}-coFseofScore{i,9}{1,2})/(biomassWT.f-coFseofScore{i,k}{j,4}); % Score B           
                        coFseofScore{i,k}{j,7} = (coFseofScore{i,k}{j,5}+ coFseofScore{i,k}{j,6}); % Score A+B
                        if j == size(coFseofScore{i,k},1) 
                            header = {'Intervention','Mutant flux of ProductA','Mutant flux of ProductB','Mutant biomass flux','Score A','Score B','Score A+B'};
                            coFseofScore{i,k} = [header;coFseofScore{i,k}];
                        end
                end
            end
        end   
    end
%%
    %collate non-empty entries in coFseof_Score to construct table
    z = coFseofScore(:,6);
    maxcols = max(cellfun('size', z, 2));  %get the number of columns of the widest array
    padded = cellfun(@(m) [m, zeros(size(m, 1), maxcols - size(m, 2))], z, 'UniformOutput', false);  %pad each array
    z_Reference  = vertcat(padded{:});
    temp=0;
    for i=1: length(coFseofScore)
    if ~cellfun(@isempty,coFseofScore(i,5)) || ~cellfun(@isempty,coFseofScore(i,6))|| ~cellfun(@isempty,coFseofScore(i,7)) || ~cellfun(@isempty,coFseofScore(i,8))
    temp=temp+1;
    coFseofScoreTable(temp,:) = coFseofScore(i,:);
    end
    end
    %insert header
    header = {'Product A','Product B','Potential amplification targets','Potential deletion targets', 'Improvement in Min flux through amplification','Improvement in Max flux through amplification','Improvement in Min flux through deletion','Improvement in Max flux through deletion','Wild type min flux of products A and B','Wild type max flux of products A and B'};
    coFseofScoreTable = [header;coFseofScoreTable];
end
 
