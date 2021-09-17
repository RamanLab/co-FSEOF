% to collate the table and score the successful intervention strategies
function[multiFseofScoreTable] = multiFseof_Score(model,multiFseof,Amp,KO)

%%%% input and output parameters
%model: the GSMM with appropriate medium bounds applied
%multiFseof: A matrix of reaction for product A,reaction for product B, amplification targets, deletion targets
%Amp: data matrix with intervention, min and max flux improvement in both products,mutant biomass and wild-type fluxes of products for Amplification targets
%KO: data matrix with intervention, min and max flux improvement in both products,mutant biomass and wild-type fluxes of products for deletion targets
%multiFseofScoreTable: output table with all successful intervention strategies along with their scores
    
multiFseofScore = multiFseof;
    %rearranging the table
    for i=1:size(multiFseof,1)
       multiFseofScore{i,5}=vertcat(Amp{i,1:4})';
       multiFseofScore{i,6}=vertcat(Amp{i,6:9})';
       multiFseofScore{i,7}=vertcat(KO{i,1:4})';
       multiFseofScore{i,8}=vertcat(KO{i,6:9})';
       multiFseofScore{i,9}=Amp{i,5};
       multiFseofScore{i,10}=Amp{i,10};
    end
    biomassWT = optimizeCbModel(model);
    %Scoring flux improvement
    for k = 5:8
        for i=1:size(multiFseofScore)
            %Min Flux improvement through amplification
            if ~isempty(multiFseofScore{i,k})
                for j=1:size(multiFseofScore{i,k},1)      
                        multiFseofScore{i,k}{j,5} = (multiFseofScore{i,k}{j,2}-multiFseofScore{i,9}{1,1})/(biomassWT.f-multiFseofScore{i,k}{j,4}); % Score A     
                        multiFseofScore{i,k}{j,6} = (multiFseofScore{i,k}{j,3}-multiFseofScore{i,9}{1,2})/(biomassWT.f-multiFseofScore{i,k}{j,4}); % Score B           
                        multiFseofScore{i,k}{j,7} = (multiFseofScore{i,k}{j,5}+ multiFseofScore{i,k}{j,6}); % Score A+B
                        if j == size(multiFseofScore{i,k},1) 
                            header = {'Intervention','Mutant flux of ProductA','Mutant flux of ProductB','Mutant biomass flux','Score A','Score B','Score A+B'};
                            multiFseofScore{i,k} = [header;multiFseofScore{i,k}];
                        end
                end
            end
        end   
    end

    %collate non-empty entries in multiFseof_Score to construct table
    z = multiFseofScore(:,6);
    maxcols = max(cellfun('size', z, 2));  %get the number of columns of the widest array
    padded = cellfun(@(m) [m, zeros(size(m, 1), maxcols - size(m, 2))], z, 'UniformOutput', false);  %pad each array
    z_Reference  = vertcat(padded{:});
    temp=0;
    for i=1: length(multiFseofScore)
    if ~cellfun(@isempty,multiFseofScore(i,5)) || ~cellfun(@isempty,multiFseofScore(i,6))|| ~cellfun(@isempty,multiFseofScore(i,7)) || ~cellfun(@isempty,multiFseofScore(i,8))
    temp=temp+1;
    multiFseofScoreTable(temp,:) = multiFseofScore(i,:);
    end
    end
    %insert header
    header = {'Product A','Product B','Potential amplification targets','Potential deletion targets', 'Improvement in Min flux through amplification','Improvement in Max flux through amplification','Improvement in Min flux through deletion','Improvement in Max flux through deletion','Wild type min flux of products A and B','Wild type max flux of products A and B'};
    multiFseofScoreTable = [header;multiFseofScoreTable];
end
      
            