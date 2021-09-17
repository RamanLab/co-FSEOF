%FSEOF %to perform FSEOF for a single product
function [FluxMatrix, IncFlux, DecFlux] = FSEOF (model,TargetRxn,eliListInc,eliListDec)
%%% input and output parameters
% %output
% FluxMatrix: the entire matrix of flux change when product flux is increased from wild type to theoretical maximum
% IncFlux: FluxMatrix of reactions that undergo increase in flux : potential amplification targets
% DecFlux: FluxMatrix of reactions that undergo decrease in flux : potential deletion targets
% %input
% model: the GSMM with appropriate medium bounds applied
% TargetRxn: The reaction for the target product for which FSEOF has to be performed
% eliListInc: the list of reactions to be excluded when finidng amplification targets
% eliListDec: the list of reactions to be excluded when finidng deletion targets

%Maximum biomass
FBAsolnWT = optimizeCbModel(model);
Vbiomass = FBAsolnWT.f;

%% Initial product and maximum product
Viniprdt = FBAsolnWT.x(strcmp(model.rxns,TargetRxn));  

model_target = changeObjective(model, TargetRxn);
FBAsoln_target = optimizeCbModel(model_target);
    if FBAsoln_target.stat == 1
        if FBAsoln_target.f > 1e-5
            Vmaxprdt = FBAsoln_target.f;
        else
            Vmaxprdt = 0;
            fprintf('metabolite cannot be produced, only gets consumed in');
            disp(TargetRxn);
        end
    else
        disp('Infeasible solution')
        Vmaxprdt = 0;
    end

 
%% Setup of flux matrix
steps = randperm(10);
FluxMatrix = zeros(length(model.rxns),length(steps)+1);
FluxMatrix(:,1) = FBAsolnWT.x;
IncFlux = [];
DecFlux = [];
 
    %% Enforcing flux
if Vmaxprdt >=1e-5
    for j =1:length(steps)
        Venfprdt = (steps(j)/10)*Vmaxprdt;
        model_mut = changeRxnBounds(model, TargetRxn, Venfprdt-1e-3, 'l');
        model_mut = changeRxnBounds(model_mut, TargetRxn, Venfprdt+1e-3, 'u');
        FBAsoln = optimizeCbModel(model_mut,'max','one'); 
        if FBAsoln.stat == 1
            FluxMatrix(:,steps(j)+1) = FBAsoln.x;
        else
            FluxMatrix(:,steps(j)+1) = NaN;
            fprintf('%02d: infeasible solution at Venfprdt = %f \n\n',steps(j),Venfprdt);
        end
    end
    
    %% identification of amplification and knockout targets
    for RxnIndex = 1:length(model.rxns)
        if abs(FluxMatrix (RxnIndex, 1)) < abs(FluxMatrix (RxnIndex, 3)) && abs(FluxMatrix (RxnIndex, 4)) < abs(FluxMatrix (RxnIndex, 6)) && abs(FluxMatrix (RxnIndex, 7)) < abs(FluxMatrix (RxnIndex, 9)) && ~ismember(model.rxns(RxnIndex),eliListInc)
            NewRowI = FluxMatrix(RxnIndex,:);
            Score = 1/((Vmaxprdt-Viniprdt)/abs(abs(FluxMatrix (RxnIndex, 9)) - abs(FluxMatrix (RxnIndex, 2))));
            NewRowI = [Score, NewRowI , RxnIndex, model.rxns(RxnIndex)];
            IncFlux = [IncFlux; NewRowI];

        elseif abs(FluxMatrix (RxnIndex, 1)) > abs(FluxMatrix (RxnIndex, 3)) && abs(FluxMatrix (RxnIndex, 4)) > abs(FluxMatrix (RxnIndex, 6)) && abs(FluxMatrix (RxnIndex, 7)) > abs(FluxMatrix (RxnIndex, 9)) && ~ismember(model.rxns(RxnIndex),eliListDec)
            NewRowII = FluxMatrix(RxnIndex,:);
            Score = 1/((Vmaxprdt-Viniprdt)/abs(abs(FluxMatrix (RxnIndex, 9)) - abs(FluxMatrix (RxnIndex, 2))));
            NewRowII = [Score, NewRowII , RxnIndex, model.rxns(RxnIndex)];
            DecFlux = [DecFlux ; NewRowII];

        end
    end

    if ~isempty(IncFlux) 
     IncFlux = sortrows(IncFlux, 1,'descend');
    end
    if ~isempty(DecFlux)
     DecFlux = sortrows(DecFlux, 1,'descend');
    end
     %Score in first column ; rows arranged according to inc order of scores
     %Flux ranging across steps 0(biomass obj) to 10(product obj) in columns 2-12
     %RxnNo in column 13
end
end

