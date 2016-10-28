function [model LL] = EM_tensorSTGMM(Data, model)
%%% EM procedure for training of Task-Parameterized Semi-Tied Gaussian Mixture Model
%
% If you found this work useful, please cite:
%
% @article{Tanwani16RAL,
%   author="Tanwani, A. K. and Calinon, S.",
%   title="Learning Robot Manipulation Tasks with Task-Parameterized Semi-Tied Hidden Semi-{M}arkov Model",
%   journal="{IEEE} Robotics and Automation Letters ({RA-L})",
%   year="2016",
%   month="January",
%   volume="1",
%   number="1",
%   pages="235--242",
% 	doi="10.1109/LRA.2016.2517825"
% }
%
% Copyright (c) 2016 Idiap Research Institute
% Written by Ajay Tanwani <http://www.ajaytanwani.com>
%%
%Parameters of the EM algorithm
nbMinSteps = 5; %Minimum number of iterations allowed
nbMaxSteps = 100; %Maximum number of iterations allowed
maxDiffLL = 1E-4; %Likelihood increase threshold to stop the algorithm
diagRegularizationFactor = 1E-5;
B_sf = 1E-1;
alpha = 0.99;

%%% overwrite parameters if defined in model structure
if isfield(model,'nbMinSteps'), nbMinSteps = model.nbMinSteps; disp(['nbMinSteps = ', num2str(nbMinSteps)]); end
if isfield(model,'nbMaxSteps'), nbMaxSteps = model.nbMaxSteps; disp(['nbMaxSteps = ', num2str(nbMaxSteps)]); end
if isfield(model,'maxDiffLL'), maxDiffLL = model.maxDiffLL; disp(['maxDiffLL = ', num2str(maxDiffLL)]); end
if isfield(model,'diagRegularizationFactor'), diagRegularizationFactor = model.diagRegularizationFactor; disp(['diagRegularizationFactor = ', num2str(diagRegularizationFactor)]); end

if isfield(model,'B_sf'), B_sf = model.B_sf; disp(['B_sf = ', num2str(B_sf)]);  end
if isfield(model,'alpha'), alpha = model.alpha; disp(['alpha = ', num2str(alpha)]); end

nbData = size(Data,3);

if ~isfield(model,'B')
    Itmp = eye(model.nbVar)*B_sf;
    model.B = repmat(Itmp, [1 1 model.nbFrames]);

    for m=1:model.nbFrames
        model.InitH(:,:,m) = pinv(model.B(:,:,m)) + eye(model.nbVar) * diagRegularizationFactor;
    end

    % model.Sigma = repmat(GenerateRandomSymmetricMatrix(model.nbVar),[1 1 model.nbFrames model.nbStates]);
    for i=1:model.nbStates
        for m=1:model.nbFrames
    %         model.InitSigmaDiag(:,:,m,i) = diag(diag(model.B(:,:,m)*model.Sigma(:,:,m,i)*model.B(:,:,m)'));  
            model.InitSigmaDiag(:,:,m,i) = eig(squeeze(model.Sigma(:,:,m,i)));  
        end
    end
end

for nbIter=1:nbMaxSteps    
    disp([num2str(nbIter) '. ']);
	
	%E-step
	[L, GAMMA] = computeGamma(Data, model); %See 'computeGamma' function below
	GAMMA2 = GAMMA ./ repmat(sum(GAMMA,2),1,nbData);
	model.Pix = GAMMA2;

    %M-step 
	for i=1:model.nbStates
		%Update Priors
		model.Priors(i) = sum(GAMMA(i,:)) / nbData;
		
        for m=1:model.nbFrames
            %Matricization/flattening of tensor
			DataMat(:,:) = Data(:,m,:);
            
            %Update Mu
            model.Mu(:,m,i) = DataMat * GAMMA2(i,:)';

            %Update Sigma
            DataTmp = DataMat - repmat(model.Mu(:,m,i),1,nbData);
            model.S(:,:,m,i) = DataTmp * diag(GAMMA2(i,:)) * DataTmp' + eye(model.nbVar) * diagRegularizationFactor;
        end
    end
    

    % Update B, Sigma Diagonal matrix    
    GAMMA_T = sum(GAMMA,2);      
    
	%Update SigmaDiag
	for i=1:model.nbStates
		for m=1:model.nbFrames
			model.SigmaDiag(:,:,m,i) = diag(diag(model.B(:,:,m)*model.S(:,:,m,i)*model.B(:,:,m)'));    %Eq.(9)
		end
	end
	
	for j=1:model.nbVar
		for m=1:model.nbFrames
			C = pinv(model.B(:,:,m)')*det(model.B(:,:,m));                 %C = cof(model.B);  %Eq.(6)
			% 				G = zeros(model.nbVar);    // readable computation of G
			% 				for i=1:model.nbStates
			% 					G = G + model.S(:,:,i) * sum(GAMMA(i,:),2) / model.SigmaDiag(k,k,i); %Eq.(7)
			% 				end
			G = sum(reshape(kron(squeeze(1/(model.SigmaDiag(j,j,m,:)))',ones(model.nbVar)).* reshape(squeeze(model.S(:,:,m,:)),[model.nbVar model.nbVar*model.nbStates]).* kron(GAMMA_T',ones(model.nbVar)),[model.nbVar model.nbVar model.nbStates]),3);
			model.B(j,:,m) = C(j,:)*pinv(G)*(sqrt(sum(sum(GAMMA,2)/(C(j,:)*pinv(G)*C(j,:)'))));  %Eq.(5)
		end
	end
    
    for m=1:model.nbFrames
        model.H(:,:,m) = pinv(model.B(:,:,m)) + eye(model.nbVar) * diagRegularizationFactor;
    end
        
    for k=1:model.nbStates
        for m=1:model.nbFrames
            model.Sigma(:,:,m,k) = alpha*(model.H(:,:,m)*model.SigmaDiag(:,:,m,k)*model.H(:,:,m)') + (1 - alpha)*model.S(:,:,m,k);  %Eq.(10)
        end
    end
    
	%Compute average log-likelihood
	LL(nbIter) = sum(log(sum(L,1))) / nbData;
	%Stop the algorithm if EM converged (small change of LL)
	if nbIter>nbMinSteps
		if LL(nbIter)-LL(nbIter-1)<maxDiffLL || nbIter==nbMaxSteps-1
%             figure;plot(LL); title('Likelihood Plot'); xlabel('Iterations'); ylabel('log-likelihood')
        	disp(['EM converged after ' num2str(nbIter) ' iterations.']);
			return;
		end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lik, GAMMA, GAMMA0] = computeGamma(Data, model)

nbData = size(Data, 3);
Lik = ones(model.nbStates, nbData);
GAMMA0 = zeros(model.nbStates, model.nbFrames, nbData);
for i=1:model.nbStates
	for m=1:model.nbFrames
		DataMat(:,:) = Data(:,m,:); %Matricization/flattening of tensor
		GAMMA0(i,m,:) = gaussPDF(DataMat, model.Mu(:,m,i), model.Sigma(:,:,m,i));
		Lik(i,:) = Lik(i,:) .* squeeze(GAMMA0(i,m,:))';
	end
	Lik(i,:) = Lik(i,:) * model.Priors(i);
end
GAMMA = Lik ./ repmat(sum(Lik,1)+realmin, size(Lik,1), 1);
end
