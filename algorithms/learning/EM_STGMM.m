function [model LL] = EM_STGMM(Data, model)
%%% EM procedure for training of Semi-Tied Gaussian Mixture Model
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
diagRegularizationFactor = 1E-6;
B_sf = 1E-1;
alpha = 0.99;

%%% overwrite parameters if defined in model structure
if isfield(model,'nbMinSteps'), nbMinSteps = model.nbMinSteps; disp(['nbMinSteps = ', num2str(nbMinSteps)]); end
if isfield(model,'nbMaxSteps'), nbMaxSteps = model.nbMaxSteps; disp(['nbMaxSteps = ', num2str(nbMaxSteps)]); end
if isfield(model,'maxDiffLL'), maxDiffLL = model.maxDiffLL; disp(['maxDiffLL = ', num2str(maxDiffLL)]); end
if isfield(model,'diagRegularizationFactor'), diagRegularizationFactor = model.diagRegularizationFactor; disp(['diagRegularizationFactor = ', num2str(diagRegularizationFactor)]); end

if isfield(model,'B_sf'), B_sf = model.B_sf; disp(['B_sf = ', num2str(B_sf)]);  end
if isfield(model,'alpha'), alpha = model.alpha; disp(['alpha = ', num2str(alpha)]); end
    
nbData = size(Data,2);

if ~isfield(model,'B')
    model.B = eye(model.nbVar)*B_sf;
    model.InitH = pinv(model.B) + eye(model.nbVar) * diagRegularizationFactor;

    % model.Sigma = repmat(GenerateRandomSymmetricMatrix(model.nbVar),[1 1 model.nbStates]);
    for i=1:model.nbStates
    %     model.InitSigmaDiag(:,:,i) = diag(diag(model.B*squeeze(model.Sigma(:,:,i))*model.B'));
        [xx model.InitSigmaDiag(:,:,i)] = eig(squeeze(model.Sigma(:,:,i)));
    end
end

for nbIter=1:nbMaxSteps        	
	%E-step
	[L, GAMMA] = computeGamma(Data, model); %See 'computeGamma' function below
	GAMMA2 = GAMMA ./ repmat(sum(GAMMA,2),1,nbData);
    model.Pix = GAMMA2;

    %M-step - Update Priors, Mu, Sample Covariance
    for i=1:model.nbStates
        %Update Priors
        model.Priors(i) = sum(GAMMA(i,:)) / nbData;
        
        %Update Mu
        model.Mu(:,i) = Data * GAMMA2(i,:)';
        
        %Update Sigma
        DataTmp = Data - repmat(model.Mu(:,i),1,nbData);
        model.S(:,:,i) = DataTmp * diag(GAMMA2(i,:)) * DataTmp' + eye(model.nbVar) * diagRegularizationFactor;
    end

    % Update B matrix    
    GAMMA_T = sum(GAMMA,2);  
	
	for i=1:model.nbStates
		model.SigmaDiag(:,:,i) = diag(diag(model.B*squeeze(model.S(:,:,i))*model.B')); %Eq.(9)
	end
	
	for j=1:model.nbVar
		C = pinv(model.B')*det(model.B);                 %C = cof(model.B); Eq.(6)
		% 				G = zeros(model.nbVar);    // readable computation of G %Eq.(7)
		% 				for i=1:model.nbStates
		% 					G = G + model.S(:,:,i) * sum(GAMMA(i,:),2) / model.SigmaDiag(k,k,i);
		% 				end
		G = sum(reshape(kron(squeeze(1/(model.SigmaDiag(j,j,:)))',ones(model.nbVar)).* reshape(model.S,[model.nbVar model.nbVar*model.nbStates]).* kron(GAMMA_T',ones(model.nbVar)),[model.nbVar model.nbVar model.nbStates]),3);
		
		model.B(j,:) = C(j,:)*pinv(G)*(sqrt(sum(sum(GAMMA,2)/(C(j,:)*pinv(G)*C(j,:)'))));   %Eq.(5)
	end
    
    model.H = pinv(model.B) + eye(model.nbVar) * diagRegularizationFactor;
    for k=1:model.nbStates
        model.Sigma(:,:,k) = alpha*(model.H*model.SigmaDiag(:,:,k)*model.H') + (1 - alpha)*model.S(:,:,k); %Eq.(10)
    end
    
	%Compute average log-likelihood
	LL(nbIter) = sum(log(sum(L,1))) / nbData;
    
    disp([num2str(nbIter) '. ' num2str(LL(nbIter))]);
    
	%Stop the algorithm if EM converged (small change of LL)
	if nbIter>nbMinSteps
		if LL(nbIter)-LL(nbIter-1)<maxDiffLL || nbIter==nbMaxSteps-1
            figure;plot(LL); title('Likelihood Plot'); xlabel('Iterations'); ylabel('log-likelihood')
        	disp(['EM converged after ' num2str(nbIter) ' iterations.']);
			return;
		end
    end
end
end

function [L, GAMMA] = computeGamma(Data, model)
L = zeros(model.nbStates,size(Data,2));

for i=1:model.nbStates
    L(i,:) = model.Priors(i) * gaussPDF(Data, model.Mu(:,i), model.Sigma(:,:,i));
end
    
GAMMA = L ./ repmat(sum(L,1)+realmin, model.nbStates, 1);
end
