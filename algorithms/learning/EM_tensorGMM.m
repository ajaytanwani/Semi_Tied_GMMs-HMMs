function [model, GAMMA0] = EM_tensorGMM(Data, model)
% Training of a task-parameterized Gaussian mixture model (GMM) with an expectation-maximization (EM) algorithm.
% The approach allows the modulation of the centers and covariance matrices of the Gaussians with respect to
% external parameters represented in the form of candidate coordinate systems.
%
% This source code is given for free! In exchange, I would be grateful if you cite
% the following reference in any academic publication that uses this code or part of it:
%
% @inproceedings{Calinon14ICRA,
%   author="Calinon, S. and Bruno, D. and Caldwell, D. G.",
%   title="A task-parameterized probabilistic model with minimal intervention control",
%   booktitle="Proc. {IEEE} Intl Conf. on Robotics and Automation ({ICRA})",
%   year="2014",
%   month="May-June",
%   address="Hong Kong, China",
%   pages=""
% }
% Copyright (c) 2016 Idiap Research Institute
% Written by Sylvain Calinon

%%
%Parameters of the EM algorithm
nbMinSteps = 5; %Minimum number of iterations allowed
nbMaxSteps = 100; %Maximum number of iterations allowed
maxDiffLL = 1E-4; %Likelihood increase threshold to stop the algorithm
nbData = size(Data,3);

diagRegularizationFactor = 1E-4; %Regularization term is optional

for nbIter=1:nbMaxSteps
	fprintf('.');
	
	%E-step
	[L, GAMMA, GAMMA0] = computeGamma(Data, model); %See 'computeGamma' function below
	GAMMA2 = GAMMA ./ repmat(sum(GAMMA,2),1,nbData);
	model.Pix = GAMMA2;
    model.PiK = GAMMA;
	
	%M-step
	for i=1:model.nbStates
		
		%Update Priors
		model.Priors(i) = sum(sum(GAMMA(i,:))) / nbData;
		
		for m=1:model.nbFrames
			%Matricization/flattening of tensor
			DataMat(:,:) = Data(:,m,:);
			
			%Update Mu
			model.Mu(:,m,i) = DataMat * GAMMA2(i,:)';
			
			%Update Sigma (regularization term is optional)
			DataTmp = DataMat - repmat(model.Mu(:,m,i),1,nbData);
			model.Sigma(:,:,m,i) = DataTmp * diag(GAMMA2(i,:)) * DataTmp' + eye(model.nbVar) * diagRegularizationFactor;
		end
	end
	
	%Compute average log-likelihood
	LL(nbIter) = sum(log(sum(L,1))) / size(L,2);
	%Stop the algorithm if EM converged (small change of LL)
	if nbIter>nbMinSteps
		if LL(nbIter)-LL(nbIter-1)<maxDiffLL || nbIter==nbMaxSteps-1
            figure;plot(LL); title('Likelihood Plot'); xlabel('Iterations'); ylabel('log-likelihood')
			disp(['EM converged after ' num2str(nbIter) ' iterations.']);
			return;
		end
	end
end
disp(['The maximum number of ' num2str(nbMaxSteps) ' EM iterations has been reached.']);
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