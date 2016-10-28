function [expData, expSigma, H] = GMR(model, DataIn, in, out)
% Gaussian mixture regression (GMR)
%
% Copyright (c) 2016 Idiap Research Institute
% Written by Sylvain Calinon

nbData = size(DataIn,2);
nbVarOut = length(out);

diagRegularizationFactor = 1E-8; %Regularization term is optional

MuTmp = zeros(nbVarOut,model.nbStates);
expData = zeros(nbVarOut,nbData);
expSigma = zeros(nbVarOut,nbVarOut,nbData);
for t=1:nbData
	
	% Compute activation weight
	for i=1:model.nbStates
		H(i,t) = model.Priors(i) * gaussPDF(DataIn(:,t), model.Mu(in,i), model.Sigma(in,in,i));
	end
	H(:,t) = H(:,t) / sum(H(:,t)+realmin);
	
	% Compute conditional means
	for i=1:model.nbStates
		MuTmp(:,i) = model.Mu(out,i) + model.Sigma(out,in,i)/model.Sigma(in,in,i) * (DataIn(:,t)-model.Mu(in,i));
		expData(:,t) = expData(:,t) + H(i,t) * MuTmp(:,i);
	end
	
	% Compute conditional covariances
	for i=1:model.nbStates
		SigmaTmp = model.Sigma(out,out,i) - model.Sigma(out,in,i)/model.Sigma(in,in,i) * model.Sigma(in,out,i);
		expSigma(:,:,t) = expSigma(:,:,t) + H(i,t) * (SigmaTmp + MuTmp(:,i)*MuTmp(:,i)');
	end
	expSigma(:,:,t) = expSigma(:,:,t) - expData(:,t)*expData(:,t)' + eye(nbVarOut) * diagRegularizationFactor; 
end

