function [model, GAMMA] = EM_tensorHSMM(s, model)
%Estimation of TP-HSMM parameters with an EM algorithm
%
% If you found this work useful, please cite
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
% Written by Ajay Tanwani and Sylvain Calinon
%%
%Parameters of the EM algorithm
nbMinSteps = 5; %Minimum number of iterations allowed
nbMaxSteps = 50; %Maximum number of iterations allowed
maxDiffLL = 1E-4; %Likelihood increase threshold to stop the algorithm

diagRegularizationFactor = 1E-4; %Regularization term is optional, see Eq. (2.1.2) in doc/TechnicalReport.pdf

%Initialization of the parameters
nbSamples = length(s);

Data=[];
for n=1:nbSamples
	Data = cat(3, Data, s(n).Data);
end
nbData = size(Data,3);

minSigmaPd = 1E-6;
% model.Trans = rand(model.nbStates,model.nbStates);
% model.Trans = model.Trans ./ repmat(sum(model.Trans,2),1,model.nbStates);
% model.StatesPriors = rand(model.nbStates,1);
% model.StatesPriors = model.StatesPriors/sum(model.StatesPriors);

%Left-right model initialization
model.Trans = zeros(model.nbStates);
for i=1:model.nbStates-1
	model.Trans(i,i) = 1-(model.nbStates/nbData);
	model.Trans(i,i+1) = model.nbStates/nbData;
end
model.Trans(model.nbStates,model.nbStates) = 1.0;
model.StatesPriors = zeros(model.nbStates,1);
model.StatesPriors(1) = 1;


for nbIter=1:nbMaxSteps
	fprintf('.');
	
	%E-step
	for n=1:nbSamples
		
		%Emission probabilities
		s(n).B = computeGamma(s(n).Data, model); %See 'computeGamma' function below
		
		%Forward variable ALPHA
		s(n).ALPHA(:,1) = model.StatesPriors .* s(n).B(:,1);
		%Scaling to avoid underflow issues
		s(n).c(1) = 1 / sum(s(n).ALPHA(:,1)+realmin);
		s(n).ALPHA(:,1) = s(n).ALPHA(:,1) * s(n).c(1);
		for t=2:s(n).nbData
			s(n).ALPHA(:,t) = (s(n).ALPHA(:,t-1)'*model.Trans)' .* s(n).B(:,t); 
			%Scaling to avoid underflow issues
			s(n).c(t) = 1 / sum(s(n).ALPHA(:,t)+realmin);
			s(n).ALPHA(:,t) = s(n).ALPHA(:,t) * s(n).c(t);
		end
		
		%Backward variable BETA
		s(n).BETA(:,s(n).nbData) = ones(model.nbStates,1) * s(n).c(end); %Rescaling
		for t=s(n).nbData-1:-1:1
			s(n).BETA(:,t) = model.Trans * (s(n).BETA(:,t+1) .* s(n).B(:,t+1));
			s(n).BETA(:,t) = min(s(n).BETA(:,t) * s(n).c(t), realmax); %Rescaling
		end
		
		%Intermediate variable GAMMA
		s(n).GAMMA = (s(n).ALPHA.*s(n).BETA) ./ repmat(sum(s(n).ALPHA.*s(n).BETA)+realmin, model.nbStates, 1); 
		
		%Intermediate variable XI (fast version, by considering scaling factor)
		for i=1:model.nbStates
			for j=1:model.nbStates
				s(n).XI(i,j,:) = model.Trans(i,j) * (s(n).ALPHA(i,1:end-1) .* s(n).B(j,2:end) .* s(n).BETA(j,2:end)); 
			end
		end
	end
	
	%Concatenation of HMM intermediary variables
	GAMMA=[]; GAMMA_TRK=[]; GAMMA_INIT=[]; XI=[];
	for n=1:nbSamples
		GAMMA = [GAMMA s(n).GAMMA];
		GAMMA_INIT = [GAMMA_INIT s(n).GAMMA(:,1)];
		GAMMA_TRK = [GAMMA_TRK s(n).GAMMA(:,1:end-1)];
		XI = cat(3,XI,s(n).XI);
	end
	GAMMA2 = GAMMA ./ repmat(sum(GAMMA,2)+realmin, 1, size(GAMMA,2));
    model.Pix = GAMMA2;
    model.PiK = GAMMA;
    
	%M-step
	for i=1:model.nbStates	
		for m=1:model.nbFrames
			%Matricization/flattening of tensor
			DataMat(:,:) = Data(:,m,:);
			
			%Update Mu
			model.Mu(:,m,i) = DataMat * GAMMA2(i,:)';
			
			%Update Sigma (regularization term is optional), see Eqs (2.5.9) and (6.0.4) in doc/TechnicalReport.pdf
			DataTmp = DataMat - repmat(model.Mu(:,m,i),1,nbData);
			model.Sigma(:,:,m,i) = DataTmp * diag(GAMMA2(i,:)) * DataTmp' + eye(model.nbVar) * diagRegularizationFactor;
		end
	end
	
	%Update initial state probability vector
	model.StatesPriors = mean(GAMMA_INIT,2); 
	
	%Update transition probabilities
	model.Trans = sum(XI,3)./ repmat(sum(GAMMA_TRK,2)+realmin, 1, model.nbStates); 
	
	%Compute the average log-likelihood through the ALPHA scaling factors
	LL(nbIter)=0;
	for n=1:nbSamples
		LL(nbIter) = LL(nbIter) - sum(log(s(n).c));
	end
	LL(nbIter) = LL(nbIter)/nbSamples;
	%Stop the algorithm if EM converged
	if nbIter>nbMinSteps
		if LL(nbIter)-LL(nbIter-1)<maxDiffLL
			disp(['EM converged after ' num2str(nbIter) ' iterations.']);
			break;
		end
	end
end
% disp(['The maximum number of ' num2str(nbMaxSteps) ' EM iterations has been reached.']);


%Removal of self-transition and normalization (for HSMM representation) 
model.Trans = model.Trans - diag(diag(model.Trans)) + eye(model.nbStates)*1E-12;
model.Trans = model.Trans ./ repmat(sum(model.Trans,2),1,model.nbStates);

for i=1:model.nbStates
	st(i).d=[];
end
[~,hmax] = max(GAMMA);
currState = hmax(1);
cnt = 1;
for t=1:length(hmax)
	if (hmax(t)==currState)
		cnt = cnt+1;
	else
		st(currState).d = [st(currState).d cnt];
		cnt = 1;
		currState = hmax(t);
	end
end
st(currState).d = [st(currState).d cnt];

%Compute state duration as Gaussian distribution
for i=1:model.nbStates
	model.Mu_Pd(1,i) = mean(st(i).d);
	model.Sigma_Pd(1,1,i) = cov(st(i).d) + minSigmaPd;
end

%Compute resulting state duration density
nbData2 = round(0.5 * nbData/model.nbStates); %Number of maximum duration step to consider in the HSMM (2 is a safety factor)
for i=1:model.nbStates
	model.Pd(i,:) = gaussPDF([1:nbData2], model.Mu_Pd(:,i), model.Sigma_Pd(:,:,i));
	%model.Pd(i,:) = model.Pd(i,:) / max(model.Pd(i,:));
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
