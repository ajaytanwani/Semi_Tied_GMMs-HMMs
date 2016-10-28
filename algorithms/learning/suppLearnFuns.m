classdef suppLearnFuns
%%% This class defines supplementary learning functions for sampling
%%% from HSMMs, and product of Gaussians
% Copyright (c) 2016 Idiap Research Institute
% Written by Sylvain Calinon, Ajay Tanwani

    methods (Static)
	%%
	function [bmx, ALPHA, S, alpha] = hsmm_fwd_init_hsum(Data, m)
		%Forward variable initialization
		
		% Check if we need to evaluate marginal probability:
		if nargin ==2
			in = 1:size(m.Mu,1);
		else
			in = varargin{1};
		end
		
		for i=1:size(m.Sigma,3)
			Btmp(i,1) = gaussPDF(Data(in,:), m.Mu(in,i), m.Sigma(in,in,i)) +1e-12;
		end
		Btmp = Btmp / sum(Btmp);
		
		ALPHA = repmat(m.StatesPriors, 1, size(m.Pd,2)) .* m.Pd;	
		r = Btmp' * sum(ALPHA,2);    
		bmx(:,1) = Btmp ./ r;        
		E = bmx .* ALPHA(:,1);       
		S = m.Trans' * E;        
		
		alpha = Btmp .* sum(ALPHA,2);		
	end
	
	%%
	function [ALPHA, S, alpha] = hsmm_fwd_init_ts(model)
		%Forward variable by using only temporal/sequential information from Trans and Pd
		
		ALPHA = repmat(model.StatesPriors, 1, size(model.Pd,2)) .* model.Pd; %Equation (13)
		S = model.Trans' * ALPHA(:,1);	%Equation (6)
		alpha = sum(ALPHA,2); %Forward variable
	end
	%%
	function [bmx, ALPHA, S, alpha] = hsmm_fwd_step_hsum(Data, m, bmx, ALPHA, S)
		%Forward variable iteration
		
		% Check if we need to evaluate marginal probability:
		if nargin ==5
			%  Input for evaluation of GMM not specified, use all variables
			in = 1:size(m.Mu,1);
		else
			in =varargin{1};
		end
		
		nbD = size(m.Pd, 2);
		nbStates = size(m.Sigma, 3);
		
		for i=1:nbStates
			Btmp(i,1) = gaussPDF(Data(in), m.Mu(in,i), m.Sigma(in,in,i))+1e-12;
		end
		Btmp = Btmp/sum(Btmp);
		
		%Fast computation
		ALPHA = [repmat(S(:,end), 1, nbD-1) .* m.Pd(:,1:nbD-1) + repmat(bmx(:,end), 1, nbD-1) .* ALPHA(:,2:nbD), ...
			S(:,end) .* m.Pd(:,nbD)];	%Equation (12)
		
		% %Slow computation (but easier to read)
		% for i=1:nbD-1
		%   ALPHA(:,i) = S(:,end) .* m.Pd(:,i) + bmx(:,end) .* ALPHA(:,i+1);	%Equation (12)
		% end
		% ALPHA(:,nbD) = S(:,end) .* m.Pd(:,nbD);
		%
		% ALPHA_SUM = zeros(nbStates,1);
		% for i=1:nbStates
		%   ALPHA_SUM(i) = sum(ALPHA(i,:));
		% end
		
		r = Btmp' * sum(ALPHA,2);     %Equation (3)
		bmx = [bmx, Btmp ./ r];	      %Equation (2)
		E = bmx(:,end) .* ALPHA(:,1); %Equation (5)
		S = [S,m.Trans' * E];         %Equation (6)
		
		alpha = Btmp .* sum(ALPHA,2); %Forward variable
		alpha = alpha / sum(alpha);
	end
	
	%%
	function [ALPHA, S, alpha] = hsmm_fwd_step_ts(model, ALPHA, S)
		%Forward variable by using only temporal/sequential information from Trans and Pd
		
		t = size(ALPHA,2) + 1;
		nbD = size(model.Pd, 2);
		
		%Fast computation
		ALPHA = [repmat(S(:,end), 1, nbD-1) .* model.Pd(:,1:nbD-1) + ALPHA(:,2:nbD), ...
			S(:,end) .* model.Pd(:,nbD)]; %Equation (12)
		
		% %Slow computation (but more readable)
		% for i=1:nbD-1
		%   ALPHA(:,i) = S(:,end) .* model.Pd(:,i) + ALPHA(:,i+1);	%Equation (12)
		% end
		% ALPHA(:,nbD) = S(:,end) .* model.Pd(:,nbD);
		
		S = [S, model.Trans' * ALPHA(:,1)];	%Equation (6)
		alpha = sum(ALPHA, 2); %Forward variable
	end
	
	%%
	function H = Likelihood(x,model,in)
		for i=1:model.nbStates
			H(i) = model.Priors(i) * gaussPDF(x, model.Mu(in,i), model.Sigma(in,in,i));
		end
		H = H/sum(H);
	end
	%%
	function [Mu, Sigma] = productTPGMM0(model, p)
		%Compute the product of Gaussians for a TP-GMM, where the set of task parameters are stored in the variable 'p'.
		%Sylvain Calinon, 2015
		
		for i=1:model.nbStates
			% Reallocating
			SigmaTmp = zeros(model.nbVar);
			MuTmp = zeros(model.nbVar,1);
			% Product of Gaussians
			for m=1:model.nbFrames
				MuP = p(m).A * model.Mu(:,m,i) + p(m).b;
				SigmaP = p(m).A * model.Sigma(:,:,m,i) * p(m).A';
				SigmaTmp = SigmaTmp + inv(SigmaP);
				MuTmp = MuTmp + SigmaP\MuP;
			end
			Sigma(:,:,i) = inv(SigmaTmp);
			Mu(:,i) = Sigma(:,:,i) * MuTmp;
		end
	end
	%%
	function [q,r, alphaEst] = stateSequence(model,r,Observation,nbData,nbHist)
		
		% Input:
		% Model 	: HMM Model
		% r			: record of changes over time (include Alpha, bmx, S)
		% Observation: Current observation
		% nbData	: Number of samples to create
		
		% Outputs:
		% q			: Generated State sequence		[1 X nbData]
		% r			: HMM track record				[struct]
		% postH		: Posterior value of H for this time step
		% h			: Stochastic sampled posteriors [nbStates X nbData]
		
		if isfield(r,'ALPHA')==0
			% Initialize
			[r.bmx, r.ALPHA, r.S, r.alpha]= suppLearnFuns.hsmm_fwd_init_hsum(Observation, model);
		else
			% Otherwise
			[r.bmx, r.ALPHA, r.S, alphaNew] = suppLearnFuns.hsmm_fwd_step_hsum(Observation, model, r.bmx, r.ALPHA, r.S);
			r.alpha = [r.alpha,alphaNew];
		end
		
		%% First trial Stochastic sampling
		% Set current state distribution as the priors to start generation
		% postH = postH./sum(postH);  % Posterior mixing component
		% model.StatesPriors = postH; % Initialize state priors for stochastic sampling
		
		% % Perform stochastic sampling for the next q states
		% nbSt=0; currTime=0; iList=[];
		% h = zeros(model.nbStates,nbData);
		% while currTime<nbData
		% 		nbSt = nbSt+1;
		% 		if nbSt==1
		% 			[~,iList(1)] = max(model.StatesPriors.*rand(model.nbStates,1));
		% 			h1 = ones(1,nbData);
		% 		else
		% 			h1 = [zeros(1,currTime), cumsum(model.Pd(iList(end-1),:)), ones(1,max(nbData-currTime-nbD,0))];
		% 			currTime = currTime + round(model.Mu_Pd(1,iList(end-1)));
		% 		end
		% 		h2 = [ones(1,currTime), 1-cumsum(model.Pd(iList(end),:)), zeros(1,max(nbData-currTime-nbD,0))];
		% 		h(iList(end),:) = h(iList(end),:) + min([h1(1:nbData); h2(1:nbData)]);
		% 		[~,iList(end+1)] = max(model.Trans(iList(end),:).*rand(1,model.nbStates));
		% end
		% h = h ./ repmat(sum(h,1),model.nbStates,1); % Stochastic Mixing components
		
		%% Second trial: Forward variable without observations
		% Generate future observations:
		
		% Copy record values to temporary arrays. We the alpha's we are calculating
		% are future estimates not based on position data.
		ALPHATmp = r.ALPHA; % Copy Alpha to tmp, we don't want to include the future predicitons in our record
		Stmp	 = r.S;		% Copy S to tmp, we don't want to include the future predicitons in our record
		
		alphaEst = zeros(model.nbStates,nbData);
		for n =1:nbData
			[ALPHATmp, Stmp,alphaEst(:,n)]	= suppLearnFuns.hsmm_fwd_step_ts(model, ALPHATmp, Stmp);
		end
		
		alphaEst = alphaEst./repmat(sum(alphaEst,1),model.nbStates,1);
		
		%% Third trial: Forward variable with stochastic sampling
		
		%% Construct state sequence
		[~,Hsize] = size(r.alpha);
		nbHistVal = min([nbHist,(Hsize)]);			% Number of states from history used
		indH	  = (Hsize-nbHistVal+1):Hsize;		% History indices to use
		indF	  = 1:(nbData-length(indH));		% Future indices to use
		[~,q]	  = max([r.alpha(:,indH), alphaEst(:,indF)],[],1);			% construct state sequence
		
	end
	
	end	
end