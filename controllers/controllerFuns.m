classdef controllerFuns
%%% This class is a wrapper for calling decoding controllers.
% Copyright (c) 2016 Idiap Research Institute
% Written by Ajay Tanwani

    methods (Static)
	%%
	function [rep] = GMR_controller(Data, model)
		%%% Generate output trajectory by decoding the model with GMR
		
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
		
		disp('Reproductions with GMR...');
		
		for n=1:1
			[rep(n).Data, rep(n).expSigma] = GMR(model, Data(1,:), 1, 2:size(Data,1));
			AbsErr = abs(rep(n).Data - Data(2:end,:)).^2;
			rep(n).MSE = sum(AbsErr(:))/numel(rep(n).Data);
		end
	end
	
	%%
	
	function [r] = TP_GMR_controller(s, model)
		%%% Generate output trajectory by decoding the model with GMR
		
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
		
		disp('Reproductions with GMR...');
		DataIn = s(n).Data0(1,:);
		for n=1:nbSamples
			[r(n).Mu, r(n).Sigma] = productTPGMM0(model, s(n).p); %See Eq. (6.0.5), (6.0.6) and (6.0.7) in doc/TechnicalReport.pdf
			r(n).Priors = model.Priors;
			r(n).nbStates = model.nbStates;
			r(n).Data = [DataIn; GMR(r(n), DataIn, 1, 2:model.nbVar)]; %See Eq. (3.0.2) to (3.0.5) in doc/TechnicalReport.pdf
		end
	end
	
	%%
	function rep = LQR_controller(Data, model)
		%%% Generate output trajectory by tracking the stepwise sequence with LQR
		
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
		
		%
		model.rfactor = 9E1;
		% nbD = size(Data,2);
		for n=1:model.nbSamples
			
			rep(n).Priors = model.Priors;
			rep(n).nbStates = model.nbStates;
			rep(n).nbVarPos = model.nbVarPos;
			rep(n).nbVar = model.nbVar;
			rep(n).nbD = model.nbD;
			rep(n).dt = model.dt;
			
			%Get state sequence
			[~,rep(n).q] = max(model.Pix(:,(n-1)*model.nbD+1:n*model.nbD),[],1); %works also for nbStates=1
			
			%Build stepwise reference trajectory
			a.currTar = model(n).Mu(:,rep(n).q);
			a.currSigma = model(n).Sigma(:,:,rep(n).q);
			
			X = a.currTar(:,1);
			% 	%LQR tracking (continuous version)
			
			%     traj = TV_LQR_discrete(model, a, X, model.rfactor);
			%     traj = TV_LQR_continuous(model, a, X, model.rfactor);
			%     traj = TV_LQR_continuous_ff(rep(n), a, X, model.rfactor);
			traj = TV_LQT_continuous(rep(n), a, X, model.rfactor);
			%     traj = Infinite_LQR_continuous_ff(rep(n), a, X, model.rfactor);
			
			
			rep(n).currTar = a.currTar;
			rep(n).currSigma = a.currSigma;
			rep(n).Data = traj.Data;
			rep(n).ddxNorm = traj.ddxNorm;
			rep(n).kpDet = traj.kpDet;
			rep(n).kvDet = traj.kvDet;
			rep(n).P = traj.P;
			
		end
		
	end
	
	%%
	function rep = TP_LQR_controller(s, model)
		%%% Generate output trajectory by decoding the model with LQR. %Reconstruct GMM for each demonstration and track the stepwise sequence
		%with LQR
		
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
				
		model.rfactor = 9E1;
		
		for n=1:model.nbSamples
			
			[rep(n).Mu, rep(n).Sigma] = productTPGMM0(model, s(n).p);
			
			rep(n).Priors = model.Priors;
			rep(n).nbStates = model.nbStates;
			rep(n).nbVarPos = model.nbVarPos;
			rep(n).nbVar = model.nbVar;
			rep(n).nbD = model.nbD;
			rep(n).dt = model.dt;
			rep(n).p = s(n).p;
			
			
			%Get state sequence
			[~,rep(n).q] = max(model.Pix(:,(n-1)*model.nbD+1:n*model.nbD),[],1); %works also for nbStates=1
			
			%Build stepwise reference trajectory
			a.currTar = rep(n).Mu(:,rep(n).q);
			a.currSigma = rep(n).Sigma(:,:,rep(n).q);
			
			X = a.currTar(:,1);
			% 	%LQR tracking (continuous version)
			%     traj = TV_LQR_discrete(model, a, X, model.rfactor);
			%     traj = TV_LQR_continuous(model, a, X, model.rfactor);
			traj = TV_LQR_continuous_ff(rep(n), a, X, model.rfactor);
			
			rep(n).currTar = a.currTar;
			rep(n).currSigma = a.currSigma;
			rep(n).Data = traj.Data;
			rep(n).ddxNorm = traj.ddxNorm;
			rep(n).kpDet = traj.kpDet;
			rep(n).kvDet = traj.kvDet;
			rep(n).P = traj.P;
			
		end
	end
	end	
end