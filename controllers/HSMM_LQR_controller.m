function rep = HSMM_LQR_controller(s, model)
%%Reconstruct GMM for each demonstration, sample the HSMM sequence and track the stepwise sequence
%with LQR

%%% Author: Ajay Tanwani, 2016 <http://www.ajaytanwani.com>
%
% If you found this work useful, please cite
%
% @ARTICLE{Tanwani2016a,
% author={Tanwani, A.K. and Calinon, S.},
% journal={Robotics and Automation Letters, IEEE},
% title={Learning Robot Manipulation Tasks With Task-Parameterized Semitied Hidden Semi-Markov Model},
% year={2016},
% volume={1},
% number={1},
% pages={235-242}}

% Copyright (c) 2016 Idiap Research Institute
% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>

%%

for n=1:size(s,2)
	[rep(n).Mu, rep(n).Sigma] = suppLearnFuns.productTPGMM0(model, s(n).p);

    rep(n).Priors = model.Priors; 
    rep(n).nbStates = model.nbStates; 
    rep(n).nbVarPos = model.nbVarPos;
    rep(n).nbVar = model.nbVar;
    rep(n).StatesPriors = model.StatesPriors; 
    rep(n).Trans = model.Trans;
    rep(n).Mu_Pd = model.Mu_Pd; 
    rep(n).Sigma_Pd = model.Sigma_Pd; 
    rep(n).Pd = model.Pd; 
    rep(n).nbD = model.nbD; 
    rep(n).dt = model.dt;
    rep(n).p = s(n).p;
    
	%Generate starting point
	ind = randi([1,rep(n).nbStates]); % Select random state
% 	X = [model.Mu(1:rep(n).nbVarPos,ind)+randn(rep(n).nbVarPos,1)*5E0; zeros(rep(n).nbVarPos,1)];
%     X = [rep(n).Mu(1:rep(n).nbVar,ind)+randn(rep(n).nbVar,1)*0E0];
    X = [rep(n).Mu(1:rep(n).nbVar,1)+randn(rep(n).nbVar,1)*0E0];
    rep(n).X = X;
    rep(n).ind = ind;
%     X = s(n).Data1(:,1);
    
	%Find maximum likelihood state
	rep(n).StatesPriors = suppLearnFuns.Likelihood(X(1:rep(n).nbVar),rep(n),1:rep(n).nbVar)';
	
	%Get state sequence
	[rep(n).q tmp rep(n).alpha] = suppLearnFuns.stateSequence(rep(n),[],X,rep(n).nbD,1);		
        
	%Build stepwise reference trajectory    
    a.currTar = rep(n).Mu(:,rep(n).q);
    a.currSigma = rep(n).Sigma(:,:,rep(n).q);
	
    traj = TV_LQR_continuous_ff(rep(n), a, X, model.rfactor);
     
    rep(n).currTar = a.currTar;
    rep(n).currSigma = a.currSigma;
    rep(n).Data = traj.Data;
    rep(n).dData = traj.dData;
    rep(n).ddxNorm = traj.ddxNorm;
    rep(n).kpDet = traj.kpDet;
    rep(n).kvDet = traj.kvDet;   
    rep(n).P = traj.P;      
    rep(n).expSigma = traj.K;
end


