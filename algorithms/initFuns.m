classdef initFuns
%%% This class defines different initialization methods based on either time-based
%%% or kmeans based clustering
%
% Copyright (c) 2016 Idiap Research Institute
% Written by Ajay Tanwani

    methods (Static)
	%%
	function model = init_GMM_kmeans(Data, model)
		%
		% This function initializes the parameters of a Gaussian Mixture Model
		% (GMM) by using k-means clustering algorithm.
		%
		% Inputs -----------------------------------------------------------------
		%   o Data:     D x N array representing N datapoints of D dimensions.
		%   o nbStates: Number K of GMM components.
		% Outputs ----------------------------------------------------------------
		%   o Priors:   1 x K array representing the prior probabilities of the
		%               K GMM components.
		%   o Mu:       D x K array representing the centers of the K GMM components.
		%   o Sigma:    D x D x K array representing the covariance matrices of the
		%               K GMM components.
		% Comments ---------------------------------------------------------------
		%   o This function uses the 'kmeans' function from the MATLAB Statistics
		%     toolbox. If you are using a version of the 'netlab' toolbox that also
		%     uses a function named 'kmeans', please rename the netlab function to
		%     'kmeans_netlab.m' to avoid conflicts.
		%
		% Copyright (c) 2016, Idiap Research Institute
		% Written by Sylvain Calinon, 2006
		
		[nbVar, nbData] = size(Data);
		diagRegularizationFactor = 1E-2;
		
		[Data_id, model.Mu] = kmeansClustering(Data, model.nbStates);
		
		for i=1:model.nbStates
			idtmp = find(Data_id==i);
			model.Priors(i) = length(idtmp);
			model.Sigma(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
			%Regularization term to avoid numerical instability
			model.Sigma(:,:,i) = model.Sigma(:,:,i) + eye(nbVar)*diagRegularizationFactor;
		end
		model.Priors = model.Priors / sum(model.Priors);
	end
	
	%%	
	function model = init_GMM_timeBased(Data, model)
		%
		% This function initializes the parameters of a Gaussian Mixture Model
		% (GMM) by using k-means clustering algorithm.
		%
		% Inputs -----------------------------------------------------------------
		%   o Data:     D x N array representing N datapoints of D dimensions.
		%   o nbStates: Number K of GMM components.
		% Outputs ----------------------------------------------------------------
		%   o Priors:   1 x K array representing the prior probabilities of the
		%               K GMM components.
		%   o Mu:       D x K array representing the centers of the K GMM components.
		%   o Sigma:    D x D x K array representing the covariance matrices of the
		%               K GMM components.
		% Comments ---------------------------------------------------------------
		%   o This function uses the 'kmeans' function from the MATLAB Statistics
		%     toolbox. If you are using a version of the 'netlab' toolbox that also
		%     uses a function named 'kmeans', please rename the netlab function to
		%     'kmeans_netlab.m' to avoid conflicts.
		%
		% Copyright (c) 2016, Idiap Research Institute
		% Written by Sylvain Calinon, 2006
		
		[nbVar, nbData] = size(Data);
		diagRegularizationFactor = 1E-2;
		
		TimingSep = linspace(min(Data(1,:)), max(Data(1,:)), model.nbStates+1);
		
		for i=1:model.nbStates
			idtmp = find( Data(1,:)>=TimingSep(i) & Data(1,:)<TimingSep(i+1));
			model.Priors(i) = length(idtmp);
			model.Mu(:,i) = mean(Data(:,idtmp)');
			model.Sigma(:,:,i) = cov(Data(:,idtmp)');
			%Regularization term to avoid numerical instability
			model.Sigma(:,:,i) = model.Sigma(:,:,i) + eye(nbVar)*diagRegularizationFactor;
		end
		model.Priors = model.Priors / sum(model.Priors);
	end
	%%
	
	function model = init_tensorGMM_kmeans(Data, model)
		% Initialization of the model with k-means.
		% Copyright (c) 2016, Idiap Research Institute
		% Written by Sylvain Calinon, Tohid Alizadeh, 2013
		
		diagRegularizationFactor = 1E-4;
		
		DataAll = reshape(Data, size(Data,1)*size(Data,2), size(Data,3)); %Matricization/flattening of tensor
		
		%k-means clustering
		[Data_id, Mu] = kmeansClustering(DataAll, model.nbStates);
		
		for i=1:model.nbStates
			idtmp = find(Data_id==i);
			model.Priors(i) = length(idtmp);
			Sigma(:,:,i) = cov([DataAll(:,idtmp) DataAll(:,idtmp)]') + eye(size(DataAll,1))*diagRegularizationFactor;
		end
		model.Priors = model.Priors / sum(model.Priors);
		
		%Reshape GMM parameters into a tensor
		for m=1:model.nbFrames
			for i=1:model.nbStates
				model.Mu(:,m,i) = Mu((m-1)*model.nbVar+1:m*model.nbVar,i);
				model.Sigma(:,:,m,i) = Sigma((m-1)*model.nbVar+1:m*model.nbVar,(m-1)*model.nbVar+1:m*model.nbVar,i);
			end
		end
	end
	
	%%
	function model = init_tensorGMM_timeBased(Data, model)
		% Copyright (c) 2016, Idiap Research Institute
		% Written by Sylvain Calinon, 2014
		
		diagRegularizationFactor = 1E-4;
		
		DataAll = reshape(Data, size(Data,1)*size(Data,2), size(Data,3)); %Matricization/flattening of tensor
		
		TimingSep = linspace(min(DataAll(1,:)), max(DataAll(1,:)), model.nbStates+1);
		Mu = zeros(model.nbFrames*model.nbVar, model.nbStates);
		Sigma = zeros(model.nbFrames*model.nbVar, model.nbFrames*model.nbVar, model.nbStates);
		for i=1:model.nbStates
			idtmp = find( DataAll(1,:)>=TimingSep(i) & DataAll(1,:)<TimingSep(i+1));
			Mu(:,i) = mean(DataAll(:,idtmp),2);
			Sigma(:,:,i) = cov(DataAll(:,idtmp)') + eye(size(DataAll,1))*diagRegularizationFactor;
			model.Priors(i) = length(idtmp);
		end
		model.Priors = model.Priors / sum(model.Priors);
		
		%Reshape GMM parameters into a tensor
		for m=1:model.nbFrames
			for i=1:model.nbStates
				model.Mu(:,m,i) = Mu((m-1)*model.nbVar+1:m*model.nbVar,i);
				model.Sigma(:,:,m,i) = Sigma((m-1)*model.nbVar+1:m*model.nbVar,(m-1)*model.nbVar+1:m*model.nbVar,i);
			end
		end
	end
	
	%%
	
	function [idList, Mu] = kmeansClustering(Data, nbStates)
		% Initialization of the model with k-means.
		% Copyright (c) 2016, Idiap Research Institute
		% Written by Sylvain Calinon, Tohid Alizadeh, 2013
		
		%Criterion to stop the EM iterative update
		cumdist_threshold = 1e-10;
		maxIter = 100;
		
		%Initialization of the parameters
		[nbVar, nbData] = size(Data);
		cumdist_old = -realmax;
		nbStep = 0;
		
		idTmp = randperm(nbData);
		Mu = Data(:,idTmp(1:nbStates));
		
		%k-means iterations
		while 1
			%E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			for i=1:nbStates
				%Compute distances
				distTmp(:,i) = sum((Data-repmat(Mu(:,i),1,nbData)).^2);
			end
			[vTmp,idList] = min(distTmp,[],2);
			cumdist = sum(vTmp);
			%M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			for i=1:nbStates
				%Update the centers
				Mu(:,i) = mean(Data(:,idList==i),2);
			end
			%Stopping criterion %%%%%%%%%%%%%%%%%%%%
			if abs(cumdist-cumdist_old) < cumdist_threshold
				break;
			end
			cumdist_old = cumdist;
			nbStep = nbStep+1;
			%   if nbStep>maxIter
			%     disp(['Maximum number of iterations, ' num2str(maxIter) 'is reached']);
			%     break;
			%   end
		end
	end

	end
	
end
