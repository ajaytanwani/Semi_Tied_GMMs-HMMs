function Main_ST_GMM_ChickenDance()

%%%% Semi-Tied Gaussian Mixture Model (ST_GMM)- tie the covariance matrices of the 
%%%% Gaussian mixture model with common synergistic directions/basis vectors and 
%%%% component specific diagonal matrix for the chicken dance movement. 
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
% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>

clc
close all

addpath('./utils/');
addpath('./controllers')
addpath(genpath('./algorithms'))
%% Select the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.algo = 'STGMM'; 
model.controller = 'GMR';

model.nbStates = 75;
model.time_dim = true; 
model.nbSamples = 1;
model.nbD = 500;

%%% algorithm parameters
model.alpha = 1; %% SET model.alpha = 0 for GMM
                 %%% SET model.alpha = 1 for semi-tied GMM
				 %%% SET model.alpha = (0,1) for intermediate tying of basis vectors

model.diagRegularizationFactor = 1E-8;
model.nbMaxSteps = 10;
model.nbVariationSteps = 10;

model.col = [0.8 0 0; 0 0 0.8; 0 0.8 0; rand(20,3)];

%% Load the Dataset

load('data/ChickenDance.mat');

%Resampling
Data = []; %%%% size is njointsx3
for m=1:mot(2).njoints
    Data = [Data; mot(2).jointTrajectories{m}];
end
Data = spline(1:size(Data,2), Data , linspace(1,size(Data,2),model.nbD));
Data = [[1:model.nbD]; Data];

model.nbVar = size(Data,1);


%% Initialize Mu, and Sigma for GMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = initFuns.init_GMM_timeBased(Data, model);

if isfield(model,'time_dim')
    if ~model.time_dim
        model.Mu = model.Mu(2:model.nbVar,:); 
        model.Sigma = model.Sigma(2:model.nbVar,2:model.nbVar,:); 
        Data = Data(2:model.nbVar,:);
        model.nbVar= model.nbVar - 1;
    end
end

DataRef = Data;
modelRef= model;
%% Learn the probablistic model with: EM_STGMM

[model, LL] = EM_STGMM(Data,model);

%% Decode the trajectory model with: GMR_controller

r = controllerFuns.GMR_controller(Data,model);
           
%% Comparison with GMM 

model0 = modelRef;
[model0, LL0] = EM_GMM(Data, model0);
r0 = controllerFuns.GMR_controller(Data,model0);

%% Results 

%%% Correlation Analysis of Semi-Tied Gaussian Mixture Model and Gaussian Mixture Model

for i=1:model.nbStates
    STMM_CorrMat(:,i) = reshape(model.Sigma(:,:,i),1,size(Data,1)*size(Data,1));
    GMM_CorrMat(:,i) = reshape(model0.Sigma(:,:,i),1,size(Data,1)*size(Data,1));
end
model.STMM_CorrMetric = corr(STMM_CorrMat);
model.GMM_CorrMetric = corr(GMM_CorrMat);

%%%% figure for the covariance matrices correlation
figure('color',[1 1 1],'position',[400 400 900 300]);axis square, box off,
cmap = flipud(colormap('hot'));

minC = min(min(model.GMM_CorrMetric)); maxC = max(max(model.GMM_CorrMetric));
plotFuns.subaxis(1,2,1,'spacing',0);
imagesc(model.GMM_CorrMetric);colormap(cmap);
h1 = colorbar;
%set(gca,'clim',[minC,maxC]);

set(h1, 'XTick', [minC, (minC+maxC)/2, maxC]);
L=cellfun(@(x)sprintf('%2.1f',x),num2cell(get(h1,'xtick')),'Un',0);
set(h1,'xticklabel',L)
title('GMM')
set(gca,'clim',[minC,maxC],'xtick',[20 40 60],'ytick',[20 40 60]);
axis square;
minC = min(min(model.STMM_CorrMetric)); maxC = max(max(model.STMM_CorrMetric));
plotFuns.subaxis(1,2,2,'spacing',0);
imagesc(model.STMM_CorrMetric);colormap(cmap);
h2 = colorbar;

set(h2, 'XTick', [minC, (minC+maxC)/2, maxC]);
L=cellfun(@(x)sprintf('%2.2f',x),num2cell(get(h2,'xtick')),'Un',0);
set(h2,'xticklabel',L)

set(gca,'clim',[minC,maxC],'xtick',[20 40 60],'ytick',[20 40 60]);

axis square;
title('Semi-Tied GMM')
%%
%%`% figure for the motion sequence
skelTmp = skel;
motTmp = mot;

skelTmp(3) = skelTmp(2);
motTmp(3) = motTmp(2);
Data1 = [];
Data2 = [];
% model.mot = model.mot2;
for m=1:mot(2).njoints
    Data1 = [Data1; mot(1).jointTrajectories{m}];
    Data2 = [Data2; mot(2).jointTrajectories{m}];
end
Data1 = spline(1:size(Data1,2), Data1 , linspace(1,size(Data1,2),size(Data,2)));
Data2 = spline(1:size(Data2,2), Data2 , linspace(1,size(Data2,2),size(Data,2)));

plotFuns.DrawModel(skelTmp, motTmp, Data1, Data2, r.Data);

end
