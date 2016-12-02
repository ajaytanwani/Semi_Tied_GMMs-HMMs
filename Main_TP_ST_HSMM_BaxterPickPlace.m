function Main_TP_ST_HSMM_BaxterPickPlace()

%%%% Task-Parameterized Semi-Tied Hidden Semi-Markov Model - TP_ST_HSMM.
%%%% The algorithm exploits the spatial and temporal correlation in the
%%%% demonstrations by sharing parameters in the output state distribution of an HSMM
%%%% in a semi-tied manner. The model parameters are adopted in accordance with the
%%%% changing environmental situations in an autonomous manner. Here, we
%%%% show its example to teach Baxter the task of picking an object from
%%%% different initial configurations, avoiding an obstacle of variable
%%%% height and placing it on a fixed target.
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
addpath(genpath('./algorithms'))
%% Load the dataset and the parameters

data_path = 'data/Baxter_PickPlace.mat';
load(data_path);

model.algo = 'TP_STHSMM';                        
model.controller = 'HSMM_LQR';              
model.data_gen = data_path;

model.alpha = 1; %%% SET model.alpha = 0 for task-parameterized hidden semi-Markov model 
                 %%% SET model.alpha = 1 for task-parameterized semi-tied hidden semi-Markov model 
				 %%% SET model.alpha = (0,1) for intermediate tying of basis vectors

model.nbD = 200;
model.dt = 0.1;
model.nbFrames = 2;
        
model.nbVarPos = 9; %Dimension of position data 
model.nbDeriv = 2; %Number of static & dynamic features (D=2 for [x,dx])
model.nbVar = model.nbVarPos * model.nbDeriv - 1; %Dimension of state vector

model.nbSamples = 4;
model.nbStates = 7; %Number of components in the GMM
model.time_dim = false;
model.hsmm_left_right_init = true;
model.rfactor = 7E0;
model.varnames = {'x', 'y', 'z', 'qw', 'qx', 'qy', 'qz','grip'};

Data = [];
for n=1:model.nbSamples    
    Data = [Data s(n).Data0];
end
DataRef = Data;
%% Regenerate the Dataset observed with respect to different frames

%Create 3rd order tensor data with XHAT instead of X
Data = zeros(model.nbVar, model.nbFrames, model.nbD);
for n=1:model.nbSamples
    DataTmp = [s(n).Data0; s(n).dData0]; %Add derivatives    
    for m=1:model.nbFrames
        Data(:,m,(n-1)*model.nbD+1:n*model.nbD) = s(n).p(m).A \ (DataTmp - repmat(s(n).p(m).b, 1, model.nbD));
    end
end

%% Initialize the parameters of the tensor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters initialization of tensor model with EM:');

model = initFuns.init_tensorGMM_timeBased(Data, model); %Initialization with time
model.Mu = model.Mu(2:end,:,:);
model.Sigma = model.Sigma(2:end,2:end,:,:);
[Data s model] = auxFuns.RemoveTimeDim(Data, s, model);

% VisualizeDataFrames(Data, s, model);

%% Tensor Model learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters estimation of tensor model with EM:');

model = EM_tensorSTGMM(Data, model);      %%Task Parameterized Semi-Tied Gaussian Mixture Model
model = EM_tensorSTHSMM(s, model);       %%Task Parameterized Semi-Tied Hidden Semi-Markov Model
%%% model = EM_tensorHSMM(s, model);       %%Task Parameterized Hidden Semi-Markov Model

%% Reproduction for the task parameters used to train the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = HSMM_LQR_controller(s, model);              %%State-based Trajectory Retrieval

MSE = 0;
for n=1:model.nbSamples
    AbsErr = abs(s(n).Data0 - r(n).Data).^2;
    MSE = MSE + sum(AbsErr(:))/numel(s(n).Data0);
end
model.MSE = MSE/model.nbSamples;

%% Reproduction for the Test Set

if ~isempty(stest)    
    rtest = HSMM_LQR_controller(stest, model);              %%State-based Trajectory Retrieval

    MSE = 0;
    for n=1:size(stest,2)
        AbsErr = abs(stest(n).Data0 - rtest(n).Data).^2;
        MSE = MSE + sum(AbsErr(:))/numel(stest(n).Data0);
    end
    model.MSE_Test = MSE/size(stest,2);
end
model
%% Reproduction for new environmnental situation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
disp('New reproductions with controller...');

for n=1:model.nbSamples
    for m=1:model.nbFrames
        %Random generation of new task parameters
        id=ceil(rand(2,1)*model.nbSamples);
        w=rand(2); w=w/sum(w);
        rnew(n).p(m).b = s(id(1)).p(m).b * w(1) + s(id(2)).p(m).b * w(2);
        rnew(n).p(m).A = s(id(1)).p(m).A * w(1) + s(id(2)).p(m).A * w(2);
        rnew(n).Data1 = s(n).Data1;
    end
end

r_new = HSMM_LQR_controller(rnew, model);          %%State-based Trajectory Retrieval


%% Analysis and draw Plots

DrawPlotsStateDriven(model, Data, [s stest], [r rtest], r_new);

end
