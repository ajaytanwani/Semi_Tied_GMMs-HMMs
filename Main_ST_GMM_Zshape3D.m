function Main_ST_GMM_Zshape3D()

%%%% Semi-Tied Gaussian Mixture Model (ST_GMM)- tie the covariance matrices of the 
%%%% Gaussian mixture model with common synergistic directions/basis vectors and 
%%%% component specific diagonal matrix for the Z-shaped 3D movement data. 
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
%% Select the Dataset and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('data/Zshape3D.mat');
model.nbVar = size(Data,1);

model.algo = 'STGMM'; 

model.nbStates = 3;
model.time_dim = false; 
model.nbSamples = 1;
model.nbD = 100;

%%% algorithm parameters
model.alpha = 1.0;	 %%% SET model.alpha = 0 for GMM
	                 %%% SET model.alpha = 1 for semi-tied GMM
					 %%% SET model.alpha = (0,1) for intermediate tying of basis vectors
model.B_sf = 5E-2;

model.col = [0.8 0 0; 0 0 0.8; 0 0.8 0; rand(20,3)];
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
           
%% Comparison with GMM 

model0 = modelRef;
[model0, LL0] = EM_GMM(Data, model0);

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
title('GMM')
set(h1, 'XTick', [minC, (minC+maxC)/2, maxC]);
L=cellfun(@(x)sprintf('%2.1f',x),num2cell(get(h1,'xtick')),'Un',0);
set(h1,'xticklabel',L)

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
%%`% figure for the Zshaped data

figure('color',[1 1 1],'Position',[400 400 800 400]);
xx = round(linspace(1,64,9));
xx2 = round(linspace(1,64,3));
clrmap = colormap('jet');
clrmap = min(clrmap(xx,:),.95);

clrmap2 = colormap('jet');
clrmap2 = min(clrmap2(xx2,:),.95);

for i=1:5
    plot3(Data(1,(i-1)*300 + 1:300*i), Data(2,(i-1)*300 + 1:300*i),Data(3,(i-1)*300 + 1:300*i),'-','linewidth',1.5,'color',[0.5625 0.5625 0.5625]);
    hold on
end
hold on

model.L = model.H;

%%% comparison with basis vectors obtained with PCA
% [coeff tmp lt1]= pca(Data');
% X = Data';
% meanX = mean(X);
% X = X - repmat(meanX,size(X,1),1);
% C = X'*X;
% [U D] = svd(C);
% lt = diag(D);
% model.L = U;

clrlist1 = [2,3,1];
for i=1:model.nbVar
    plotFuns.mArrow3(zeros(model.nbVar,1),model.L(:,i),'color',clrmap2(clrlist1(i),:),'stemWidth',0.75, 'tipWidth',1.0, 'facealpha',0.75);
end

clrlist = [2,7,4];
for i=1:model.nbStates
    plotFuns.plotGMM3D(model.Mu(:,i), model.Sigma(:,:,i), clrmap(clrlist(i),:), .5);
    
    for j=1:model.nbVar
        plotFuns.mArrow3(model.Mu(:,i), model.Mu(:,i) + model.H(:,j).*(model.SigmaDiag(j,j,i).^0.5),'color',clrmap2(clrlist1(j),:),'stemWidth',0.75, 'tipWidth',1.25, 'facealpha',1);
    end
end
axis off;
view(-40,6);
box off;
axis equal

i=1;text(model.Mu(1,i)+9,model.Mu(2,i) + 45,model.Mu(3,i),num2str(i),'FontSize',10);
i=2;text(model.Mu(1,i)+ 15,model.Mu(2,i),model.Mu(3,i),num2str(i),'FontSize',10);
i=3;text(model.Mu(1,i),model.Mu(2,i),model.Mu(3,i)+6,num2str(i),'FontSize',10);
hold on;

ax=get(gca,'Position');
ax(1) = ax(1) + 0.63;
ax(2) = ax(2) + 0.15;
ax(4)=ax(4)*0.35; 
ax(3)=ax(3)*0.35; 
axes('Position',ax);

imagesc(model.STMM_CorrMetric);
%     cmap = flipud(colormap('hot'));
colormap(cmap);
h = colorbar;
minC = min(min(model.STMM_CorrMetric)); maxC = max(max(model.STMM_CorrMetric));

set(h, 'XTick', [minC, (minC+maxC)/2, maxC]);

L=cellfun(@(x)sprintf('%2.1f',x),num2cell(get(h,'xtick')),'Un',0);
set(h,'xticklabel',L)

set(gca,'clim',[minC,maxC],'xtick',[1 2 3],'ytick',[1 2 3]);
axis square;

end
