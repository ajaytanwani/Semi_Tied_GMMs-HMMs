function DrawPlotsStateDriven(model, Data, s, r, r_new)
% This function contains plotting features for data analysis of task-parameterized hidden semi-Markov model
%
% Copyright (c) 2016 Idiap Research Institute
% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>

figure('Name', 'Demos and Repros','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

xx = round(linspace(1,64,model.nbVarPos));
clrmap = colormap('jet');
clrmap = min(clrmap(xx,:),.95);

MSE = 0;
nbSamples = model.nbSamples;
nbD = model.nbD;
for n=1:nbSamples
    AbsErr = abs(s(n).Data0 - r(n).Data).^2;
    MSE = MSE + sum(AbsErr(:))/numel(s(n).Data0);
    
    for i=1:model.nbVarPos
        subplot(2,nbSamples, n); hold on; plot(s(n).Data0(i,:)','color',clrmap(i,:));
        subplot(2,nbSamples,n + nbSamples); hold on; plot(r(n).Data(i,:)','color',clrmap(i,:));
    end
    title(num2str(n));
%     pause;
end
legend(gca,'x','y','z','qw','qx','qy','qz','grip');
suptitle('Demonstrations (above) and Reproductions (below)')

%%

figure('Name', 'HSMM Sequence 2','color',[1 1 1],'Position',[200 200 700 400]);
nbD = model.nbD;
nbSamples = model.nbSamples;
clrmap = lines(model.nbStates);

%Timeline plot of the state sequence probabilities
plotFuns.subaxis(1, 2, 2,'spacing',0.1); hold on;
for i=1:model.nbStates
	patch([1, 1:nbD, nbD], [0, r(1).alpha(i,:), 0], clrmap(i,:), 'EdgeColor', 'none', 'facealpha', .6);
	plot(1:nbD, r(1).alpha(i,:), 'linewidth', 2, 'color', clrmap(i,:));
end
set(gca,'xtick',[0:100:nbD],'ytick',[0 0.5 1],'fontsize',12); axis([1 nbD 0 1]);
xlabel('$t$','fontsize',16,'interpreter','latex'); 
ylabel('$h^{\mathrm{HSMM}}$','fontsize',16,'interpreter','latex');
axis square;

plotFuns.subaxis(1,2,1,'spacing',0);hold on;
plotFuns.plotHSMM(model.Trans, model.StatesPriors, model.Pd);
axis off
axis square;

%%
figure('Name', 'Frame Observation','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

xx = round(linspace(1,64,nbSamples));
clrmap = colormap('jet');
clrmap = min(clrmap(xx,:),.95);
clrmap2 = lines(model.nbStates);

for n=1:nbSamples
    for j=1:model.nbVarPos
        subplot(model.nbVarPos,model.nbFrames,model.nbFrames*(j-1)+1);hold on;plot(squeeze(Data(j,1,(n-1)*nbD + 1:n*nbD))','color',clrmap(n,:),'LineWidth',1.5);xlabel('t'); str = sprintf('x%d',j);ylabel(str);
        
        limAxes = [1, nbD, min(Data(j,1,:))-0.3, max(Data(j,1,:))+0.3];
        msh=[]; x0=[];
        for t=1:nbD-1
            if size(msh,2)==0
                msh(:,1) = [t; model.Mu(j,1,r(n).q(t))];
            end
            if t==nbD-1 || r(n).q(t+1)~=r(n).q(t)
                msh(:,2) = [t+1; model.Mu(j,1,r(n).q(t))];
                sTmp = model.Sigma(j,j,1, r(n).q(t))^.5;
                msh2 = [msh(:,1)+[0;sTmp], msh(:,2)+[0;sTmp], msh(:,2)-[0;sTmp], msh(:,1)-[0;sTmp], msh(:,1)+[0;sTmp]];
				
				patch(msh2(1,:), msh2(2,:), clrmap2(r(n).q(t),:),'edgecolor',clrmap2(r(n).q(t),:),'FaceAlpha',0.1);
				
                plot(msh(1,:), msh(2,:), '-','linewidth',1,'color',[.7 .7 .7]);
                plot([msh(1,1) msh(1,1)], limAxes(3:4),':','linewidth',1,'color',[.7 .7 .7]);
                x0 = [x0 msh];
                msh=[];
            end
        end
        axis(limAxes);
        
        
        subplot(model.nbVarPos,model.nbFrames,model.nbFrames*j);hold on; plot(squeeze(Data(j,2,(n-1)*nbD + 1:n*nbD))','color',clrmap(n,:),'LineWidth',1.5);xlabel('t'); str = sprintf('x%d',j);ylabel(str);
        limAxes = [1, nbD, min(Data(j,2,:))-0.3, max(Data(j,2,:))+0.3];
        msh=[]; x0=[];
        for t=1:nbD-1
            if size(msh,2)==0
                msh(:,1) = [t; model.Mu(j,2,r(n).q(t))];
            end
            if t==nbD-1 || r(n).q(t+1)~=r(n).q(t)
                msh(:,2) = [t+1; model.Mu(j,2,r(n).q(t))];
                sTmp = model.Sigma(j,j,2, r(n).q(t))^.5;
                msh2 = [msh(:,1)+[0;sTmp], msh(:,2)+[0;sTmp], msh(:,2)-[0;sTmp], msh(:,1)-[0;sTmp], msh(:,1)+[0;sTmp]];
				
				patch(msh2(1,:), msh2(2,:), clrmap2(r(n).q(t),:),'edgecolor',clrmap2(r(n).q(t),:),'FaceAlpha',0.1);
				                
                plot(msh(1,:), msh(2,:), '-','linewidth',1,'color',[.7 .7 .7]);
                plot([msh(1,1) msh(1,1)], limAxes(3:4),':','linewidth',1,'color',[.7 .7 .7]);
                x0 = [x0 msh];
                msh=[];
            end
        end
        axis(limAxes);
        
    end
    
end
suptitle('Frame 1 (left) vs Frame 2 (right) Analysis')

%%
figure('Name','Variance in Demos with Frames','color',[1 1 1],'Position',[400 400 1650 450]);
nbSamples = model.nbSamples;
nbD = model.nbD;
xx = round(linspace(1,64,nbSamples));
clrmap = colormap('jet');
clrmap = min(clrmap(xx,:),.95);

clrmap2 = lines(model.nbStates);
str1 = {'x', 'y', 'z', 'q_{w}', 'q_{x}', 'q_{y}', 'q_{z}','grip'};
for n=1:nbSamples
    for j=1:model.nbVarPos
        
        plotFuns.subaxis(model.nbFrames,model.nbVarPos,j,'SpacingHoriz',0.0,'SpacingVert',0.075);hold on;
        p1 = plot(squeeze(Data(j,1,(n-1)*nbD + 1:n*nbD))','color',clrmap(n,:),'LineWidth',1.5); str = sprintf('%s^{(%d)}',str1{j},1);title(str);
        
        a1 = round(min(Data(j,1,:))-0.2,1);        
        a2 = round(max(Data(j,1,:))+0.2,1);
       
        a1 = -1.2; a2 = 1.2;            
        
        a3 = round(mean([a1 a2]),1);
        
        if j == 1
            set(gca, 'xTick',[],'yTick', [a1, a3, a2]);
        else
            set(gca, 'xTick',[],'yTick', []);
        end
%         set(gca, 'xTick',[],'yTick', [-1, 0, 1]);
%         L=cellfun(@(x)sprintf('%2.1f',x),num2cell(get(gca,'xtick')),'Un',0);
%         set(gca,'xticklabel',L)
        
        limAxes = [1, nbD, a1, a2];
        msh=[]; x0=[];
        for t=1:nbD-1
            if size(msh,2)==0
                msh(:,1) = [t; model.Mu(j,1,r(n).q(t))];
            end
            if t==nbD-1 || r(n).q(t+1)~=r(n).q(t)
                msh(:,2) = [t+1; model.Mu(j,1,r(n).q(t))];
                sTmp = model.Sigma(j,j,1, r(n).q(t))^.5;
                msh2 = [msh(:,1)+[0;sTmp], msh(:,2)+[0;sTmp], msh(:,2)-[0;sTmp], msh(:,1)-[0;sTmp], msh(:,1)+[0;sTmp]];
				
                patch(msh2(1,:), msh2(2,:), clrmap2(r(n).q(t),:),'edgecolor',clrmap2(r(n).q(t),:),'FaceAlpha',0.1);
				
                if r(n).q(t) == 3 || r(n).q(t) == 4
                    plot([msh(1,1) msh(1,1)], limAxes(3:4),'-.','linewidth',1,'color',clrmap(1,:));
                else
                    plot([msh(1,1) msh(1,1)], limAxes(3:4),':','linewidth',1,'color',[.7 .7 .7]);
                end
                x0 = [x0 msh];
                msh=[];
            end
        end
        axis(limAxes);
        axis square;
        
        plotFuns.subaxis(model.nbFrames,model.nbVarPos,j+model.nbVarPos,'SpacingHoriz',0.0,'SpacingVert',0.075);hold on; plot(squeeze(Data(j,2,(n-1)*nbD + 1:n*nbD))','color',clrmap(n,:),'LineWidth',1.5);
        xlabel('$t$','fontsize',16,'interpreter','latex'); 
        str = sprintf('%s^{(%d)}',str1{j},2);title(str);
        
        a1 = round(min(Data(j,2,:))-0.2,1);
        a2 = round(max(Data(j,2,:))+0.2,1);
		a1 = -1.2; a2 = 1.2;
        a3 = round(mean([a1 a2]),1);
        if j== 1
            set(gca, 'xTick',[1 100 200],'yTick', [a1, a3, a2]);
        else
            set(gca, 'xTick',[1 100 200],'yTick', []);
        end 
        limAxes = [1, nbD, a1, a2];
        msh=[]; x0=[];
        for t=1:nbD-1
            if size(msh,2)==0
                msh(:,1) = [t; model.Mu(j,2,r(n).q(t))];
            end
            if t==nbD-1 || r(n).q(t+1)~=r(n).q(t)
                msh(:,2) = [t+1; model.Mu(j,2,r(n).q(t))];
                sTmp = model.Sigma(j,j,2, r(n).q(t))^.5;
                msh2 = [msh(:,1)+[0;sTmp], msh(:,2)+[0;sTmp], msh(:,2)-[0;sTmp], msh(:,1)-[0;sTmp], msh(:,1)+[0;sTmp]];
				
				patch(msh2(1,:), msh2(2,:), clrmap2(r(n).q(t),:),'edgecolor',clrmap2(r(n).q(t),:),'FaceAlpha',0.1);
                if r(n).q(t) == 5 || r(n).q(t) == 6
                    plot([msh(1,1) msh(1,1)], limAxes(3:4),'-.','linewidth',1,'color',clrmap(1,:));
                else
                    plot([msh(1,1) msh(1,1)], limAxes(3:4),':','linewidth',1,'color',[.7 .7 .7]);
                end                
                x0 = [x0 msh];
                msh=[];
            end
        end
        axis(limAxes);
        axis square;

    end

end

%%
  
switch model.data_gen
	case 'data/Baxter_PickPlace.mat'
	%%
	figure('Name','Trajectory Reproductions','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
	% figure('color',[1 1 1]);
	nbD = model.nbD;
	nbSamples = model.nbSamples;
	% xx3 = round(linspace(1,64,model.nbStates));
	% clrmap3 = colormap('jet');
	% clrmap3 = min(clrmap3(xx3,:),.95);
	clrmap3 = lines(model.nbStates);
	
	xx2 = round(linspace(1,64,12));
	clrmap2 = colormap('jet');
	clrmap2 = min(clrmap2(xx2,:),.95);
	
	xx = round(linspace(1,64,3));
	clrmap = colormap('jet');
	clrmap = min(clrmap(xx,:),.95);
	
	clrlist = [4;5;9];
	framecol = reshape(clrmap2(clrlist,:),[1 3 3]);
	
	clrlist2 = [5;5;12];
	framecol3 = reshape(clrmap2(clrlist2,:),[1 3 3]);
	
	col2 = [0 0 0.8;0 0.8 0; 0.8 0 0];
	framecol2 = reshape(col2,[1 3 3]);
	
	object_width = 0.15;
	object_length = 0.1;
	
	X = [-object_width/2 object_width/2 object_width/2 -object_width/2 -object_width/2];
	Y = [object_length/2 object_length/2 -object_length/2 -object_length/2 object_length/2];
	P = [X;Y;zeros(1,5)];
	% 		F = [1, 2, 3, 4];
	for i=1:size(r,2)
		plotFuns.subaxis(2,4,i,'spacing',0);
		plot3(r(i).Data(1,:), r(i).Data(2,:), r(i).Data(3,:), '-','linewidth',2,'Color',[0.5 0.5 0.5]);% title(num2str(i));
		hold on; plot3(r(i).Data(1,1), r(i).Data(2,1), r(i).Data(3,1), 's','Color', clrmap(2,:),'MarkerSize',10,'MarkerFaceColor',clrmap(2,:));
		hold on; plot3(r(i).Data(1,end), r(i).Data(2,end), r(i).Data(3,end), 'x','Color', clrmap(1,:),'MarkerSize',10,'LineWidth',5,'MarkerFaceColor',clrmap(1,:));
		hold on;
		
		hold on; plotFuns.plot3Dframe2(s(i).p(2).A(1:3,1:3),s(i).p(2).b(1:3),0.1, squeeze(framecol(1,:,:))','-',0.005,0.008,[1,3]);
		hold on; ha = plotFuns.plot3Dframe2(s(i).p(1).A(1:3,1:3),s(i).p(1).b(1:3),0.1,squeeze(framecol2(1,:,:))','-',0.005,0.008,[1,3]);
		
		
		orgPt1 = s(i).p(2).b(1:3);
		endPt = 0.25*[-1;0;0];
		plotFuns.mArrow3(orgPt1-endPt,orgPt1+endPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.005, 'facealpha',0.2);hold on;
		hold on;
		% 			obj_cords = s(i).p(1).A(1:3,1:3)*[P(1,:) + s(i).p(1).b(1); P(2,:) + s(i).p(1).b(2); P(3,:) + s(i).p(1).b(3)];
		obj_cords = s(i).p(1).A(1:3,1:3)*P + repmat(s(i).p(1).b(1:3), 1, size(P,2));%[P(1,:) + s(i).p(1).b(1); P(2,:) + s(i).p(1).b(2); P(3,:) + s(i).p(1).b(3)];
		
		patch(obj_cords(1,:), obj_cords(2,:),obj_cords(3,:),'EdgeColor','green','FaceColor','none','LineWidth',2);
		% 			[cent_X centr_Y] = s(i).p(1).b(1:2);
		%     plotFuns.mArrow3(orgPt1,orgPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.001, 'facealpha',0.2);hold on;
		
		for j=1:10:nbD
			Orient_Frame = Quaternion(r(i).Data(4:7,j));
			hold on;h3 = plotFuns.plot3Dframe2(Orient_Frame.R,r(i).Data(1:3,j),0.06,squeeze(framecol3(1,:,:))','-', 0.0025,0.004,[1,3]);
		end
		
		text(s(i).p(1).b(1),s(i).p(1).b(2) + 0.075,s(i).p(1).b(3)+0.03, '$$  \bf{\textit{A}_{1}} $$','Interpreter','latex', 'FontSize', 20)
		text(s(i).p(2).b(1),s(i).p(2).b(2) + 0.075,s(i).p(2).b(3)+0.025, '$$  \bf{\textit{A}_{2}} $$','Interpreter','latex', 'FontSize', 20)
		
		for j=1:model.nbStates
			hold on;plotFuns.plotGMM3D(r(i).Mu(1:3,j),r(i).Sigma(1:3,1:3,j),clrmap3(j,:),0.5);
		end
		grid on;box on
		set(gca,'xtick',[],'ytick',[],'ztick',[]);
		
		% 			view(-161,34);
		view(66.5,30);
		axis equal;
		av1 = axis;
	end
	case 'data/BaxterValveDemos2.mat'
	%%
	nbD = model.nbD; nbSamples = model.nbSamples;
	figure('Name','Trajectory Reproductions','color',[1 1 1],'Position',[400 400 1000 600]);
	nbD = model.nbD;

	clrmap3 = lines(model.nbStates);

	xx2 = round(linspace(1,64,12));
	clrmap2 = colormap('jet');
	clrmap2 = min(clrmap2(xx2,:),.95);

	xx = round(linspace(1,64,3));
	clrmap = colormap('jet');
	clrmap = min(clrmap(xx,:),.95);

	clrlist = [4;5;9];
	framecol = reshape(clrmap2(clrlist,:),[1 3 3]);

	clrlist2 = [5;5;12];
	framecol3 = reshape(clrmap2(clrlist2,:),[1 3 3]);

	col2 = [0 0 0.8;0 0.8 0; 0.8 0 0];
	framecol2 = reshape(col2,[1 3 3]);

	plotFuns.subaxis(1,2,1,'spacing',0);
	i = 1;
	plot3(r(i).Data(1,:), r(i).Data(2,:), r(i).Data(3,:), '-','linewidth',2,'Color',[0.5 0.5 0.5]);
	hold on; plot3(r(i).Data(1,1), r(i).Data(2,1), r(i).Data(3,1), 's','Color', clrmap(2,:),'MarkerSize',10,'MarkerFaceColor',clrmap(2,:));
	hold on; plot3(r(i).Data(1,end), r(i).Data(2,end), r(i).Data(3,end), 'x','Color', clrmap(1,:),'MarkerSize',10,'LineWidth',5,'MarkerFaceColor',clrmap(1,:));
	hold on;

	hold on; plotFuns.plot3Dframe2(s(1).p(2).A(1:3,1:3),s(1).p(2).b(1:3),0.1, squeeze(framecol(1,:,:))','-',0.005,0.008,[1,3]);

	for i=1:1
		hold on; ha = plotFuns.plot3Dframe2(s(i).p(1).A(1:3,1:3),s(i).p(1).b(1:3),0.1,squeeze(framecol2(1,:,:))','-',0.005,0.008,[1,3]);
	end
	orgPt1 = mean([s(1).p(1).b(1:3), s(1).p(2).b(1:3)],2);
	endPt = s(1).p(2).A(1:3,1:3)*0.1*[0;0;1];
	plotFuns.mArrow3(orgPt1,orgPt1+endPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.005, 'facealpha',0.2);hold on;
	orgPt = s(1).p(1).b(1:3);
	plotFuns.mArrow3(orgPt1,orgPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.001, 'facealpha',0.2);hold on;

	for j=1:10:nbD
		Orient_Frame = Quaternion(r(i).Data(4:7,j));
		hold on;h3 = plotFuns.plot3Dframe2(Orient_Frame.R,r(i).Data(1:3,j),0.06,squeeze(framecol3(1,:,:))','-', 0.0025,0.004,[1,3]);
	end

	text(s(1).p(1).b(1) + 0.05,s(1).p(1).b(2) - 0.02,s(1).p(1).b(3)+0.03, '$$  \bf{\textit{A}_{1}} $$','Interpreter','latex', 'FontSize', 20)
	text(s(1).p(2).b(1) + 0.09,s(1).p(2).b(2) - 0.02,s(1).p(2).b(3)+0.01, '$$  \bf{\textit{A}_{2}} $$','Interpreter','latex', 'FontSize', 20)

	for i=1:model.nbStates
		hold on;plotFuns.plotGMM3D(r(1).Mu(1:3,i),r(1).Sigma(1:3,1:3,i),clrmap3(i,:),0.5);
	end

	% axis([0.01 1.5 0.01 1 -0.4 -0.05]);
	grid on;box on
	set(gca,'xtick',[],'ytick',[],'ztick',[]);
	view(-161,34);
	axis equal;
	av1 = axis;

	plotFuns.subaxis(1,2,2,'spacing',0.01);
	i = 1;
	plot3(r_new(i).Data(1,:), r_new(i).Data(2,:), r_new(i).Data(3,:), '-','linewidth',2,'Color',[0.5 0.5 0.5]);
	hold on; plot3(r_new(i).Data(1,1), r_new(i).Data(2,1), r_new(i).Data(3,1), 's','Color', clrmap(2,:),'MarkerSize',10,'MarkerFaceColor',clrmap(2,:));
	hold on; plot3(r_new(i).Data(1,end), r_new(i).Data(2,end), r_new(i).Data(3,end), 'x','Color', clrmap(1,:),'MarkerSize',10,'LineWidth',5,'MarkerFaceColor',clrmap(1,:));
	%         hold on; plotFuns.plotGMM3D(r_new(i).Data(1:3,1:10:end), r_new(i).expSigma(1:3,1:3,1:10:end)*0.5, clrmap(2,:), .15, 2);
	hold on; plotFuns.plot3Dframe2(r_new(1).p(2).A(1:3,1:3),r_new(1).p(2).b(1:3),0.1,squeeze(framecol(1,:,:))','-',0.005,0.008,[1,3]);

	for i=1:1
		hold on; ha = plotFuns.plot3Dframe2(r_new(i).p(1).A(1:3,1:3),r_new(i).p(1).b(1:3),0.1,squeeze(framecol2(1,:,:))','-',0.005,0.008,[1,3]);
	end

	orgPt1 = mean([r_new(1).p(1).b(1:3), r_new(1).p(2).b(1:3)],2);
	endPt = r_new(1).p(2).A(1:3,1:3)*0.1*[0;0;1];
	plotFuns.mArrow3(orgPt1,orgPt1+endPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.005, 'facealpha',0.2);hold on;
	orgPt = r_new(i).p(1).b(1:3);
	plotFuns.mArrow3(orgPt1,orgPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.001, 'facealpha',0.2);hold on;

	for j=1:10:nbD
		Orient_Frame = Quaternion(r_new(i).Data(4:model.nbVarPos,j));
		hold on;h3 = plotFuns.plot3Dframe2(Orient_Frame.R,r_new(i).Data(1:3,j),0.06,squeeze(framecol3(1,:,:))','--',0.0025,0.004,[1,3]);
	end

	for i=1:model.nbStates
		plotFuns.plotGMM3D(r_new(1).Mu(1:3,i),r_new(1).Sigma(1:3,1:3,i),clrmap3(i,:),0.5);
	end

	text(r_new(1).p(1).b(1) - 0.03,r_new(1).p(1).b(2) - 0.02,r_new(1).p(1).b(3)+0.03, '$$   \bf{\tilde{\textit{A}}}_1 $$','Interpreter','latex', 'FontSize', 20)
	text(r_new(1).p(2).b(1) + 0.09,r_new(1).p(2).b(2) - 0.02,r_new(1).p(2).b(3)+0.01, '$$  \bf{\tilde{\textit{A}}}_2 $$','Interpreter','latex', 'FontSize', 20)

	view(-161,34);
	% plot3(r_new(i).Data(2,end), r_new(i).Data(3,end), r_new(i).Data(4,end), 'o','Color', clrmap(3,:),'MarkerSize',10,'MarkerFaceColor',clrmap(3,:));
	grid on;box on
	set(gca,'xtick',[],'ytick',[],'ztick',[]);
	axis equal;
	av2 = axis;
	%         axis([av1(1) av1(2) av1(3) av1(4) av1(5) av1(6)])

	ax2 = get(gca,'position');
	ax2(1) = ax2(1) - 0.15;
	ax2(2) = ax2(2) + 0.42;
	ax2(3) = ax2(3)*0.7;
	ax2(4) = ax2(4)*0.7;
	a2 = axes('position',ax2);

	plotFuns.plot3Dframe2(s(1).p(2).A(1:3,1:3),s(1).p(2).b(1:3),0.1, squeeze(framecol(1,:,:))','-',0.005,0.009,[1,3]);
	hold on;
	for i=1:nbSamples/2
		ha1 = plotFuns.plot3Dframe2(s(i).p(1).A(1:3,1:3),s(i).p(1).b(1:3),0.1,squeeze(framecol2(1,:,:))','-',0.005,0.009,[1,3]);hold on;
	end
	view(-161,34);
	box on;
	set(gca,'xtick',[],'ytick',[],'ztick',[]);
	text(s(1).p(1).b(1) - 0.04,s(1).p(1).b(2) - 0.02,s(1).p(1).b(3)+0.02, '$$ \bf{\textit{A}_{1}^{(1)}} $$','Interpreter','latex', 'FontSize', 18)
	text(s(2).p(1).b(1)-0.04,s(2).p(1).b(2) - 0.02,s(2).p(1).b(3)-0.025, '$$ \bf{\textit{A}_{1}^{(2)}} $$','Interpreter','latex', 'FontSize', 18)
	text(s(3).p(1).b(1)- 0.01,s(3).p(1).b(2) - 0.01,s(3).p(1).b(3)-0.05, '$$ \bf{\textit{A}_{1}^{(3)}} $$','Interpreter','latex', 'FontSize', 18)
	text(s(4).p(1).b(1) + 0.01,s(4).p(1).b(2) + 0.015,s(4).p(1).b(3)-0.025, '$$ \bf{\textit{A}_{1}^{(4)}} $$','Interpreter','latex', 'FontSize', 18)
	text(s(1).p(2).b(1) + 0.06,s(1).p(2).b(2) + 0.015,s(1).p(2).b(3)-0.06, '$$ \bf{\textit{A}_{2}^{(1 \ldots 4)}} $$','Interpreter','latex', 'FontSize', 18)

	orgPt = mean([s(1).p(1).b(1:3), s(1).p(2).b(1:3)],2);
	endPt = s(1).p(2).A(1:3,1:3)*0.1*[0;0;1];
	plotFuns.mArrow3(orgPt,orgPt+endPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.005, 'facealpha',0.2);hold on;
	endPt = s(1).p(2).b(1:3);
	plotFuns.mArrow3(orgPt,endPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.0005, 'facealpha',0.2);hold on;
	view(-161,34);
	axis equal;
	a2v = axis;	
end
