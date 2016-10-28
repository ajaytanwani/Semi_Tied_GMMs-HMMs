function VisualizeDataFrames(Data, s, model)
% This function contains data analysis features of the datasets attached
%
% Copyright (c) 2016 Idiap Research Institute
% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>

nbSamples = model.nbSamples;
nbD = model.nbD;


%%
figure('Name', 'Orientation Comparison','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

xx = round(linspace(1,64,4));
clrmap2 = colormap('jet');
clrmap2 = min(clrmap2(xx,:),.95);

for n=1:model.nbSamples + 1
    for i=1:4
        if n == model.nbSamples + 1            
            subplot(2,nbSamples+1,n+nbSamples+1);hold on;plot(s(1).Frame2.rorient(:,i),'color',clrmap2(i,:));title('Orient Frame 2');        
        else
            subplot(2,nbSamples+1,n);hold on;plot(s(n).Data0(3+i,:),'color',clrmap2(i,:));str = sprintf('%s (%d)','Orient Data' ,n);title(str); 
            subplot(2,nbSamples+1,n+nbSamples+1);hold on;plot(s(n).Frame.rorient(:,i),'color',clrmap2(i,:));str = sprintf('%s (%d)','Orient Frame' ,n);title(str); 
        end
    end    
end

suptitle('Orientation of Data and Frame Comparison')
%%

if isfield(s(1),'dDataVis')
    nbSamples = model.nbSamples; nbD = model.nbD;
    figure('Name', 'rxDot Comparison','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

    xx = round(linspace(1,64,3));
    clrmap = colormap('jet');
    clrmap = min(clrmap(xx,:),.95);
    for n = 1:nbSamples %%%%%compare rdpos
        for i=1:3
            subplot(2,nbSamples,n); hold on; plot(s(n).dDataVis(i,:),'color', clrmap(i,:));str = sprintf('%s (%d)','Baxter xdot' ,n);title(str);     
            subplot(2,nbSamples,n+nbSamples); hold on; plot(s(n).dData0(i,:),'color', clrmap(i,:));str = sprintf('%s (%d)','Estimated xdot' ,n);title(str); 
        end
    end   
    suptitle('rxDot Comparison - Measured (above) vs Estimated (below)')
    
    figure('Name', 'rQuaternionDot Comparison','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

    xx = round(linspace(1,64,4));
    clrmap2 = colormap('jet');
    clrmap2 = min(clrmap2(xx,:),.95);
    for n = 1:nbSamples %%%%%compare rquat_dot
        for i=1:3
            subplot(2,nbSamples,n); hold on; plot(s(n).dDataVis(3+i,:),'color', clrmap2(i,:));str = sprintf('%s (%d)','QuatDot \omega' ,n);title(str);    
            subplot(2,nbSamples,n+nbSamples); hold on; plot(s(n).dData0(3+i,:),'color', clrmap2(i,:));str = sprintf('%s (%d)','QuatDot FD' ,n);title(str);     
        end
    end
    suptitle('rQuaternionDot Comparison - Raw Estimation (above) vs FiniteDiff (below)')
end    

%%

if isfield(s(1).DataRaw,'rdorient_hat')
    figure('Name', '\omega Comparison','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
    nbSamples = model.nbSamples; nbD = model.nbD;

    xx = round(linspace(1,64,3));
    clrmap = colormap('jet');
    clrmap = min(clrmap(xx,:),.95);
    for n = 1:nbSamples %%%%%compare romega vs romega_hat
        for i=1:3
            subplot(2,nbSamples,n); hold on; plot(s(n).DataRaw.rdorient(:,i),'color', clrmap(i,:));str = sprintf('%s (%d)','Baxter \omega' ,n);title(str);     
            subplot(2,nbSamples,n+nbSamples); hold on; plot(s(n).DataRaw.rdorient_hat(:,i),'color', clrmap(i,:));str = sprintf('%s (%d)','Estimated \omega' ,n);title(str);     
        end
    end
    suptitle('\omega Comparison - Measured (above) vs Estimated from quat (below)')
end
%%
figure('Name', 'Frame Observation','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

xx = round(linspace(1,64,nbSamples));
clrmap = colormap('jet');
clrmap = min(clrmap(xx,:),.95);

for n=1:nbSamples
    for j=1:model.nbVarPos
        subplot(model.nbVarPos,model.nbFrames,model.nbFrames*(j-1)+1);hold on;plot(squeeze(Data(j,1,(n-1)*nbD + 1:n*nbD))','color',clrmap(n,:),'LineWidth',1.5);xlabel('t'); str = sprintf('x%d',j);ylabel(str);       
        
        limAxes = [1, nbD, min(Data(j,1,:))-0.3, max(Data(j,1,:))+0.3];
        msh=[]; x0=[];
        axis(limAxes);
        
        
        subplot(model.nbVarPos,model.nbFrames,model.nbFrames*j);hold on; plot(squeeze(Data(j,2,(n-1)*nbD + 1:n*nbD))','color',clrmap(n,:),'LineWidth',1.5);xlabel('t'); str = sprintf('x%d',j);ylabel(str);
        limAxes = [1, nbD, min(Data(j,2,:))-0.3, max(Data(j,2,:))+0.3];
        msh=[]; x0=[];
        axis(limAxes);
    end
end

suptitle('Frame 1 (left) vs Frame 2 (right) Analysis')

%%

figure('Name', 'Compact Frame Observation','color',[1 1 1],'Position',[400 400 1650 450]);

xx = round(linspace(1,64,nbSamples));
clrmap = colormap('jet');
clrmap = min(clrmap(xx,:),.95);

str1 = {'x', 'y', 'z', 'q_{w}', 'q_{x}', 'q_{y}', 'q_{z}'};
for n=1:nbSamples
    for j=1:model.nbVarPos
        
        subaxis(model.nbFrames,model.nbVarPos,j,'SpacingHoriz',0.0,'SpacingVert',0.075);hold on;
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
        axis(limAxes);
        axis square;
                
        subaxis(model.nbFrames,model.nbVarPos,j+model.nbVarPos,'SpacingHoriz',0.0,'SpacingVert',0.075);hold on; plot(squeeze(Data(j,2,(n-1)*nbD + 1:n*nbD))','color',clrmap(n,:),'LineWidth',1.5);
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
%         set(gca, 'xTick',[1 100 200],'yTick', [-1, 0, 1]);
        limAxes = [1, nbD, a1, a2];
        axis(limAxes);
        axis square;

    end

end

%%
figure('Name', 'Position Orientaion with Frames','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

xx = round(linspace(1,64,3));
clrmap = colormap('jet');
clrmap = min(clrmap(xx,:),.95);

xx = round(linspace(1,64,4));
clrmap2 = colormap('jet');
clrmap2 = min(clrmap2(xx,:),.95);


for n=1:nbSamples
    for i=1:3
        subplot(2,2,1);hold on;plot(squeeze(Data(i,1,(n-1)*nbD + 1:n*nbD))','color',clrmap(i,:));
        subplot(2,2,2);hold on; plot(squeeze(Data(i,2,(n-1)*nbD + 1:n*nbD))','color',clrmap(i,:));
    end

    for i=1:4
        subplot(2,2,3);hold on;plot(squeeze(Data(3+i,1,(n-1)*nbD + 1:n*nbD))','color',clrmap2(i,:));
        subplot(2,2,4);hold on; plot(squeeze(Data(3+i,2,(n-1)*nbD + 1:n*nbD))','color',clrmap2(i,:));    
%         pause;
    end
end
subplot(2,2,1);hold on; title('Position Frame 1')
subplot(2,2,2);hold on; title('Position Frame 2')
subplot(2,2,3);hold on; title('Orientation Frame 1')
subplot(2,2,4);hold on; title('Orientation Frame 2')
suptitle('Position and Orientation Observaion: Frame 1 (left) vs Frame 2 (right)')

%%

figure('Name','Demonstrations','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
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

for i=1:nbSamples
    subaxis(2,4,i,'spacing',0);
    plot3(s(i).Data0(1,:), s(i).Data0(2,:), s(i).Data0(3,:), '-','linewidth',2,'Color',[0.5 0.5 0.5]); title(num2str(i));
    hold on; plot3(s(i).Data0(1,1), s(i).Data0(2,1), s(i).Data0(3,1), 's','Color', clrmap(2,:),'MarkerSize',10,'MarkerFaceColor',clrmap(2,:));
    hold on; plot3(s(i).Data0(1,end), s(i).Data0(2,end), s(i).Data0(3,end), 'x','Color', clrmap(1,:),'MarkerSize',10,'LineWidth',5,'MarkerFaceColor',clrmap(1,:));
    hold on;
    
    hold on; plot3Dframe2(s(i).p(2).A(1:3,1:3),s(i).p(2).b(1:3),0.1, squeeze(framecol(1,:,:))','-',0.005,0.008,[1,3]);
    hold on; ha = plot3Dframe2(s(i).p(1).A(1:3,1:3),s(i).p(1).b(1:3),0.1,squeeze(framecol2(1,:,:))','-',0.005,0.008,[1,3]);
    
    
    orgPt1 = mean([s(1).p(1).b(1:3), s(1).p(2).b(1:3)],2);
    endPt = s(1).p(2).A(1:3,1:3)*0.1*[0;0;1];
    mArrow3(orgPt1,orgPt1+endPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.005, 'facealpha',0.2);hold on;
    orgPt = s(i).p(1).b(1:3);
    mArrow3(orgPt1,orgPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.001, 'facealpha',0.2);hold on;
    
    for j=1:10:nbD
        Orient_Frame = Quaternion(s(i).Data0(4:model.nbVarPos,j));
        hold on;h3 = plot3Dframe2(Orient_Frame.R,s(i).Data0(1:3,j),0.06,squeeze(framecol3(1,:,:))','-', 0.0025,0.004,[1,3]);
    end
    
    text(s(i).p(1).b(1),s(i).p(1).b(2) + 0.075,s(i).p(1).b(3)+0.03, '$$  \bf{\textit{A}_{1}} $$','Interpreter','latex', 'FontSize', 20)
    text(s(i).p(2).b(1),s(i).p(2).b(2) + 0.075,s(i).p(2).b(3)+0.025, '$$  \bf{\textit{A}_{2}} $$','Interpreter','latex', 'FontSize', 20)
    
    grid on;box on
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
    
    view(-161,34);
    
    axis equal;
    av1 = axis;
end

                
        
%%

figure('Name','Frames','color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

plot3Dframe2(s(1).p(2).A(1:3,1:3),s(1).p(2).b(1:3),0.1, squeeze(framecol(1,:,:))','-',0.005,0.009,[1:3]);
hold on;
for i=1:nbSamples
    ha1 = plot3Dframe2(s(i).p(1).A(1:3,1:3),s(i).p(1).b(1:3),0.1,squeeze(framecol2(1,:,:))','-',0.005,0.009,[1:3]);hold on;
end
%         view(-120,40);
%         view(-108,46);
%         view(-71,18);
view(-112,46)
box on;
set(gca,'xtick',[],'ytick',[],'ztick',[]);
text(s(1).p(1).b(1),s(1).p(1).b(2) - 0.02,s(1).p(1).b(3)+0.05, '$$ \bf{\textit{A}_{1}^{(1)}} $$','Interpreter','latex', 'FontSize', 16)
text(s(2).p(1).b(1)+0.025,s(2).p(1).b(2) - 0.03,s(2).p(1).b(3)-0.025, '$$ \bf{\textit{A}_{1}^{(2)}} $$','Interpreter','latex', 'FontSize', 16)
text(s(3).p(1).b(1),s(3).p(1).b(2) - 0.035,s(3).p(1).b(3)-0.05, '$$ \bf{\textit{A}_{1}^{(3)}} $$','Interpreter','latex', 'FontSize', 16)
text(s(4).p(1).b(1) - 0.09,s(4).p(1).b(2) - 0.015,s(4).p(1).b(3)+0.0, '$$ \bf{\textit{A}_{1}^{(4)}} $$','Interpreter','latex', 'FontSize', 16)
text(s(1).p(2).b(1) - 0.03,s(1).p(2).b(2) - 0.03,s(1).p(2).b(3)-0.05, '$$ \bf{\textit{A}_{2}^{(1 \ldots 4)}} $$','Interpreter','latex', 'FontSize', 16)

orgPt = mean([s(1).p(1).b(1:3), s(1).p(2).b(1:3)],2);
endPt = s(1).p(2).A(1:3,1:3)*0.1*[-1;0;0];
mArrow3(orgPt,orgPt+endPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.005, 'facealpha',0.2);hold on;
%         endPt = s(1).p(2).A(1:3,1:3)*0.2*[1;0;0];
endPt = s(1).p(2).b(1:3);
mArrow3(orgPt,endPt,'color', 'k','stemWidth',0.014, 'tipWidth',0.0005, 'facealpha',0.2);hold on;

axis equal;

%         a2v = axis;
%
%         axis(a2v + [-0.0 0.0 -0.08 0.0 0 0]);
view(-112,46)
view(-161,34);
