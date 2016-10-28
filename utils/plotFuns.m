classdef plotFuns
%%% This class contains plotting functions used in other parts of the
%%% code. Feel free to modify the functions as per your convenience.
%
% Copyright (c) 2016 Idiap Research Institute
% Written by Sylvain Calinon and Ajay Tanwani

    methods (Static)
	%%
	function plotHSMM(Trans, StatesPriors, Pd)
		% Copyright (c) 2015 Idiap Research Institute
		% Sylvain Calinon, 2015
		
		nbStates = size(Trans,1);
		nbD = round(size(Pd,2)/1.15);
		valAlpha = 1;
		szPd = .5;
		graphRadius = 1;
		nodeRadius = .2;
		nodePts = 40;
		nodeAngle = linspace(pi/2,2*pi+pi/2,nbStates+1);
		for i=1:nbStates
			nodePos(:,i) = [cos(nodeAngle(i)); sin(nodeAngle(i))] * graphRadius;
		end
		clrmap = lines(nbStates);
		
		for i=1:nbStates
			%Plot StatesPriors
			posTmp = [cos(nodeAngle(i)); sin(nodeAngle(i))] * graphRadius + ...
				[cos(nodeAngle(i)+pi/3); sin(nodeAngle(i)+pi/3)] * nodeRadius*2;
			dirTmp = nodePos(:,i) - posTmp;
			dirTmp = (norm(dirTmp)-nodeRadius) * dirTmp/norm(dirTmp);
			plotFuns.plot2DArrow(posTmp, dirTmp, [.8 .8 .8]-StatesPriors(i)*.8);
			
			for j=i+1:nbStates
				%Plot Trans
				dirTmp = nodePos(:,j)-nodePos(:,i);
				dirTmp = (norm(dirTmp)-2*nodeRadius) * dirTmp/norm(dirTmp);
				offTmp = [dirTmp(2); -dirTmp(1)] / norm(dirTmp);
				posTmp = nodePos(:,i) + nodeRadius * dirTmp/norm(dirTmp) + offTmp*0.05;
				plotFuns.plot2DArrow(posTmp, dirTmp, [.8 .8 .8]-Trans(i,j)*.8);
			end
			for j=1:i
				%Plot Trans
				dirTmp = nodePos(:,j)-nodePos(:,i);
				dirTmp = (norm(dirTmp)-2*nodeRadius) * dirTmp/norm(dirTmp);
				offTmp = [dirTmp(2); -dirTmp(1)] / norm(dirTmp);
				posTmp = nodePos(:,i) + nodeRadius * dirTmp/norm(dirTmp) + offTmp*0.05;
				plotFuns.plot2DArrow(posTmp, dirTmp, [.8 .8 .8]-Trans(i,j)*.8);
			end
		end
		
		%Plot nodes
		for i=1:nbStates
			a = linspace(0,2*pi,nodePts);
			meshTmp = [cos(a); sin(a)] * nodeRadius + repmat(nodePos(:,i),1,nodePts);
			patch(meshTmp(1,:), meshTmp(2,:), clrmap(i,:),'edgecolor',clrmap(i,:)*0.5, 'facealpha', valAlpha,'edgealpha', valAlpha);
		end
		
		%Plot Pd
		for i=1:nbStates
			posTmp = [cos(nodeAngle(i)); sin(nodeAngle(i))] * graphRadius * 1.6;
			yTmp = Pd(i,1:nbD) / max(Pd(i,1:nbD));
			
			meshTmp = ([[0, linspace(0,1,nbD), 1]; [0, yTmp, 0]]-0.5) * szPd + repmat(posTmp,1,nbD+2);
			patch(meshTmp(1,:), meshTmp(2,:), clrmap(i,:), 'edgecolor',clrmap(i,:)*0.5, 'facealpha', valAlpha,'edgealpha', valAlpha);
			meshTmp = ([[0, 0, 1]; [1, 0, 0]]-0.5) * szPd + repmat(posTmp,1,3);
			plot(meshTmp(1,:), meshTmp(2,:), 'color',[0 0 0]);
		end
		
		axis equal;
	end
	%%
	function plotHMM(Trans, StatesPriors)
		% Copyright (c) 2015 Idiap Research Institute
		% Sylvain Calinon, 2015
		
		nbStates = size(Trans,1);
		valAlpha = 1;
		graphRadius = 1;
		nodeRadius = .2;
		nodePts = 40;
		nodeAngle = linspace(pi/2,2*pi+pi/2,nbStates+1);
		for i=1:nbStates
			nodePos(:,i) = [cos(nodeAngle(i)); sin(nodeAngle(i))] * graphRadius;
		end
		clrmap = lines(nbStates);
		
		for i=1:nbStates
			%Plot StatesPriors
			posTmp = [cos(nodeAngle(i)); sin(nodeAngle(i))] * graphRadius + ...
				[cos(nodeAngle(i)+pi/3); sin(nodeAngle(i)+pi/3)] * nodeRadius*2;
			dirTmp = nodePos(:,i) - posTmp;
			dirTmp = (norm(dirTmp)-nodeRadius) * dirTmp/norm(dirTmp);
			plot2DArrow(posTmp, dirTmp, [.8 .8 .8]-StatesPriors(i)*.8);
			
			%Plot self-transitions
			d = nodeRadius*1.2;
			posTmp = [cos(nodeAngle(i)); sin(nodeAngle(i))] * (graphRadius+d);
			R = nodeRadius;
			r = nodeRadius*0.5;
			aTmp = asin((4*d^2*R^2-(d^2-r^2+R^2)^2)^.5/(2*d*r));
			a = linspace(nodeAngle(i)+pi-aTmp, nodeAngle(i)-pi+aTmp, nodePts);
			meshTmp = [cos(a); sin(a)] * r + repmat(posTmp,1,nodePts);
			plot(meshTmp(1,:), meshTmp(2,:), 'color',[.8 .8 .8]-Trans(i,i)*.8,'linewidth',2);
			plot2DArrow(meshTmp(:,end-10), meshTmp(:,end)-meshTmp(:,end-10), [.8 .8 .8]-Trans(i,i)*.8);
			
			for j=i+1:nbStates
				%Plot Trans
				dirTmp = nodePos(:,j)-nodePos(:,i);
				dirTmp = (norm(dirTmp)-2*nodeRadius) * dirTmp/norm(dirTmp);
				offTmp = [dirTmp(2); -dirTmp(1)] / norm(dirTmp);
				posTmp = nodePos(:,i) + nodeRadius * dirTmp/norm(dirTmp) + offTmp*0.05;
				plot2DArrow(posTmp, dirTmp, [.8 .8 .8]-Trans(i,j)*.8);
			end
			for j=1:i
				%Plot Trans
				dirTmp = nodePos(:,j)-nodePos(:,i);
				dirTmp = (norm(dirTmp)-2*nodeRadius) * dirTmp/norm(dirTmp);
				offTmp = [dirTmp(2); -dirTmp(1)] / norm(dirTmp);
				posTmp = nodePos(:,i) + nodeRadius * dirTmp/norm(dirTmp) + offTmp*0.05;
				plot2DArrow(posTmp, dirTmp, [.8 .8 .8]-Trans(i,j)*.8);
			end
		end
		
		%Plot nodes
		for i=1:nbStates
			a = linspace(0,2*pi,nodePts);
			meshTmp = [cos(a); sin(a)] * nodeRadius + repmat(nodePos(:,i),1,nodePts);
			patch(meshTmp(1,:), meshTmp(2,:), clrmap(i,:),'edgecolor',clrmap(i,:)*0.5, 'facealpha', valAlpha,'edgealpha', valAlpha);
		end
		
		axis equal;
	end
	%%
	function h=plotGMM3D(Mu, Sigma, col1,alpha, dispOpt)
		
		if nargin<4
			alpha=1;
		end
		if nargin<5
			dispOpt=1;
		end
		
		nbData = size(Mu,2);
		nbPoints = 20; %nb of points to form a circular path
		nbRings = 10; %Number of circular paths following the principal direction
		
		pts0 = [cos(linspace(0,2*pi,nbPoints)); sin(linspace(0,2*pi,nbPoints))];
		
		h=[];
		for n=1:nbData
			[V0,D0] = eigs(Sigma(:,:,n));
			U0 = V0*D0^.5;
			
			ringpts0 = [cos(linspace(0,pi,nbRings+1)); sin(linspace(0,pi,nbRings+1))];
			ringpts = zeros(3,nbRings);
			ringpts([2,3],:) = ringpts0(:,1:nbRings);
			U = zeros(3);
			U(:,[2,3]) = U0(:,[2,3]);
			ringTmp = U*ringpts;
			
			%Compute touching circular paths
			for j=1:nbRings
				U = zeros(3);
				U(:,1) = U0(:,1);
				U(:,2) = ringTmp(:,j);
				pts = zeros(3,nbPoints);
				pts([1,2],:) = pts0;
				xring(:,:,j) = U*pts + repmat(Mu(:,n),1,nbPoints);
			end
			
			%Plot filled ellispoid
			xringfull = xring;
			xringfull(:,:,end+1) = xringfull(:,end:-1:1,1); %Close the ellipsoid
			for j=1:size(xringfull,3)-1
				for i=1:size(xringfull,2)-1
					xTmp = [xringfull(:,i,j) xringfull(:,i+1,j) xringfull(:,i+1,j+1) xringfull(:,i,j+1) xringfull(:,i,j)];
					if dispOpt==1
						%Version 1 (for GMM plots)
						h = [h patch(xTmp(1,:),xTmp(2,:),xTmp(3,:), min(col1+0.3,1),'edgecolor',col1,'linewidth',1,'facealpha',alpha,'edgealpha',alpha)]; %,'facealpha',0.5
					else
						%Version 2 (for GMR plots)
						patch(xTmp(1,:),xTmp(2,:),xTmp(3,:), min(col1+0.35,1),'linestyle','none','facealpha',alpha);
					end
				end
			end
		end
	end		
	%%
	function h = plotGMM(Mu, Sigma, color, valAlpha)
		% This function displays the parameters of a Gaussian Mixture Model (GMM).
		% Inputs -----------------------------------------------------------------
		%   o Mu:           D x K array representing the centers of K Gaussians.
		%   o Sigma:        D x D x K array representing the covariance matrices of K Gaussians.
		%   o color:        3 x 1 array representing the RGB color to use for the display.
		%   o valAlpha:     transparency factor (optional).
		%
		% Copyright (c) 2015 Idiap Research Institute
		% Sylvain Calinon, 2014

		nbStates = size(Mu,2);
		nbDrawingSeg = 35;
		lightcolor = min(color+0.5,1);
		t = linspace(-pi, pi, nbDrawingSeg);

		h=[];
		for i=1:nbStates
			R = real(sqrtm(1.0.*Sigma(:,:,i)));
			X = R * [cos(t); sin(t)] + repmat(Mu(:,i), 1, nbDrawingSeg);
			if nargin>3 %Plot with alpha transparency
				h = [h patch(X(1,:), X(2,:), lightcolor, 'lineWidth', 1, 'EdgeColor', color, 'facealpha', valAlpha,'edgealpha', valAlpha)];
				MuTmp = [cos(t); sin(t)] * 0.006 + repmat(Mu(:,i),1,nbDrawingSeg);
				h = [h patch(MuTmp(1,:), MuTmp(2,:), color, 'LineStyle', 'none', 'facealpha', valAlpha)];
			else %Plot without transparency
				h = [h patch(X(1,:), X(2,:), lightcolor, 'lineWidth', 1, 'EdgeColor', color)];
				h = [h plot(Mu(1,:), Mu(2,:), '.', 'lineWidth', 2, 'markersize', 6, 'color', color)];
			end
		end
	end
		
	%%
	function h = plot3Dframe2( rotAxis, orgPt, sf, colMat, linestyle, stem, tip, components)
		
		% This function allows to graph a nice 3D frame of reference using
		% different colors for each axis.
		%
		%   rotAxis:    The rotation matrix representing the 3D frame
		%   orgPt:      The origin of the frame
		%
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
		
		if nargin < 3
			sf = 1;
		end
		
		if nargin<4
			colMat = eye(3);
		end
		
		if nargin < 5
			linestyle = '-';
		end
		
		if nargin < 6
			stem = 0.01;
		end
		
		if nargin < 7
			tip = 0.015;
		end
		
		if nargin < 8
			components = 1:3;
		end
		
		endPt = rotAxis*(sf.*eye(3));
		
		for i = 1:length(components)
			h(i) = plotFuns.mArrow3(orgPt,orgPt+endPt(:,components(i)),'color', colMat(:,components(i)),'stemWidth',stem, 'tipWidth',tip, 'facealpha',0.75);hold on;
		end
	end
	%%
	function h = plot3Dframe( rotAxis, orgPt, colMat, linestyle )
		% This function allows to graph a nice 3D frame of reference using
		% different colors for each axis.
		%
		%   rotAxis:    The rotation matrix representing the 3D frame
		%   orgPt:      The origin of the frame
		%
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>

		if nargin<3
			colMat = eye(3);
		end

		if nargin < 4
			linestyle = '-';
		end
		h(1) = quiver3( orgPt(1), orgPt(2), orgPt(3), rotAxis(1,1), rotAxis(2,1), rotAxis(3,1), 0.2, 'linewidth', 2, 'LineStyle', linestyle,'color', colMat(:,1)); hold on;
		% h(2) = quiver3( orgPt(1), orgPt(2), orgPt(3), rotAxis(1,2), rotAxis(2,2), rotAxis(3,2), 0.2, 'linewidth', 2, 'LineStyle', linestyle, 'color', colMat(:,2)); hold on;
		h(3) = quiver3( orgPt(1), orgPt(2), orgPt(3), rotAxis(1,3), rotAxis(2,3), rotAxis(3,3), 0.2, 'linewidth', 2, 'LineStyle', linestyle, 'color', colMat(:,3));

		% h(1) = quiver3( orgPt(1), orgPt(2), orgPt(3), rotAxis(1,1), rotAxis(2,1), rotAxis(3,1), 0.2, 'linewidth', 2, 'color', colMat(1,:));
		% h(2) = quiver3( orgPt(1), orgPt(2), orgPt(3), rotAxis(1,2), rotAxis(2,2), rotAxis(3,2), 0.2, 'linewidth', 2, 'color', colMat(2,:));
		% h(3) = quiver3( orgPt(1), orgPt(2), orgPt(3), rotAxis(1,3), rotAxis(2,3), rotAxis(3,3), 0.2, 'linewidth', 2, 'color', colMat(3,:));		
	end
	%%
	function h = plot2DArrow(pos,dir,col)
		% Copyright (c) 2015 Idiap Research Institute
		% Sylvain Calinon, 2015
		h(1) = plot([pos(1) pos(1)+dir(1)], [pos(2) pos(2)+dir(2)], '-','linewidth',2,'color',col);
		sz = 8E-2;
		pos = pos+dir;
		if norm(dir)>sz
			dir = dir/norm(dir);
			prp = [dir(2); -dir(1)];
			dir = dir*sz;
			prp = prp*sz;
			msh = [pos pos-dir-prp/2 pos-dir+prp/2 pos];
			h(2) = patch(msh(1,:),msh(2,:),col,'edgecolor',col);
		end
	end
	%%
	function drawFrame(skel, mot, current_frame, color)
		% This function is used for plotting the skeleton trajectories
		
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
		
		npaths = size(skel.paths,1);
		skel_lines = cell(npaths,1);
		for k = 1:npaths
			path = skel.paths{k};
			nlines = length(path)-1;
			px = zeros(2,1); py = zeros(2,1); pz = zeros(2,1);
			px(1) = mot.jointTrajectories{path(1)}(1,current_frame);
			py(1) = mot.jointTrajectories{path(1)}(2,current_frame);
			pz(1) = mot.jointTrajectories{path(1)}(3,current_frame);
			for j = 2:nlines
				px(2) = mot.jointTrajectories{path(j)}(1,current_frame);
				py(2) = mot.jointTrajectories{path(j)}(2,current_frame);
				pz(2) = mot.jointTrajectories{path(j)}(3,current_frame);
				skel_lines{k}(j-1) = line(px,py,pz,'color',color,'linewidth',2);
				px(1) = px(2);
				py(1) = py(2);
				pz(1) = pz(2);
			end
			px(2) = mot.jointTrajectories{path(nlines+1)}(1,current_frame);
			py(2) = mot.jointTrajectories{path(nlines+1)}(2,current_frame);
			pz(2) = mot.jointTrajectories{path(nlines+1)}(3,current_frame);
			skel_lines{k}(nlines) = line(px,py,pz,'color',color,'linewidth',2);
			hold on
		end
	end
	
	%% 
	function DrawModel(skel, mot, Data1, Data2, DataOut)
		% This function is used for plotting the skeleton trajectories
		
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
		
		for m=1:mot(1).njoints
			mot(1).jointTrajectories{m} = Data1((m-1)*3+1:m*3,:);
			mot(2).jointTrajectories{m} = Data2((m-1)*3+1:m*3,:);
			mot(3).jointTrajectories{m} = DataOut((m-1)*3+1:m*3,:);
		end
		
		TimeDemos = mot(2).nframes/mot(2).samplingRate;
		
		figure('PaperPosition',[0 0 5 3],'position',[0,0,1000,600],'color',[1 1 1]);
		set(0,'DefaultAxesLooseInset',[0,0,0,0]);
		xx = round(linspace(1,64,3));
		clrmap = colormap('jet');
		clrmap = min(clrmap(xx,:),.95);
		
		nbData = size(Data1,2);
		pi=1;
		for t=round(linspace(1,nbData,15))
			plotFuns.subaxis(3,5,pi,'spacing',0); hold on; axis off;
			pi = pi+1;
			plotFuns.drawFrame(skel(1), mot(1), t, clrmap(1,:)); % blue
			
			%                 hold on
			plotFuns.drawFrame(skel(2), mot(2), t, clrmap(2,:)); % green
			%                 hold on
			plotFuns.drawFrame(skel(3), mot(3), t, clrmap(3,:)); % red
			%                 hold on
			axis equal;
			% pause(0.1);
			axis([-20 60 0 60 -40 40]);
			
			h = title(['$$t= $$' num2str((TimeDemos/(nbData-1))*(t-1),2)],'Interpreter','latex','fontsize',12);
			pos = get(h,'position');
			
			set(h,'position',[pos(1) pos(2)-100 pos(3)],  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
			
			
			view(gca,-128,43);
			camup(gca,[0 1 0]);
			
		end
		
	end
	%%	
	function h = mArrow3(p1,p2,varargin)
		%mArrow3 - plot a 3D arrow as patch object (cylinder+cone)
		%
		% syntax:   h = mArrow3(p1,p2)
		%           h = mArrow3(p1,p2,'propertyName',propertyValue,...)
		%
		% with:     p1:         starting point
		%           p2:         end point
		%           properties: 'color':      color according to MATLAB specification
		%                                     (see MATLAB help item 'ColorSpec')
		%                       'stemWidth':  width of the line
		%                       'tipWidth':   width of the cone
		%
		%           Additionally, you can specify any patch object properties. (For
		%           example, you can make the arrow semitransparent by using
		%           'facealpha'.)
		%
		% example1: h = mArrow3([0 0 0],[1 1 1])
		%           (Draws an arrow from [0 0 0] to [1 1 1] with default properties.)
		%
		% example2: h = mArrow3([0 0 0],[1 1 1],'color','red','stemWidth',0.02,'facealpha',0.5)
		%           (Draws a red semitransparent arrow with a stem width of 0.02 units.)
		%
		% hint:     use light to achieve 3D impression
		%
		% Source: http://ch.mathworks.com/matlabcentral/fileexchange/25372-marrow3-m-easy-to-use-3d-arrow
		
		propertyNames = {'edgeColor'};
		propertyValues = {'none'};
		
		% evaluate property specifications
		for argno = 1:2:nargin-2
			switch varargin{argno}
				case 'color'
					propertyNames = {propertyNames{:},'facecolor'};
					propertyValues = {propertyValues{:},varargin{argno+1}};
				case 'stemWidth'
					if isreal(varargin{argno+1})
						stemWidth = varargin{argno+1};
					else
						warning('mArrow3:stemWidth','stemWidth must be a real number');
					end
				case 'tipWidth'
					if isreal(varargin{argno+1})
						tipWidth = varargin{argno+1};
					else
						warning('mArrow3:tipWidth','tipWidth must be a real number');
					end
				otherwise
					propertyNames = {propertyNames{:},varargin{argno}};
					propertyValues = {propertyValues{:},varargin{argno+1}};
			end
		end
		
		% default parameters
		if ~exist('stemWidth','var')
			ax = axis;
			if numel(ax)==4
				stemWidth = norm(ax([2 4])-ax([1 3]))/300;
			elseif numel(ax)==6
				stemWidth = norm(ax([2 4 6])-ax([1 3 5]))/300;
			end
		end
		if ~exist('tipWidth','var')
			tipWidth = 3*stemWidth;
		end
		tipAngle = 22.5/180*pi;
		tipLength = tipWidth/tan(tipAngle/2);
		ppsc = 50;  % (points per small circle)
		ppbc = 250; % (points per big circle)
		
		% ensure column vectors
		p1 = p1(:);
		p2 = p2(:);
		
		% basic lengths and vectors
		x = (p2-p1)/norm(p2-p1); % (unit vector in arrow direction)
		y = cross(x,[0;0;1]);    % (y and z are unit vectors orthogonal to arrow)
		if norm(y)<0.1
			y = cross(x,[0;1;0]);
		end
		y = y/norm(y);
		z = cross(x,y);
		z = z/norm(z);
		
		% basic angles
		theta = 0:2*pi/ppsc:2*pi; % (list of angles from 0 to 2*pi for small circle)
		sintheta = sin(theta);
		costheta = cos(theta);
		upsilon = 0:2*pi/ppbc:2*pi; % (list of angles from 0 to 2*pi for big circle)
		sinupsilon = sin(upsilon);
		cosupsilon = cos(upsilon);
		
		% initialize face matrix
		f = NaN([ppsc+ppbc+2 ppbc+1]);
		
		% normal arrow
		if norm(p2-p1)>tipLength
			% vertices of the first stem circle
			for idx = 1:ppsc+1
				v(idx,:) = p1 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
			end
			% vertices of the second stem circle
			p3 = p2-tipLength*x;
			for idx = 1:ppsc+1
				v(ppsc+1+idx,:) = p3 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
			end
			% vertices of the tip circle
			for idx = 1:ppbc+1
				v(2*ppsc+2+idx,:) = p3 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
			end
			% vertex of the tiptip
			v(2*ppsc+ppbc+4,:) = p2;
			
			% face of the stem circle
			f(1,1:ppsc+1) = 1:ppsc+1;
			% faces of the stem cylinder
			for idx = 1:ppsc
				f(1+idx,1:4) = [idx idx+1 ppsc+1+idx+1 ppsc+1+idx];
			end
			% face of the tip circle
			f(ppsc+2,:) = 2*ppsc+3:(2*ppsc+3)+ppbc;
			% faces of the tip cone
			for idx = 1:ppbc
				f(ppsc+2+idx,1:3) = [2*ppsc+2+idx 2*ppsc+2+idx+1 2*ppsc+ppbc+4];
			end
			
			% only cone v
		else
			tipWidth = 2*sin(tipAngle/2)*norm(p2-p1);
			% vertices of the tip circle
			for idx = 1:ppbc+1
				v(idx,:) = p1 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
			end
			% vertex of the tiptip
			v(ppbc+2,:) = p2;
			% face of the tip circle
			f(1,:) = 1:ppbc+1;
			% faces of the tip cone
			for idx = 1:ppbc
				f(1+idx,1:3) = [idx idx+1 ppbc+2];
			end
		end
		
		% draw
		fv.faces = f;
		fv.vertices = v;
		h = patch(fv);
		for propno = 1:numel(propertyNames)
			try
				set(h,propertyNames{propno},propertyValues{propno});
			catch
				disp(lasterr)
			end
		end
	end
	
	%% 
	function ArgStruct=parseArgs(args,ArgStruct,varargin)
		% Helper function for parsing varargin.
		%
		%
		% ArgStruct=parseArgs(varargin,ArgStruct[,FlagtypeParams[,Aliases]])
		%
		% * ArgStruct is the structure full of named arguments with default values.
		% * Flagtype params is params that don't require a value. (the value will be set to 1 if it is present)
		% * Aliases can be used to map one argument-name to several argstruct fields
		%
		%
		% example usage:
		% --------------
		% function parseargtest(varargin)
		%
		% %define the acceptable named arguments and assign default values
		% Args=struct('Holdaxis',0, ...
		%        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
		%        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
		%        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
		%        'rows',[],'cols',[]);
		%
		% %The capital letters define abrreviations.
		% %  Eg. parseargtest('spacingvertical',0) is equivalent to  parseargtest('sv',0)
		%
		% Args=parseArgs(varargin,Args, ... % fill the arg-struct with values entered by the user
		%           {'Holdaxis'}, ... %this argument has no value (flag-type)
		%           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
		%
		% disp(Args)
		%
		% Copyright (c) 2014, Aslak Grinsted
		% All rights reserved (BSD license).
		
		Aliases={};
		FlagTypeParams='';
		
		if (length(varargin)>0)
			FlagTypeParams=strvcat(varargin{1});
			if length(varargin)>1
				Aliases=varargin{2};
			end
		end
		
		
		%---------------Get "numeric" arguments
		NumArgCount=1;
		while (NumArgCount<=size(args,2))&(~ischar(args{NumArgCount}))
			NumArgCount=NumArgCount+1;
		end
		NumArgCount=NumArgCount-1;
		if (NumArgCount>0)
			ArgStruct.NumericArguments={args{1:NumArgCount}};
		else
			ArgStruct.NumericArguments={};
		end
		
		
		%--------------Make an accepted fieldname matrix (case insensitive)
		Fnames=fieldnames(ArgStruct);
		for i=1:length(Fnames)
			name=lower(Fnames{i,1});
			Fnames{i,2}=name; %col2=lower
			AbbrevIdx=find(Fnames{i,1}~=name);
			Fnames{i,3}=[name(AbbrevIdx) ' ']; %col3=abreviation letters (those that are uppercase in the ArgStruct) e.g. SpacingHoriz->sh
			%the space prevents strvcat from removing empty lines
			Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams)); %Does this parameter have a value? (e.g. not flagtype)
		end
		FnamesFull=strvcat(Fnames{:,2});
		FnamesAbbr=strvcat(Fnames{:,3});
		
		if length(Aliases)>0
			for i=1:length(Aliases)
				name=lower(Aliases{i,1});
				FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
				if isempty(FieldIdx)
					FieldIdx=strmatch(name,FnamesFull); %&??????? exact or not?
				end
				Aliases{i,2}=FieldIdx;
				AbbrevIdx=find(Aliases{i,1}~=name);
				Aliases{i,3}=[name(AbbrevIdx) ' ']; %the space prevents strvcat from removing empty lines
				Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
			end
			%Append aliases to the end of FnamesFull and FnamesAbbr
			FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1}));
			FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3}));
		end
		
		%--------------get parameters--------------------
		l=NumArgCount+1;
		while (l<=length(args))
			a=args{l};
			if ischar(a)
				paramHasValue=1; % assume that the parameter has is of type 'param',value
				a=lower(a);
				FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
				if isempty(FieldIdx)
					FieldIdx=strmatch(a,FnamesFull);
				end
				if (length(FieldIdx)>1) %shortest fieldname should win
					[mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));
					FieldIdx=FieldIdx(mxi);
				end
				if FieldIdx>length(Fnames) %then it's an alias type.
					FieldIdx=Aliases{FieldIdx-length(Fnames),2};
				end
				
				if isempty(FieldIdx)
					error(['Unknown named parameter: ' a])
				end
				for curField=FieldIdx' %if it is an alias it could be more than one.
					if (Fnames{curField,4})
						val=args{l+1};
					else
						val=1; %parameter is of flag type and is set (1=true)....
					end
					ArgStruct.(Fnames{curField,1})=val;
				end
				l=l+1+Fnames{FieldIdx(1),4}; %if a wildcard matches more than one
			else
				error(['Expected a named parameter: ' num2str(a)])
			end
		end
	end
	
	%%
	function h=subaxis(varargin)
		%SUBAXIS Create axes in tiled positions. (just like subplot)
		%   Usage:
		%      h=subaxis(rows,cols,cellno[,settings])
		%      h=subaxis(rows,cols,cellx,celly[,settings])
		%      h=subaxis(rows,cols,cellx,celly,spanx,spany[,settings])
		%
		% SETTINGS: Spacing,SpacingHoriz,SpacingVert
		%           Padding,PaddingRight,PaddingLeft,PaddingTop,PaddingBottom
		%           Margin,MarginRight,MarginLeft,MarginTop,MarginBottom
		%           Holdaxis
		%
		%           all units are relative (e.g from 0 to 1)
		%
		%           Abbreviations of parameters can be used.. (Eg MR instead of MarginRight)
		%           (holdaxis means that it wont delete any axes below.)
		%
		%
		% Example:
		%
		%   >> subaxis(2,1,1,'SpacingVert',0,'MR',0);
		%   >> imagesc(magic(3))
		%   >> subaxis(2,'p',.02);
		%   >> imagesc(magic(4))
		%
		% Copyright (c) 2014, Aslak Grinsted
		% All rights reserved (BSD license).
		f=gcf;				
		Args=[];
		UserDataArgsOK=0;
		Args=get(f,'UserData');
		if isstruct(Args)
			UserDataArgsOK=isfield(Args,'SpacingHorizontal')&isfield(Args,'Holdaxis')&isfield(Args,'rows')&isfield(Args,'cols');
		end
		OKToStoreArgs=isempty(Args)|UserDataArgsOK;
		
		if isempty(Args) %&(~UserDataArgsOK)
			Args=struct('Holdaxis',0, ...
				'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
				'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
				'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
				'rows',[],'cols',[]);
		end
		Args=plotFuns.parseArgs(varargin,Args,{'Holdaxis'},{'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
		
		if (length(Args.NumericArguments)>1)
			Args.rows=Args.NumericArguments{1};
			Args.cols=Args.NumericArguments{2};
			%remove these 2 numerical arguments
			Args.NumericArguments={Args.NumericArguments{3:end}};
		end
		
		if OKToStoreArgs
			set(f,'UserData',Args);
		end
		
		switch length(Args.NumericArguments)
			case 0
				return % no arguments but rows/cols....
			case 1
				x1=mod((Args.NumericArguments{1}-1),Args.cols)+1; x2=x1;
				y1=floor((Args.NumericArguments{1}-1)/Args.cols)+1; y2=y1;
			case 2
				x1=Args.NumericArguments{1};x2=x1;
				y1=Args.NumericArguments{2};y2=y1;
			case 4
				x1=Args.NumericArguments{1};x2=x1+Args.NumericArguments{3}-1;
				y1=Args.NumericArguments{2};y2=y1+Args.NumericArguments{4}-1;
			otherwise
				error('subaxis argument error')
		end
		
		
		cellwidth=((1-Args.MarginLeft-Args.MarginRight)-(Args.cols-1)*Args.SpacingHorizontal)/Args.cols;
		cellheight=((1-Args.MarginTop-Args.MarginBottom)-(Args.rows-1)*Args.SpacingVertical)/Args.rows;
		xpos1=Args.MarginLeft+Args.PaddingLeft+cellwidth*(x1-1)+Args.SpacingHorizontal*(x1-1);
		xpos2=Args.MarginLeft-Args.PaddingRight+cellwidth*x2+Args.SpacingHorizontal*(x2-1);
		ypos1=Args.MarginTop+Args.PaddingTop+cellheight*(y1-1)+Args.SpacingVertical*(y1-1);
		ypos2=Args.MarginTop-Args.PaddingBottom+cellheight*y2+Args.SpacingVertical*(y2-1);
		
		if Args.Holdaxis
			h=axes('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
		else
			h=subplot('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
		end
		
		
		set(h,'box','on');
		%h=axes('position',[x1 1-y2 x2-x1 y2-y1]);
		set(h,'units',get(gcf,'defaultaxesunits'));
		set(h,'tag','subaxis');
		
		
		
		%if (nargout==0) clear h; end;
	end
	end
end