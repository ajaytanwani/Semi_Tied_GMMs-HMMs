classdef auxFuns
%%% This class contains auxiliary functions used in other parts of the
%%% code. Feel free to modify the functions as per your convenience.
%
% Copyright (c) 2016 Idiap Research Institute
% Written by Ajay Tanwani

    methods (Static)
	%%
	function X = solveAlgebraicRiccati_eig(A, G, Q)
		%Solves the algebraic Riccati equation of the form A'X+XA'-XGX+Q=0, where X is symmetric with eigendecomposition.
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Danilo Bruno
		
		n = size(A,1);
		Z = [A -G; -Q -A']; %See Eq. (5.2.3) in doc/TechnicalReport.pdf
		
		[V,D] = eig(Z); %See Eq. (5.2.4) in doc/TechnicalReport.pdf
		U = [];
		for j=1:2*n
			if real(D(j,j)) < 0
				U = [U V(:,j)];
			end
		end
		
		X = U(n+1:end,1:n) / U(1:n,1:n); %See Eq. (5.2.5) in doc/TechnicalReport.pdf
		X = real(X);
	end
	%%
	function [DataR s model] = RemoveTimeDim(Data, s, model)
		% Call this function after initialization with time_based kmeans in order to remove the time dimension from the model and data
		%
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
		
		model.nbVarPos = model.nbVarPos - 1;
		model.nbVar = model.nbVar - 1; %Dimension of state vector
		
		DataR = Data(2:end,:,:);
		
		for n=1:model.nbSamples
			s(n).Data = s(n).Data(2:end,:,:);
			s(n).Data0 = s(n).Data0(2:end,:);
			s(n).Data1 = s(n).Data1(2:end,:);
			s(n).p(1).b = s(n).p(1).b(2:end);
			s(n).p(2).b = s(n).p(2).b(2:end);
			
			s(n).p(1).A = s(n).p(1).A(2:end,2:end);
			s(n).p(2).A = s(n).p(2).A(2:end,2:end);
		end
	end
	%%
	function Sigma = GenerateRandomSymmetricMatrix(N,c,dimNull,pivot)
		% This function is used for generating a random symmetric matrix of dimension N, scaled by c, with null space dimension number specified by dimNull along the directions of pivot
		
		% Copyright (c) 2016 Idiap Research Institute
		% Written by Ajay Tanwani <ajay.tanwani@idiap.ch>
		
		if nargin < 2
			c = 1;
			dimNull = 0;
			pivot = N;
		else if nargin < 3
				dimNull = 0;
				pivot = N;
			else if nargin < 4
					pivot = N;
				end
			end
		end
		
		
		% if nargin < 2
		%     c = 1;
		% end
		% dimNull
		d = c*rand(N,1); % The diagonal values
		if dimNull > 0
			%     tol = 1e-12;
			tol = 0;
			d(pivot-dimNull + 1:pivot) = tol;
		end
		t = triu(bsxfun(@min,d,d.').*rand(N),1); % The upper trianglar random values
		Sigma = diag(d)+t+t.'; % Put them together in a symmetric matrix
	end
	%%
	function [rx ry rz]= GetEulerAngles(R)
		
		% This function return the rotation along x,y and z direction from a
		% Rotation Matrix
		
		%Inputs:
		% R= 3x3 Rotation Matrix
		%Outputs:
		% rx= Rotation along x direction in radians
		% ry= Rotation along y direction in radians
		% rz= Rotation along z direction in radians
		
		%     R =
		%
		% [                           cos(ry)*cos(rz),                          -cos(ry)*sin(rz),          sin(ry)]
		% [ cos(rx)*sin(rz) + cos(rz)*sin(rx)*sin(ry), cos(rx)*cos(rz) - sin(rx)*sin(ry)*sin(rz), -cos(ry)*sin(rx)]
		% [ sin(rx)*sin(rz) - cos(rx)*cos(rz)*sin(ry), cos(rz)*sin(rx) + cos(rx)*sin(ry)*sin(rz),  cos(rx)*cos(ry)]
		
		% ry=asin(R(1,3));
		% rz=acos(R(1,1)/cos(ry));
		% rx=acos(R(3,3)/cos(ry));
		
		rx = atan2(R(3,2), R(3,3));
		ry = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
		rz = atan2(R(2,1), R(1,1));
	end
	end
end