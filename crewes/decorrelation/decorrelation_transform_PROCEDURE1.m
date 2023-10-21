function TT = decorrelation_transform_PROCEDURE1(PHIS,TT)
%
% AUTHOR
% K. A. Innanen 2022 (contact k.innanen@ucalgary.ca)
%
% Useage: TT = decorrelation_transform_PROCEDURE2(PHIS,TT);
% Input:  N x N covariance or Hessian matrix PHIS
%         N x N transformation matrix TT with lower triangular entries
%         included.
% Output: N x N transformation matrix TT with diagonal and upper triangular
%         entries filled in.

% Get sizes
[NN,jnk] = size(PHIS); 

% Column 1 
% Get vector v components
vf0 = zeros(NN,1); vf0(1) = 1; vg0 = TT(:,1); vg0(1) = 0;
% Determine quadratic equation coefficients
G2 = vf0'*PHIS*vf0; G1 = 2*PHIS(1,:)*vg0; G0 = vg0'*PHIS*vg0 - 1;
% Select positive root
t00 = (- G1 + sqrt( G1*G1-4*G2*G0 ))/( 2*G2 ); TT(1,1) = t00;

% Columns 2 - (N-1).
cc = [];
for ii=2:NN-1,
     %
     % Set up auxiliary matrices
     cv=PHIS*TT(:,ii-1); cc=[cc;cv'];                                    
     bb=cc(:,ii+1:NN)*TT(ii+1:NN,ii); 
     aa = [cc(:,1:ii),bb];  
     % Carry out row operations
     [jnk,coeffs]=decorrelation_transform_ELIMS(aa); 
     % Get vector v components
     vf=coeffs(:,1); vf=[vf;1;zeros(NN-ii,1)]; vg=coeffs(:,2); vg=[vg;0;TT(ii+1:NN,ii)];        
     % Determine quadratic equation coefficients
     G2=vf'*PHIS*vf; G1=2*vf'*PHIS*vg; G0=vg'*PHIS*vg-1; 
     % Choose positive root
     temp=(-G1+sqrt(G1*G1-4*G2*G0))/(2*G2);    
     TT(ii,ii)=temp; 
     % Fill in non-diagonal elements of column
     TT(1:ii-1,ii)=coeffs*[temp;1];  
     %
end

% Column N
% Set up auxiliary matrices
cv=PHIS*TT(:,NN-1); cc=[cc;cv'];
aa = [cc,zeros(NN-1,1)];
% Carry out row operations
[jnk,coeffs]=decorrelation_transform_ELIMS(aa);
% Gev vector v components
vf=coeffs(:,1); vf=[vf;1]; vg=coeffs(:,2); vg=[vg;0];
% Determine quadratic equation coefficients
G2=vf'*PHIS*vf; G1=2*vf'*PHIS*vg; G0=vg'*PHIS*vg-1; 
% Choose positive root
temp=(-G1+sqrt(G1*G1-4*G2*G0))/(2*G2);
% Fill in non-diagonal elements of column
TT(NN,NN)=temp; TT(1:NN-1,NN)=coeffs*[temp;1];
