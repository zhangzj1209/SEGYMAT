% Calculate transform matrix for FWI Hessian using Procedure 2.
clear; close all;

% Load Hessian matrix determined from wave simulation
load input_EXAMPLE3.txt; hh = input_EXAMPLE3; [NN,jnk] = size(hh); 

% Preselect zeros for matrix TT1
TT1 = zeros(NN,NN);

% Determine transformation matrix using Procedure 1 
TT1 = decorrelation_transform_PROCEDURE1(hh,TT1);

% Check by transforming cv
hhT1 = TT1'*hh*TT1;

% Preselect zeros for matrix TT2
TT2 = zeros(NN,NN);

% Determine transformation matrix using Procedure 2 
TT2 = decorrelation_transform_PROCEDURE2(hh,TT2);

% Check by transforming cv
hhT2 = TT2'*hh*TT2;

%% Plot results
figure,

subplot(3,3,1),
imagesc(hh); colormap('gray');
title('FWI Hessian matrix');

subplot(3,3,2),
imagesc(TT1); colormap('gray');
title('T matrix Proc 12');

subplot(3,3,3),
imagesc(hhT1); colormap('gray');
title('Xformed hh Proc 1');

subplot(3,3,4),
imagesc(hh); colormap('gray');
title('FWI Hessian matrix');

subplot(3,3,5),
imagesc(TT2); colormap('gray');
title('T matrix Proc 2');

subplot(3,3,6),
imagesc(hhT2); colormap('gray');
title('Xformed hh Proc 2');

