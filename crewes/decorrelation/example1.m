% Calculate transform matrix for MNIST covariance example using both 
% Procedure 1 and 2.
clear; close all;

% Load covariance matrix determined from MNIST data
load input_EXAMPLE1.txt; cv = input_EXAMPLE1; [NN,jnk] = size(cv);

% Preselect zeros for matrix TT1
TT1 = zeros(NN,NN);

% Determine transformation matrix using Procedure 1
TT1 = decorrelation_transform_PROCEDURE1(cv,TT1);

% Check by transforming cv
cvT1 = TT1'*cv*TT1;

% Preselect zeros for matrix TT2
TT2 = zeros(NN,NN);

% Determine transformation matrix using Procedure 2 
TT2 = decorrelation_transform_PROCEDURE2(cv,TT2);

% Check by transforming cv
cvT2 = TT2'*cv*TT2;

% Plot results
figure,

subplot(3,3,1),
imagesc(cv,[-0.1,0.1]); colormap('gray');
title('Lithofacies covariance matrix');

subplot(3,3,2),
imagesc(TT1); colormap('gray');
title('T matrix Proc 1');

subplot(3,3,3),
imagesc(cvT1,[-0.1,0.1]); colormap('gray');
title('Xformed cv Proc 1');

subplot(3,3,4),
imagesc(cv,[-0.1,0.1]); colormap('gray');
title('Lithofacies covariance matrix');

subplot(3,3,5),
imagesc(TT2); colormap('gray');
title('T matrix Proc 2');

subplot(3,3,6),
imagesc(cvT2,[-0.5,0.5]); colormap('gray');
title('Xformed cv Proc 2');
