function [out1,out2] = decorrelation_transform_ELIMS(in);
%
% AUTHOR
% K. A. Innanen 2022 (contact k.innanen@ucalgary.ca)
%
% Useage: [out1,out2] = decorrelation_transform_ELIMS(in);
% in:  [N,N+2] fully populated matrix 
% out1: [N,3] matrix, out(ii,:) = 3 nonzero bottom elements derived by
% elimination of matrix with ii'th and N'th columns exchanged, i.e.,
% coefficients for determining x_N+1 in terms of x_ii.
% out2: ratios of coefficients f and g.
%
[N,jnk] = size(in); out1 = zeros(N,3); out2 = zeros(N,2);
for ii=1:N,
    %
    % Exchange ii'th and N'th columns
    temp = in; temp(:,ii) = in(:,N); temp(:,N) = in(:,ii);
    % Decompose leftmost NxN elements of temp into upper triangular form
    [temp,pivoterror] = decorrelation_transform_A2U(temp);
    if pivoterror>0 
        %
        fprintf('Pivot error on exchange %d in row %d',ii,pivoterror);
        fprintf('\n');
    end
    % Extract and save bottom nonzero 3 components A t_i^(j<i) + B t_i^i + C = 0.
    out1(ii,:) = [temp(N,N),temp(N,N+1),temp(N,N+2)];
    % Alternatively export the 2 coefficients t_i^(j<i) = (-B/A) t_i^i + (-C/A).
    out2(ii,:) = [-temp(N,N+1)/temp(N,N),-temp(N,N+2)/temp(N,N)];
    %
end