function [outmatrix,weight] = decorrelation_transform_ROWOP(inmatrix,killrow,killelement,userow)
%
% AUTHOR
% K. A. Innanen 2022 (contact k.innanen@ucalgary.ca)
%
% Useage: [outmatrix,weight] = rowop(inmatrix,killrow,killelement,userow);
% Input:  inmatrix = input matrix
%         killrow = row # on which to zero a column element
%         killelement = column element to zero
%         userow = row to use in linear combination with killrow to kill
%         element
% Output: outmatrix = matrix after row operation
%         weight = weight used to kill element

% Initialize
outmatrix = inmatrix;
% Determine weight
weight = -inmatrix(killrow,killelement)/inmatrix(userow,killelement); 
% Replace row
outmatrix(killrow,:) = inmatrix(killrow,:) + weight*inmatrix(userow,:); 