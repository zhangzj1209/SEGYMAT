function [out,pivoterror] = decorrelation_transform_A2U(in);
%
% AUTHOR
% K. A. Innanen 2022 (contact k.innanen@ucalgary.ca)
%
% Useage: [out,pivoterror] = decorrelation_transform_A2U(in);
% in:  [N,N+2] fully populated matrix with no 0 pivots (tested)
% out: [N,N+2] with left [N,N] elements in upper triangular form
% error: if a zero pivot is found program dumps and writes offending row #
%

[N,jnk] = size(in);
out = in;
for ii=2:N,
    %
    index = ii-1;
    pivot = out(index,index);
    if ( abs(pivot)<10e-9 ) 
        pivoterror = ii;
        return; 
    end;
    for jj=ii:N,
        %
        factor = out(jj,index)/pivot;
        out(jj,:) = out(jj,:) - factor*out(index,:);
        %
    end
    %
end
pivoterror = 0;