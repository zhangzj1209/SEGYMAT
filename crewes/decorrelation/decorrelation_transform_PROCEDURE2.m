function TT = decorrelation_transform_PROCEDURE2(PHIS,TT);
%
% AUTHOR
% K. A. Innanen 2022 (contact k.innanen@ucalgary.ca)
%
% Useage: TT = decorrelation_transform_PROCEDURE2(PHIS,TT,root);
% Input:  N x N covariance or Hessian matrix PHIS
%         N x N transformation matrix TT with lower triangular entries
%         included.
% Output: N x N transformation matrix TT with diagonal and upper triangular
%         entries filled in.

% Get sizes 
[NN,jnk] = size(PHIS);

% Prep iterations
AA = []; ww = [];

% Do JJ=1 outside loop
ua = zeros(NN,1); ua(1) = 1;
ub = zeros(NN,1); ub(2:NN) = TT(2:NN,1);
alpha = ua'*PHIS*ua; beta = 2*ua'*PHIS*ub; gamma = ub'*PHIS*ub-1;
TT(1,1) = (-beta+sqrt(beta*beta-4*alpha*gamma))/(2*alpha);

% Do JJ=2 outside loop
JJ=2;
% Set up new equation info
AA = [AA;TT(:,JJ-1)'*PHIS];
BB = AA(1,JJ+1:NN)*TT(JJ+1:NN,JJ);
% Insert into core eqn system
CC = [AA(:,1:JJ),BB];
% Do row reductions 
% (none needed)
% Find new JJ-1 length vectors f, g
FF = -CC(:,2)./CC(:,1); GG = -CC(:,3)./CC(:,1);
% Assemble into new ua, ub
ua = [FF;1;zeros(NN-JJ,1)]; ub = [GG;0;TT(JJ+1:NN,JJ)];
% Find quadratic equation for T(JJ,JJ)
alpha = ua'*PHIS*ua; beta = 2*ua'*PHIS*ub; gamma = ub'*PHIS*ub-1;
TT(JJ,JJ) = (-beta+sqrt(beta*beta-4*alpha*gamma))/(2*alpha);
% Find transform elements above the diagonal
TT(1:JJ-1,JJ) = [FF,GG]*[TT(JJ,JJ);1];
%TT'*PHIS*TT
%return

% Do JJ=3 outside loop
JJ=3;
% 1. Add new bottom row to AA
AA = [AA;TT(:,JJ-1)'*PHIS];
% 2. Compute BB
BB = AA(:,JJ+1:NN)*TT(JJ+1:NN,JJ);
% 3. Adjust previous CC and augment
% a. First, cut off rightmost column of old CC
CC = CC(1:JJ-2,1:JJ-1);
% b. Then, grab two new columns, new A and new B with bottom elements
% missing, and put together in a JJ-2 x 2 matrix
newC1 = AA(1:JJ-2,JJ); newC2 = BB(1:JJ-2); newC = [newC1,newC2];
% c. Apply all previous row ops to them 
%    (none yet)
% d. Append to CC
CC = [CC,newC];
% e. Then grab new row from bottom of AA, with bottom element of BB
newR = [AA(JJ-1,1:JJ),BB(JJ-1)];
% f. Append to CC 
CC = [CC;newR];
% 4. New row ops.  Different for JJ=3!
% Eliminate first element of second row
[CC,weight] = decorrelation_transform_ROWOP(CC,2,1,1); ww = [ww;[weight,2,1]]; 
% Eliminate second element of 1st row
[CC,weight] = decorrelation_transform_ROWOP(CC,1,2,2); ww = [ww;[weight,1,2]];
% 5. Determine vectors f and g from final (JJ-1)x(JJ+1) CC system
CCd = diag(CC(:,1:JJ-1));
FF = -CC(:,JJ)./CCd; GG = -CC(:,JJ+1)./CCd;
% 6. Compute ua, ub, then solve for TT(JJ,JJ).
ua = [FF;1;zeros(NN-JJ,1)]; ub = [GG;0;TT(JJ+1:NN,JJ)];
alpha = ua'*PHIS*ua; beta = 2*ua'*PHIS*ub; gamma = ub'*PHIS*ub-1;
TT(JJ,JJ) = (-beta+sqrt(beta*beta-4*alpha*gamma))/(2*alpha);
% 7. Find transform elements above the diagonal
TT(1:JJ-1,JJ) = [FF,GG]*[TT(JJ,JJ);1];

% JJ = 4 to NN
if ( NN>3 )
    %    
    for JJ=4:NN,
        %
        % 1. Add new bottom row to AA
        AA = [AA;TT(:,JJ-1)'*PHIS];
        % 2. Compute BB
        BB = AA(:,JJ+1:NN)*TT(JJ+1:NN,JJ);
        % 3. Adjust previous CC and augment
        % a. First, cut off rightmost column of old CC
        CC = CC(1:JJ-2,1:JJ-1);
        % b. Grab two new columns, new A and new B with bottom elements missing
        newC1 = AA(1:JJ-2,JJ); newC2 = BB(1:JJ-2); newC = [newC1,newC2];
        % c. Apply all previous row ops to these
        % ww useage: [weight,killrow,userow]
        [wr,wc] = size(ww);
        for ii=1:wr,
            %
            newC(ww(ii,2),:) = newC(ww(ii,2),:) + ww(ii,1)*newC(ww(ii,3),:);
            %
        end
        % d. Append to CC
        CC = [CC,newC];
        % e. Grab new row from bottom of AA, with bottom element of BB
        newR = [AA(JJ-1,1:JJ),BB(JJ-1)];
        % f. Append to CC 
        CC = [CC;newR];
        % 4. a. Row reduce new bottom row until 3 rightmost elements remain.  Use upper rows one at a time.
        for ii=1:JJ-2,
            %
            % Useage: decorrelation_transform_ROWOP(inmatrix,killrow,killelement,userow)
            [CC,weight] = decorrelation_transform_ROWOP(CC,JJ-1,ii,ii);
            % Useage: ww = [weight,killrow,userow]
            ww = [ww;[weight,JJ-1,ii]];
            %
        end
        % 4. b. Row reduce third from right column until bottom element remains.
        for ii=1:JJ-2,
            %
            % Useage: decorrelation_transform_ROWOP(inmatrix,killrow,killelement,userow)
            [CC,weight] = decorrelation_transform_ROWOP(CC,ii,JJ-1,JJ-1);
            % Useage: ww = [weight,killrow,userow]
            ww = [ww;[weight,ii,JJ-1]];
            %
        end
        % 5. Determine vectors f and g from final (JJ-1)x(JJ+1) CC system
        CCd = diag(CC(:,1:JJ-1));
        FF = -CC(:,JJ)./CCd; GG = -CC(:,JJ+1)./CCd;
        % 6. Compute ua, ub, then solve for TT(JJ,JJ).
        ua = [FF;1;zeros(NN-JJ,1)]; ub = [GG;0;TT(JJ+1:NN,JJ)];
        alpha = ua'*PHIS*ua; beta = 2*ua'*PHIS*ub; gamma = ub'*PHIS*ub-1;
        TT(JJ,JJ) = (-beta+sqrt(beta*beta-4*alpha*gamma))/(2*alpha);
        % 7. Find transform elements above the diagonal
        TT(1:JJ-1,JJ) = [FF,GG]*[TT(JJ,JJ);1];
        %
    end
    %
end