function P = Func_vecperm(p, s)
%VECPERM    Vec-permutation matrix.
%           VECPERM(M, N) is the vec-permutation matrix, an MN-by-MN
%           permutation matrix P with the property that if A is M-by-N then
%           vec(A) = P*vec(A').
%           If N is omitted, it defaults to M.

%   P is formed by taking every n'th row from EYE(M*N), starting with
%   the first and working down - see p. 277 of the reference.

%   Reference:
%   H. V. Henderson and S. R. Searle The vec-permutation matrix,
%   the vec operator and Kronecker products: A review Linear and
%   Multilinear Algebra, 9 (1981), pp. 271-288.

if nargin == 1, s = p; end

P = zeros(p*s);
I = eye(p*s);

k = 1;
for i=1:s
    for j=i:s:p*s
        P(k,:) = I(j,:);
        k = k+1;
    end
end


% Example
% clear all
% 
% % number of stages
% s=3;
% 
% % number of patches
% p=2;
% 
% P = Func_vecperm(p,s);
% 
% % Demographic matrices for patches 
% z2 = zeros(2,2)
% z3 = zeros(3,3)
% 
% syms f11 f12 f13 s11 s12 f21 f22 f23 s21 s22
% B1 = [f11 f12 f13;
%       s11 0   0  ;
%       0   s12 0  ];
%   
% B2 = [f21 f22 f23;
%       s21 0   0  ;
%       0   s22 0  ];
%   
% B = [B1 z3;
%      z3 B2];
%     
% % Connectivity matrices for stages
% syms d11 d12 d22 d21
% M1 = [d11 d21 ;
%       d12 d22];
% M2 = eye(2);
% M3 = eye(2);
% 
% M = [M1 z2 z2 ;
%      z2 M2 z2 ;
%      z2 z2 M3];
% 
% A = P'*M*P*B