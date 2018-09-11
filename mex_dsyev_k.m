function [eigsval eigsvec] = mex_dsyev_k(A,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if(nargin<0)
    A=rand(10,10);
end
n=size(A,1);
[eigsvec eigsval]=mex_dsyev_ex(A);
[maxv ind]=sort(eigsval,'descend');

eigsval=eigsval(ind(1:k));
eigsvec=eigsvec(:,ind(1:k));

end

