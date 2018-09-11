function c = mmf_assign (S, m, v)

n = size (v, 2);
d = size (v, 1);
c = zeros (m, n, 'uint8');

% process separately each subquantizer
for i = 1:m
  
  % find the nearest centroid for each subvector
  
 vsub = (S{i}.subspace'*v);
%   S{i}.C
%   pause
%   
  [idx, dis] = yael_nn(S{i}.C, vsub, 1, 2);
 
  c(i, :) = idx - 1;
end