function [ids, dis] = mmf_search (S, cbase, vquery, k,h,m)
    
n = size (cbase, 2);
nq = size (vquery, 2);
d = size (vquery, 1);

distab  = zeros (h, m, 'single');
dis = zeros (nq, k, 'single');
ids = zeros (nq, k, 'single');

for query = 1:nq

  % pre-compute the table of squared distance to centroids
  for q = 1:m
   
      vsub = (S{q}.subspace'*vquery(:,query));
   
      distab(:,q) = yael_L2sqr (vsub, S{q}.C)';
   
  end 

  % add the tabulated distances to construct the distance estimators
  disquerybase = sumidxtab (distab, cbase, 0);
 
  [dis1, ids1] = yael_kmin ((disquerybase), k);
  
  dis(query, :) = dis1;
  ids(query, :) = ids1;
end