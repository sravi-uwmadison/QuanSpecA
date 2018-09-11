function R=mmf_test_stats(nquery, k, ids_mmf,ids_gnd)
nn_ranks_pqc = zeros (nquery, 1);
hist_pqc = zeros (k+1, 1);
for i = 1:nquery
  gnd_ids = ids_gnd(:,i);
  
    %nn_pos = find (ids_mmf(i, :) == gnd_ids);
    nn_pos=intersect(int32(ids_mmf(:,i)), int32(gnd_ids),'stable');
    
    if length (nn_pos) >= 1
      nn_ranks_pqc (i) = find(ids_mmf(:, i)==nn_pos(1));
    else
      nn_ranks_pqc (i) = k + 1; 
    end
end
nn_ranks_pqc = sort(nn_ranks_pqc);
R=[];
for i = [1 2 5 10 20 50 100 200 500 1000 2000 5000 10000]
  if i <= k
    r_at_i = length (find (nn_ranks_pqc <= i & nn_ranks_pqc <= k))/nquery * 100;
    R=[R; r_at_i];
    fprintf ('r@%3d = %.3f\n', i, r_at_i); 
  end
end