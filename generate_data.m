function [vtrain vbase vquery ids_gnd]=generate_data(dataset,p)
if strcmp (dataset, 'random')
  % synthetic dataset
  if(nargin<2)
      d=16;
  else
    d =p;
  end
    % Generate a set of unit norm vectors
    ntrain = 1000;
    nbase = 1000; 
    nquery = 1000;
 
    vtrain = single (rand (d, ntrain));
    vbase = single (rand(d, nbase));
    vquery = single (rand(d, nquery)); 

    % Compute the ground-truth
    
   [ids_gnd, dis_gnd] = yael_nn(vbase, vquery, 10,2);

end


