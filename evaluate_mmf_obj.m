function o=evaluate_mmf_obj(X,S)

[n p]=size(X);
m=size(S,2);
o=0;
for i=1:m
  u=S{i}.subspace;
  h=S{i}.H;
   
  o=o+trace(u'*X'*full(h)*full(h)'*X*u);

end
