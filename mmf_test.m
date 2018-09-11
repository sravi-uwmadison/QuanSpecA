function R=mmf_test(m,h,k,dim)

if(nargin<4)
    [Xtrain Xbase Xquery ids_gnd]=generate_data('random',16);
else
    [Xtrain Xbase Xquery ids_gnd]=generate_data('random',dim);
end
if(nargin==0)
    m=8;
    h=256; %%should be power of 2
    k=10;
end


 mu = mean(Xtrain, 2);
% 


[S objkmeans]=call_mmf(single(Xtrain)',m,h,5);


cbase=mmf_assign(S,m,Xbase);


[ids_mmf, dis_mmf] = mmf_search(S, cbase, Xquery, k, h, m);

R_mmf=mmf_test_stats(size(Xquery,2), k, ids_mmf',ids_gnd(1,:));


 