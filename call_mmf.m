function [S prevX]=call_mmf(X,m,h,maxitr,Hint)
%%Dataset X
%%m : number of subspaces
%%h : number of cluster centers
%% maxiter: maximum iteration
%%Hint: initialization of H


%finding initilizatioprev of H
if(nargin<4)
    maxitr=1;
end
[n p]=size(X);

if(nargin<5)
   %[Hint Sini]=intializeH(X,p,h,m);%% if initailized by eigen_allocation method
  [Hint Sini len0 len1]=intializeH(X,p,h,m);%% if initailized by ckmeans  method
end

% % finding intilization to subspaces
 C=X'*Hint*Hint'*X;
 C=double(C+C')/2; 

% 
 iniC=C;
count=1;
if(mod(p,m)~=0)
    display('number of subspace should divide dimensions');
    return
end
prevX=[];


 for i=1:m 
  %S{count}.subspace=Sini(:,i:i+p/m-1);  %% if initailized by eigen_allocation method
  S{count}.subspace=Sini(:,len0(i):len1(i));%% if initailized by ckmeans method
  prevX(:,:,count)=S{count}.subspace*S{count}.subspace';
  count=count+1;
  
end

% display('Initial Subspaces done');
%
 
for i=1:m
    R_this = S{i}.subspace;
  
    [H I centroids_tmp]=find_best_kmeans(X,h,R_this);
    S{i}.idx=I;
    S{i}.C=centroids_tmp;
    
    S{i}.H=H;
    
end
%

obj_bef=evaluate_mmf_obj(X,S);
  
E=[];
oldS=S;
seq=0;
Quant_err=[];
epsilon=.00000001;
% starting subsequent iterations
for iter = 1 : maxitr %itr number
  
    %finding new Hs
    obj=0;
   
    
    for i=1:m
        
        R_this = S{i}.subspace;
      
        [H I centroids_tmp]=find_best_kmeans(X,h,R_this);
        S{i}.idx=I;
        S{i}.C=centroids_tmp;
        
        S{i}.H=H;
        
    end
    
    allC=[];
    for i=1:m
        H_this=S{i}.H;
        
        C=double(X'*H_this*H_this'*X);
        iniC=C;
        allC=cat(3,allC,C);
    end
 
    %finding new subspaces
        [S_this prevX obj_this]=optimize_subspaces_blockcord(allC,p/m,seq,iniC,prevX);
    
      
    for i=1:m

        S{i}.subspace=S_this(:,:,i);
    
    end
    

    obj_this=evaluate_mmf_obj(X,S);
  
    E=[E; obj_this];
    if(iter>1)
        if(E(iter)<E(iter-1))
            S=oldS;
            E(end)=[];
            break
        end
        
        if(abs(E(iter)-E(iter-1))<epsilon)
            S=oldS;
            E(end)=[];
            break
        end
    end
        if (iter==1)&&(obj_bef>E(1))
              S=oldS;
    
            break
        end
    
    
     oldS=S;
end
 

function H=converttoH(idxs,n,h)
H=zeros(n,h);
for i=1:length(idxs)
    H(i,idxs(i))=1;
end

for j=1:size(H,2)
    if(numel(find(H(:,j))>0))
     H(:,j)=H(:,j)./norm(H(:,j),2);
    end
end



function C=clean(C)

[U S]=svd(C);

ind=find(diag(S)<=0.01);

for i=1:length(ind)
    S(ind(i),ind(i))=0;
end

C=real(real(U)*real(S)*real(U)');
C=(C+C')/2;


function [H I C]=find_best_kmeans(X,h,R)
sumds_best=0;
Xproj=X*R;
[n p]=size(X);
opts = statset('Display','off','MaxIter',10);
[indx, centroids_tmp, sumds] = kmeans(Xproj, h, 'Options', opts, 'EmptyAction', 'singleton');
centroids_tmp=centroids_tmp';
% [centroids_tmp, sumds, indx] =  yael_kmeans(Xproj', h, 'niter', 10, 'verbose', 0);
    
    H=(zeros(n,h));
    H(sub2ind([n,h],(1:n)', indx))=1;
    for j=1:size(H,2)
        if(numel(find(H(:,j))>0))
            H(:,j)=H(:,j)./norm(H(:,j),2);
        end
    end

C=centroids_tmp;
I=indx;
sumds_best=sum(sumds);
maxitr=5;
for i=1:maxitr
        [indx, centroids_tmp sumds] = kmeans(Xproj, h, 'Options', opts, 'EmptyAction', 'singleton');
         centroids_tmp=centroids_tmp';
    % [centroids_tmp, sumds, indx] = yael_kmeans(Xproj', h, 'niter', 10, 'verbose', 0);
    
    H_this=(zeros(n,h));
    H_this(sub2ind([n,h],(1:n)', indx))=1;
    for j=1:size(H_this,2)
        if(numel(find(H_this(:,j))>0))
            H_this(:,j)=H_this(:,j)./norm(H_this(:,j),2);
        end
    end
    %H_this=converttoH(indx,n,h);
    o=trace(H_this'*X*R*R'*X'*H_this);
    if(o > sumds_best)
        H=H_this;
        I=indx;
        C=centroids_tmp;
    end
end


