function [S X finalobj]=optimize_subspaces_blockcoord(C,d,seq,iniC,prevX)
m=size(C,3);
n=size(C,1);

%parameters

ratio=(n/m);

if(numel(prevX)==0)
    X=initializeX(n,m);
  % X=iniX;
else
    X=prevX;
  %  X=initializeX(n,m);
   X=initializeX(n,m,iniC,C,prevX); 
end
if(nargin<3)
    seq=1;%sequenial selection
end
I=eye(n,n);
obj=0;
constraints=[];

T=5;%20number of iterations
if(nargin<6)
sdp=0; %% eigen apporach
%sdp=1;%% solve using sdp for 2 case
%sdp=2;%solve using quadprog
%sdp=3;%solve using lagrangian
end
O=[];
Xsave=X;
index_list_orig=sort(combnk(1:m,2),1);

index_list=index_list_orig;
objbef=0;
for i=1:m
        objbef=objbef+trace(C(:,:,i)*X(:,:,i));
end


sdp=0;

if(m==2) %%2 case
    
    Cij=C(:,:,1)-C(:,:,2);
    
    
    %[Xret obj]=solve2case_sdp(Cij,eye(n),n,m,ratio);
    [Xret obj]=solve2case(Cij,n,m,ratio);
    X(:,:,1)=Xret;
    
    X(:,:,2)=I-Xret;
    
else %% more than 2 case
    
    for itr=1:T
        if(itr>1)
            prand_index=randperm(size(index_list,1));
            if(seq==0)%%random incdex selection
                index_list=index_list_orig(prand_index,:);
            end
        end
        for p=1:size(index_list,1)
            
            ind=index_list(p,:);
            i=ind(1);
            j=ind(2);
            if(i~=j)
                
                if (sdp==2)
                    Cij=C(:,:,i)-C(:,:,j);
                    Xret=solve2case_sdp_relax(Cij,eye(n),n,m,ratio);
                    
                    Xret_new=(Xret+Xret')/2;
                    X(:,:,i)=(Xret_new);
                    
                    X(:,:,j)= eye(n) - Xret_new;
                end
                if(sdp==0)
                    otherind=setdiff(1:m,ind);
                    if(numel(otherind)>1)
                        E = I - sum(X(:,:,otherind),3);
                        L=sum(X(:,:,otherind),3);
                    else
                        
                        E = I - X(:,:,otherind);
                        L=X(:,:,otherind);
                    end
                    Cij=(C(:,:,i)-C(:,:,j));
                    
                
                    %[ES EU]=mex_dsyev_k(E,2*ratio);%%can use this also
                    [EU ES]=eigs(E,2*ratio,'lm');
                    EU=EU(:,1:2*ratio);
                    
                    [Xret obj]=solve2case(EU'*Cij*EU,n,m,2*ratio);
                    
                    
                    Xret_new=EU*Xret*EU';
                    Xother=E - Xret_new;
                    Xret_new=find_rankapporx(Xret_new,ratio);
                    
                    Xret_new=(Xret_new+Xret_new')/2;
                    
                    
                    X(:,:,i)=(Xret_new);
                    Xother=find_rankapporx(Xother,ratio);
                    
                    X(:,:,j)=  (Xother+Xother')/2;
                end
                
       
            end%end if
            
        end%end for
        
        
    
    obj=0;
    for i=1:m
        obj=obj+trace(C(:,:,i)*X(:,:,i));
    end
    O=[O; [obj]];
    
    if(itr==1)
        if(O(itr)<objbef)
            X=Xsave;
            break;
        else
            Xsave=X;
       end
    end
    
    if(itr>1)
        if((abs(O(itr)-O(itr-1))<=1) ||(O(itr) < O(itr-1)))
            break;
        else
            Xsave=X;
        end
    
    end
    
    end
end
if(m>2)
    X=Xsave;
end

   
S=[];
%
finalobj=0;
for i=1:m
    finalobj=finalobj+trace(C(:,:,i)*X(:,:,i));
    X1=double(X(:,:,i));
    
    [Ux1 Sx1]=eigs(X1,ratio);
    
    S=cat(3,S,Ux1(:,1:ratio));
end

end



function [X obj]=solve2case(C,n,m,rhat)
   
    C=(C+C')/2;
  
   %[S U]=mex_dsyev_k(C,rhat/2);%rsvd
    [U S V]=rsvd(C,rhat/2);%rsvd use this or the prev line
    if(size(S,1)<rhat/2)
        [U S V]=svdsecon(C,rhat/2);
    end
%    
%    [U S]=svd(C,'econ');
    %[U S]=eigs(C,rhat,'lm');
    [l ind]=sort(diag(S),'descend');
   % size(l)
    obj=sum(l(1:rhat/2));
    X=U(:,ind(1:rhat/2))*U(:,ind(1:rhat/2))';
end
 
function C=clean(C)

    [U S]=eigs(C);

    ind=find(diag(S)<=0.01);

    for i=1:length(ind)
        S(ind(i),ind(i))=0;
    end

    C=real(real(U)*real(S)*real(U)');
    C=(C+C')/2;
end
 
 


function Etemp=find_rankapporx(E,r)
    E=(E+E')/2;
    [U1 S1]=svd(E);
    Etemp=0;
   
    for i=1:r

        Etemp=Etemp+real(U1(:,i)*S1(i,i)*U1(:,i)');
    end

end
    

function finalobj=evaluate_obj(C,X)
m=size(C,3);
finalobj=0;
for i=1:m
    finalobj=finalobj+trace(C(:,:,i)*X(:,:,i));
end
end