function X=initialize(n,m,iniC,C,prevX)
I=eye(n);
%[U S]=eigs(I,n);
U=eye(n);
r=round(n/m);
X=zeros(n,n,m);

for i=1:m
    ind=r*(i-1)+1:r*i;
    
    X(:,:,i)=U(:,ind)*U(:,ind)';
end

if(nargin>=3)
    X1=X;
    X2=prevX;
    
    [U S]=svd((iniC+iniC')/2);
    r=round(n/m);
    X=zeros(n,n,m);
    
    for i=1:m
        ind=r*(i-1)+1:r*i;
        
        X(:,:,i)=U(:,ind)*U(:,ind)';
     
    end
    obj1=0;
    obj2=0;
    obj3=0;
  
    for i=1:m
        obj1=obj1+trace(X1(:,:,i)*C(:,:,i)*X1(:,:,i)');
        obj2=obj1+trace(X(:,:,i)*C(:,:,i)*X(:,:,i)');
        obj3=obj3+trace(X2(:,:,i)*C(:,:,i)*X2(:,:,i)');
    end
    if(obj1 >obj2)
        X=X1;
    end
   
    if((obj3 > obj2) ||(obj3 > obj1))
        X=X2;
    end
end
