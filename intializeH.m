function [Hint R len0 len1]=intializeH(X,p,h,m)

[n p]=size(X);
d=p/m;
len0=0;
len1=0;
%R = eigenvalue_allocation(X, m);
[R len0 len1]=callCKMeansIni(X,p,h,m);

Xproj = X*R; % pre-projecting X
[C, sumds, idxs_best]=yael_kmeans(Xproj', h, 'niter', 10, 'verbose', 0);
sumds_best=sum(sumds);
for i = 1:1%m %earlier
    Xsub = Xproj(:, (1:d) + (i-1)*d);
    %opts = statset('Display','off','MaxIter',10);
    %[idxs, centers, sumds] = kmeans(Xsub, h, 'Options', opts, 'EmptyAction', 'singleton');
    [centroids_tmp, sumds, idxs] = yael_kmeans (Xsub', h, 'niter', 10, 'verbose', 0);%%iter =100
    
    if(sumds_best>=sum(sumds))
        idxs_best=idxs;
        sumds_best=sum(sumds);
     
    end
    
end
H=(zeros(n,h));
H(sub2ind([n,h],(1:n)', idxs_best))=1;
for j=1:size(H,2)
    if(numel(find(H(:,j))>0))
        H(:,j)=H(:,j)./norm(H(:,j),2);
    end
end

Hint=H;
end



function [R len0 len1]=callCKMeansIni(X,p,h,m)
X=X';
n=size(X,2)
p = size(X, 1);
R = eye(p, p, 'single');
RX = R'*X;
len = ones(m, 1) * floor(p / m);
len(1:mod(p, m)) = len(1:mod(p, m)) + 1;  % p = m * floor(p / m) + mod(p, m)


len0 = 1 + cumsum([0; len(1:end-1)]);
len1 = cumsum(len);
if (length(h) == 1)
  h1 = ones(m, 1) * h;
end
for (i=1:m)
    perm = 1:n;%randperm(n);
    perm=perm(1:h1(i));
    D{i} = RX(len0(i):len1(i), perm);
end
B = zeros(n, m, 'int32');
for (i=1:m)
    B(:, i) = euc_nn_mex(D{i}, RX(len0(i):len1(i), :)); %D{i} is the centers, RX are the original dataset projected on to this subspace, B(:,i) is finding nn for each point in Rx to D
    
    DB(len0(i):len1(i), :) = D{i}(:, B(:, i));% heare each point in Rx is replaced with its center
    
end

    
    tmp = R * DB;
    tmp = tmp - X;
    tmp = tmp.^2;
    obj = mean(sum(tmp, 'double'));
    
    
    % update R
    [U, S, V] = svd(X * DB', 0);
    
    R = U * V';
    
    
    % update R*X
    RX = R' * X;  % return in a structure % R is the projection matrix of all the subspace (concatanted).
    
    

X=X';
end

    function H=converttoH(idxs,n,h)
        H=zeros(n,h);
        for i=1:length(idxs)
            H(i,idxs(i))=1;
        end
        
        for j=1:size(H,2)
            if(norm(H(:,j),2)>0)
                H(:,j)=H(:,j)./norm(H(:,j),2);
            end
        end
    end