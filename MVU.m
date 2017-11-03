function [Y] = MVU(I,k)

n= size(I,2);
Idx= knnsearch(I', I', 'k', k);

cvx_begin
variable K(n,n) semidefinite
maximize trace(K)
subject to
  for i=1:n 
    for j=1:n
        if (ismember(j, Idx(i,:)))
            K(i,i)-2*K(i,j)+K(j,j) == norm(double(I(:,i)-I(:,j)))^2;
        end
    end
  end
  sum(sum(K)) == 0
cvx_end

[V,D] = eig(K); 
e=sort(diag(D), 'descend');

i=2; % Take top 2 for visualization
e= e(1:i); 
evectors= []; % dominant eigen vectors
evalues = [];
for i=1:n
    if (ismember(D(i,i),e))
        evectors = [evectors V(:,i)];
        evalues  = [evalues D(i,i)];               
    end
end

Y= evectors*sqrt(diag(evalues)); 
end
     
  
