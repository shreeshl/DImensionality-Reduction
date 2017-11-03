function [Y] = fastmvu(I,k,r,m)
tic
n= size(I,2);
Idx= knnsearch(I', I', 'k', k);

W = zeros(n,n);
for i=1:n
	G = I(:,i)*ones(1,k) - I(:, Idx(i,:));
	Gi = G'*G;
	invg = inv(Gi+eye(k,k)*150);
	wi = invg*ones(k,1)/sum(sum(invg));
	for j=1:k
		W(i,Idx(i,j)) = wi(j);
	end
end

phi = (eye(n) - W)'*(eye(n) - W);
Q = [eye(m);inv(phi(m+1:end,m+1:end))*phi(m+1:end,1:m)];
Idx2 = Idx(:,1:r);
idx = linspace(1,m,m);
left_index = linspace(m+1,n,n-m);
iterate = 1;
while(iterate==1)
	cvx_begin
		variable L(m,m) semidefinite
		M = Q*L*Q';
		maximize trace(M)
		subject to
		  for i=idx 
		    for j=idx
		        if (ismember(j, Idx2(i,:)))
		            M(i,i)-2*M(i,j)+M(j,j) <= norm(double(I(:,i)-I(:,j)))^2;
		        end
		    end
		  end
		  sum(sum(M)) == 0
	cvx_end

	curr_size = size(idx,2);
	for i=left_index 
	    for j=left_index
	    	if (ismember(j, Idx2(i,:))) && ~(M(i,i)-2*M(i,j)+M(j,j) <= norm(double(I(:,i)-I(:,j)))^2);
	    		if i==j and ~(ismember(i, idx))
	    			idx = [idx i];
	    			left_index(left_index==i) = [];
	    		
	    		elseif ~(ismember(j, idx))
	    			idx = [idx j];
	    			left_index(left_index==j) = [];
	    		
	    		elseif ~(ismember(i, idx))
	    			idx = [idx i];
	    			left_index(left_index==i) = [];
	    		end
	    		iterate = 1;
	    		%break
	    	end
	    end
	end
	
	if curr_size==size(idx,2)
		break
	end

[V,D] = eig(L); 
e=sort(diag(D), 'descend');

i=2; % Take top 2 for visualization
e= e(1:i); 
evectors= []; % dominant eigen vectors
evalues = [];
for i=1:m
    if (ismember(D(i,i),e))
        evectors = [evectors V(:,i)];
        evalues  = [evalues D(i,i)];               
    end
end

Lnew= evectors*sqrt(diag(evalues)); 
Y = Q*Lnew;
Y = Y';
toc
end    
