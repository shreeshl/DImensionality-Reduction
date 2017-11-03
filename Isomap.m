function [Y] = Isomap(I,k)
$tic
n=size(I,2);
Idx= knnsearch(I', I', 'k', k);

G= repmat(realmax,n,n);

for i=1:size(Idx,1)
    for j=1:size(Idx,2)
        G(i,Idx(i,j)) = norm(I(:,i) - I(:,Idx(i,j)))^2;
    end
end


for k=1:n
    for i=1:n
        for j=1:n
            G(i,j) = min(G(i,j), G(i,k) + G(k,j));
        end
    end
end

% Formulate tau(G) = -HSH/2 

S= G.^2; 
H = eye(n,n) + repmat(-1/n, n, n);
tau= -H*S*H/2 ; 

[V,D] = eig(tau); 
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

Y= evectors*sqrt(diag(evalues))'; 
%toc
end




