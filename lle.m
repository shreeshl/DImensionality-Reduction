function [Y] = llee(I,k)
tic
tol=1;
Idx= knnsearch(I', I', 'k', k);
n= size(I,2);
w = zeros(k,1);
W = zeros(n,n);
for i=1:n
	G = I(:,i)*ones(1,k) - I(:, Idx(i,:));
	Gi = G'*G;
	invg = inv(Gi + eye(k,k)*200);
	wi = invg*ones(k,1)/sum(sum(invg));
	for j=1:k
		W(i,Idx(i,j)) = wi(j);
	end
end


M = (eye(n) - W)'*(eye(n) - W);
[V,D] = eig(M); 
e=sort(diag(D));

i=3; % Take top 2 for visualization
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
Y = Y(:,2:end)';
toc
end

% for plotting
% idx = [1,50,100,150,200,250,300,350];
% plot(Y(1,idx),Y(2,idx),'-r')
% hold on
% for ii = 1:length(idx)
% text(Y(1,idx(ii)),Y(2,idx(ii)),num2str(ii),'Color','r')
% end