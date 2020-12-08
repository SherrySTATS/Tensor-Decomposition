
%% calculate by the tensor toolbox
clear all;
% load the data for p1
load('aminoacid.mat')
p = parafac_als(X,3)
%score matrix
p.U{1}

%% ALS CP decomposition
clear all;
load('aminoacid.mat')
% initialization
s = size(X); 
R = 3; 
for i = 1:R
    P{i} = rand(s(i), R);
end
%iteration
N = 1000; 
X1 = tenmat(X, 1);
for i = 1:N  
    for k = 1:3
        V = (P{1}'*P{1}).*(P{2}'*P{2}).*(P{3}'*P{3});
        V = V./(P{k}'*P{k});
        khat{1} = khatrirao(P{3}, P{2});
        khat{2} = khatrirao(P{3}, P{1});
        khat{3} = khatrirao(P{2}, P{1});
        P{k} = double(tenmat(X, k) * khat{k} * inv(V' * V) * V');
        P{k} = normc(P{k});  
        lambda(k) = norm(P{k}); 
    end   
 % check for convergence
    U= ktensor(lambda',P);
    U1 = tenmat(U,1);
    dif = tensor(X1-U1);
    err(i) = log(innerprod(dif,dif));
end

%score matrix
U{1}

%% compare
figure
hold on
for i = 1:3
    plot (1:5,U{1}(:,i))
    xlabel('Sample idx')
    ylabel('Score')
end
hold off

%present
figure
hold on
for i = 1:3
    plot(U{2}(:,i))
end
title('Emission Loading')
hold off

figure
hold on
for i = 1:3
    plot(U{3}(:,i))
end
title('Excitation Loading')
hold off