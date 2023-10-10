clc;clear;
set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
      'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};

rank =15;
epislon = 10^-5; %desired arrurarcy
dataroot = "data\sdplib\";
saveroot = "result\manopt\"; 

%Use manopt to solve the problem
%Remember to include manopt's path

for i = 1:length(set)
    load(dataroot+set{i});
    %load("data\G1.mat");
    name = split(set{i},'.');
    n = height(C);
    problem.M = obliquefactory(rank,n,true); %create a manifold
    problem.egrad = @(X) 2*C*X; %euclidean gradient
    problem.ehess = @(X, U) 2*C*U;%euclidean hessian
    problem.cost  = @(X) cost(X,C);%trace(X.'*C*X);
    opt.maxiter   = 30000;
    opt.tolgradnorm = epislon;
    [x, xcost, info, options] = trustregions(problem,[],opt); %manopt function
    save(saveroot+name{1}+"-result.mat",'x', 'xcost', 'info', 'options');
end

function y = cost(X,C)
    R = X*X';
    y = C(:).'*R(:);
end
