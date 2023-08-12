clc;clear;
set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
      'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};

rank =15;
epislon = 10^-6;
dataroot = "data\sdplib\";



%manopt

for i = 1:length(set)
    load(dataroot+set{i});
    name = split(set{i},'.');
    n = height(C);
    egrad = @(X) 2*C*X;
    ehess = @(X, U) 2*U;
    problem.M = obliquefactory(rank,n,true);
    problem.egrad = egrad;
    problem.ehess = ehess ;
    problem.cost  = @(X) trace(X.'*C*X);
    opt.maxiter = 30000;
    opt.tolgradnorm = epislon;
    [x, xcost, info, options] = trustregions(problem,[],opt);
    %disp('wait');
   % save("result\"+name{1}+"-result.mat",'x', 'xcost', 'info', 'options');
end
