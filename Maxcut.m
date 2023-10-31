clc;clear;
addpath("package\");
set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
      'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};

saveroot = "result\hsodm\MaxCut\"; 
dataroot = "data\sdplib\";
%load(dataroot+set{1});

para.epislon   = 10^-5;%desired gradient accuracy
para.Maxiter   = 30000;%Maximum iterations
para.beta      = 0.5;%line search parameter: reduction
para.gamma     = 2;%line search parameter: a constant
para.Threshold = 2; %This is cap delta (trigangle) in the paper
para.nu        = 0.45;
para.delta     = 2; %the button right constant (control eigenvalue)
para.eta       = 1; %%initial line search step size
para.step      = 10;
para.adp_delta = true;  %adaptively tuning delta or not
para.delta_min = 10^-3;



for i = 1:length(set)
    load(dataroot+set{i});
    name = split(set{i},'.');

    prob.n = height(C);
    prob.rank    = 15;
    prob.M       = obliquefactory(prob.rank,prob.n,true); %Create a mainfold
    prob.cost    = @(X) cost(X,C);
    prob.egrad   = @(X) 2*C*X;    %euclidean gradient
    prob.ehess   = @(X, U) 2*C*U; %euclidean hessian
    prob.routine = @routine;      %power method routine

    Out = HSODM(prob,para); 
    save(saveroot+name{1}+"-result.mat",'Out');
end


function y = routine(x,Xk,prob,para) %power method routine  %Xk is current iterate
    gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    x1       = reshape(x(1:end-1),prob.n,prob.rank);
    x2       = x(end);   
    y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
        gk(:).'*x1(:) - para.delta*x2];
end

function y = cost(X,C)
%     R = X*X';
%     y = C(:).'*R(:);
    R = C*X;
    y = X(:).'*R(:);
end

function y = routine_2(x,Xk,prob,para,delta) %power method routine  %Xk is current iterate
    gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    x1       = reshape(x(1:end-1),prob.n,prob.rank);
    x2       = x(end);   
    y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
        gk(:).'*x1(:) - delta*x2];
end
