clc;clear; close all;
addpath("package\");


%saveroot = "result\hsodm"; 
dataroot = "data\KSE\";
load(dataroot+"n1000r20"+".mat");

para.epislon   = 10^-6; %desired gradient accuracy
para.Maxiter   = 1000; %Maximum iterations
para.beta      = 0.5;   %line search parameter: reduction
para.gamma     = 2;     %line search parameter: a constant
para.Threshold = 2;     %This is cap delta (trigangle) in the paper
para.nu        = 0.45;
para.delta     = 2;     %the button right constant (control eigenvalue)
para.eta       = 1;     %initial line search step size

prob.n         = n;%height(C);
prob.rank      = r;
prob.M         = obliquefactory(prob.rank,prob.n,true); %Create a mainfold
prob.cost      = @(X) cost(X,C);
prob.egrad     = @(X) 2*C*X;                            %euclidean gradient
prob.ehess     = @(X, U) 2*U;                           %euclidean hessian
prob.routine   = @routine;                              %power method routine

tic;
Out = HSODM(prob,para);  %main function 
toc;



function y = routine(x,Xk,prob,para) %power method routine  %Xk is current iterate
    gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    x1       = reshape(x(1:end-1),prob.n,prob.rank);
    x2       = x(end);   
%     y        = zeros(prob.n*prob.rank+1,1);
%     temp     = prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+para.delta*gk;
%     y(1:end-1) = temp(:);
%     y(end)   = gk(:).'*x1(:) - para.delta*x2;
    y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+para.delta*gk,[],1);
        gk(:).'*x1(:) - para.delta*x2];
end

function y = cost(X,C)
    R = X*X';
    y = C(:).'*R(:);
end