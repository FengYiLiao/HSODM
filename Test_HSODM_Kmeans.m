clc;clear; close all;
addpath("package\");
%addpath("C:\Users\SOC-LAB\FengYi\Solver\DRSOM\Spherical_alg\Kmeans");
addpath("package\Spherical_alg\Kmeans");
%saveroot = "result\hsodm"; 
dataroot = "data\Kmeans\";
load(dataroot+"ecoli.mat");
%L = full(L);

% %testing
% n = 100;
% e = ones(n,1);
% L = spdiags([-e,2*e,-e],[-1,0,1],n,n);
% L = full(L);
% r = 10;

para.epislon   = 10^-6; %desired gradient accuracy
para.Maxiter   = 5000;  %Maximum iterations
para.beta      = 0.75;  %line search parameter: reduction
para.gamma     = 2;     %line search parameter: a constant
para.Threshold = 2;     %This is cap delta (trigangle) in the paper
para.nu        = 0.45;
para.delta     = 2;     %the button right constant (control eigenvalue)
para.eta       = 1;     %initial line search step size
para.step      = 1;
para.adp_delta = true;  %adaptively tuning delta or not



prob.n         = height(M);
prob.rank      = 10;%width(M);
Y              = M./max(abs(M));
W              = -Y*Y'; 
lambda         = 100;
prob.M         = Kmeansfactory(prob.n,prob.rank); %Create a mainfold


%prob.M.retraction()

prob.cost      = @(X) cost(X,W,lambda);

% X0             = randn(prob.n,prob.rank);
% prob.cost(X0)
%tang           = prob.M.proj(X0,X2);

%XX             = prob.M.retraction(X0,tang);
%XX'*XX
prob.egrad     = @(X) egrad(X,W,lambda);      %euclidean gradient
%prob.ehess     = @(X, U) ehess(X,L,invL,U); %euclidean hessian
prob.routine   = @routine;                  %power method routine

%% Random initialization
X0             = randn(prob.n,prob.rank);
X0             = Kmeans_retrac(X0,0,0);
para.X0 = X0;

%prob.M.retr    = @Kmeans_retrac;%(X,U,t)

tic;
opt.maxiter = 30000;
opt.tolgradnorm = para.epislon;

%[x, xcost, info, options] = trustregions(prob,[],opt); %manopt function
% [x, xcost, info, options] = steepestdescent(prob,X0,opt);
% [x, xcost, info, options] = barzilaiborwein(prob,X0,opt);
[x, xcost, info, options] = conjugategradient(prob,X0,opt);
%Out = HSODM(prob,para);  %main function 
toc;



function y = routine(X,Xk,prob,para) %power method routine  %Xk is current iterate
    gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    x1       = reshape(X(1:end-1),prob.n,prob.rank);
    x2       = X(end);
    y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
        gk(:).'*x1(:) - para.delta*x2];
end

function y = cost(X,W,lambda)
    %lambda is prechosen
    R           = W'*X;
    %d           = sum(X.*X,2);
%     y = 0.5*(X(:).'*R(:)) + (d'*invL*d)/4;
    PX          = X; %projected X 
    nonpos      = PX(PX < 0);
    y           = (X(:).'*R(:)) +lambda *norm(nonpos)^2;
    %R = X'*X;
    %y = L(:).'*R(:);
end

function y = egrad(X,W,lambda)
    %d = sum(X.*X,2);
    %R = invL*d;
    %R = L\d;
    nonpos                = X;
    nonpos(nonpos>=0)     = 0;
    y                     = 2*W*X + 2*lambda*nonpos;
end

function y = ehess(X,L,invL,V)
    d  = sum(X.*X,2);
    d2 = sum(X.*V,2);
%     invLd  = invL*d;
%     invLd2 = invL*(2*d2);
    invLd  = L\d;
    invLd2 = L\(2*d2);
    y      = L*V+invLd.*V+invLd2.*X;
end