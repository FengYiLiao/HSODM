clc;clear;
addpath("package\");
%set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
      %'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};
set = {'mcp100.mat'};
%saveroot = "result\hsodm"; 
dataroot = "data\";
%load(dataroot+set{7});
load(dataroot+"G1.mat");
%load(dataroot+"n1000r20");

para.epislon   = 10^-6; %desired gradient accuracy
para.Maxiter   = 1000; %Maximum iterations
para.beta      = 0.5;   %line search parameter: reduction
para.gamma     = 2;     %line search parameter: a constant
para.Threshold = 2;     %This is cap delta (trigangle) in the paper
para.nu        = 0; %%
para.delta     = 0.01;     %the button right constant (control eigenvalue)
para.eta       = 1;     %initial line search step size
para.step      = 10;
para.adp_delta = false;  %adaptively tuning delta or not
para.linesearch = false;
para.L         = 2;   %adaptive parameter
para.delta_min = 10^-3;



prob.n         = height(C);
prob.rank      = 15;
prob.numvar    = height(C)*prob.rank; %dimension of total variabls;
prob.M         = obliquefactory(prob.rank,prob.n,true); %Create a mainfold
prob.cost      = @(X) cost(X,C);
prob.egrad     = @(X) 2*C*X;                            %euclidean gradient
prob.ehess     = @(X, U) 2*C*U;                         %euclidean hessian
prob.routine   = @routine;                              %power method routine
prob.vec2mani  = @vec2mani;


tic;
Out = HSODM(prob,para);  %main function 
%[x, xcost, info, options] = trustregions(prob); %manopt function
%[x, xcost, info, options] = steepestdescent(prob);
toc;



function y = routine(x,Xk,prob,para) %power method routine  %Xk is current iterate
    %gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    %x1       = reshape(x(1:end-1),prob.n,prob.rank);

    gk = getGradient(prob, Xk);
    [x1,x2]  = prob.vec2mani(x,prob);
    %x2       = x(end);   
%     y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
%         gk(:).'*x1(:) - para.delta*x2];
    Hk = getHessian(prob, Xk,x1);
    y = [reshape(Hk+x2*gk,[],1);
        gk(:).'*x1(:) - para.delta*x2];
end

function y = cost(X,C)
%     R = X*X';
%     y = C(:).'*R(:);
    R = C*X;
    y = X(:).'*R(:);
end

%this will depend on the specific manifold structure
function [vk,tk] = vec2mani(v,prob)
    %This part convert a long vector output by power method procedure back
    %to the structure in tangent space and a scaler
    vk      = reshape(v(1:end-1),prob.n,prob.rank);%vector on the tangent space
    tk      = v(end);                              %the last element is scalar t
end



% function y = routine_2(x,Xk,prob,para,delta) %power method routine  %Xk is current iterate
%     gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
%     x1       = reshape(x(1:end-1),prob.n,prob.rank);
%     x2       = x(end);   
%     y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
%         gk(:).'*x1(:) - delta*x2];
% end
