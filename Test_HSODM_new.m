clc;clear;
%This version has more update!
addpath("package\");
addpath(genpath("data\"));
%cd data\
% set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
%       'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};

idx   = 2;
switch idx 
    case 1
       prob = dominant_invariant_subspace_problem([],  512, 12);
       prob.n         = 512;
       prob.rank      = 12;
       prob.Tolvar    = prob.n*prob.rank;
       prob.vec2mani  = @vec2mani;
    case 2
       %168,240,20
       prob             = truncated_svd_problem([], 168, 240, 20);
       prob.vec2mani    = @(v,prob)  vec2mani(v,prob,168, 240, 20);
       prob.mani2vec    = @mani2vec;
       prob.Tolvar      = 168*20+240*20;
       prob.nummanifold = 2;
       prob.dim{1}      = [168;20];
       prob.dim{2}      = [240;20];
       prob.numvar{1}   = prob.dim{1}(1)*prob.dim{1}(2);
       prob.numvar{2}   = prob.dim{2}(1)*prob.dim{2}(2);
%       para.X0        = prob.x0;
    case 3
       prob             = lrmc_grassmann(2000, 5000, 10, 4);
       prob.n           = 2000;
       prob.rank        = 10;
       prob.Tolvar      = prob.n*prob.rank;
       prob.vec2mani    = @vec2mani;
    case 4
       prob             = maxcut(22);
       prob.n           = 2000;
       prob.rank        = 64;
       prob.Tolvar      = prob.n*prob.rank;
       prob.vec2mani    = @vec2mani;
    case 5
       prob             = rotation_synchronization(3, 50, .75);
       prob.vec2mani    = @(v,prob) vec2mani(v,prob,50);
       prob.Tolvar      = 9*50;
       prob.nummanifold = 1;
       %prob.X0          = prob.x0;
    case 6
       prob             = shapefit_leastsquares(500, 3);
       prob.n           = 3;
       prob.rank        = 500;
       prob.Tolvar      = prob.n*prob.rank;
       prob.vec2mani    = @vec2mani;
end



%prob = dominant_invariant_subspace_problem([],  512, 12);
%prob = truncated_svd_problem([], 168, 240, 20)
%prob = lrmc_grassmann(2000, 5000, 10, 4);
%prob =  maxcut(22);
%prob  = rotation_synchronization(3, 50, .75);


% prob.n         = 512;
% prob.rank      = 12;



saveroot = "result\hsodm"; 
dataroot = "data\";
%load(dataroot+set{7});
%load(dataroot+"G1.mat");
%load(dataroot+"n1000r20");

para.epislon    = 10^-6; %desired gradient accuracy
para.Maxiter    = 1000; %Maximum iterations
para.beta       = 0.5;   %line search parameter: reduction
para.gamma      = 1;     %line search parameter: a constant
para.Threshold  = 2;     %This is cap delta (trigangle) in the paper (Dead)
para.nu         = 0; %sqrt(para.epislon)
para.delta      = 0.1;%sqrt(para.epislon);     %the button right constant (control eigenvalue)
para.eta        = 1;     %initial line search step size
para.step       = 1;     %How often to print
para.adp_delta  = false;  %adaptively tuning delta or not
para.linesearch = false;%true;
para.L          = 2;   %adaptive parameter
para.delta_min  = 10^-3; %dead
para.ck         = 1;    
para.rho2       = 2/3;
para.rho1       = 1/3;


% prob.M         = obliquefactory(prob.rank,prob.n,true); %Create a mainfold
% prob.cost      = @(X) cost(X,C);
% prob.egrad     = @(X) 2*C*X;                            %euclidean gradient
% prob.ehess     = @(X, U) 2*C*U;                         %euclidean hessian
%prob.Tolvar    = prob.n*prob.rank;
prob.routine   = @routine;                              %power method routine



%tic;
Out = HSODM(prob,para);  %main function 
%[x, xcost, info, options] = trustregions(prob); %manopt function
%[x, xcost, info, options] = steepestdescent(prob);
%toc;
%save(saveroot +'\Result_HSODM_Maxcut',"Out");


% function y = routine(x,Xk,prob,para) %power method routine  %Xk is current iterate
%     gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
%     x1       = reshape(x(1:end-1),prob.n,prob.rank);
%     x2       = x(end);   
%     y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
%         gk(:).'*x1(:) - para.delta*x2];
% end


function y = routine(x,Xk,prob,para) %power method routine  %Xk is current iterate
    %gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    %x1       = reshape(x(1:end-1),prob.n,prob.rank);

    %storedb = StoreDB(2);
    %key = storedb.getNewKey();
    

    %gk = getGradient(prob, Xk);
    
    newkey = prob.storedb.getNewKey();
    gk = getGradient(prob, Xk,prob.storedb,newkey);

    [v,t]  = prob.vec2mani(x,prob);
    %x2       = x(end);   
%     y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
%         gk(:).'*x1(:) - para.delta*x2];

    %Hk = getHessian(prob, Xk,x1);
    Hk = getHessian(prob, Xk,v,prob.storedb,newkey);
    vec_Hk = prob.mani2vec(Hk,prob);
    vec_gk = prob.mani2vec(gk,prob);
    %Hk = getHessianFD(prob, Xk,x1);
    if para.ck ~= 1
%         y = [reshape(Hk+t*gk/para.ck,[],1);
%                 gk(:).'*v(:)/para.ck - para.delta*t/(para.ck^2)];
        y = [vec_Hk+t*vec_gk/para.ck;
            prob.M.inner(Xk,gk,v)/para.ck - para.delta*t/(para.ck^2)];
    else
%         y = [reshape(Hk+t*gk,[],1);
%             gk(:).'*v(:) - para.delta*t];
        
        y = [vec_Hk+t*vec_gk;
            prob.M.inner(Xk,gk,v) - para.delta*t];
    end
end



%this will depend on the specific manifold structure
function [vk,tk] = vec2mani(v,prob,m,n,p)
    %This part convert a long vector output by power method procedure back
    %to the structure in tangent space and a scaler
    if nargin == 2
        vk      = reshape(v(1:end-1),prob.n,prob.rank);%vector on the tangent space
        tk      = v(end);                              %the last element is scalar t
    elseif nargin == 3
        vk      = zeros(3,3,50);
        for i = 1:m
            vk(:,:,i) = reshape(v((i-1)*9+1:i*9),3,3);
        end
        %vk      = reshape(v(1:end-1),prob.n,prob.rank);%vector on the tangent space
        tk      = v(end);                              %the last element is scalar t
    elseif nargin == 5
        vk.U    = reshape(v(1:m*p),m,p);
        vk.V    = reshape(v(m*p+1:end-1),n,p);
        tk      = v(end);
    end
end



function v = mani2vec(gk,prob)
    %Reture the gradient into vector form
    %This part is needed for the power method
    
    v     = zeros(prob.Tolvar,1);
    start = 1;
    v(start:start+prob.numvar{1}-1) = reshape(gk.U,[],1);
    start = start + prob.numvar{1};
    v(start:start+prob.numvar{2}-1) = reshape(gk.V,[],1);
end




% function y = cost(X,C)
% %     R = X*X';
% %     y = C(:).'*R(:);
%     R = C*X;
%     y = X(:).'*R(:);
% end

% function y = routine_2(x,Xk,prob,para,delta) %power method routine  %Xk is current iterate
%     gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
%     x1       = reshape(x(1:end-1),prob.n,prob.rank);
%     x2       = x(end);   
%     y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
%         gk(:).'*x1(:) - delta*x2];
% end
