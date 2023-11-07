clc;clear; close all;
addpath("package\");
addpath("package\Spherical_alg\SNL\")

%saveroot = "result\hsodm"; 
%dataroot = "data\KSE\";
%load(dataroot+"n2000r20"+".mat");
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
% prob.n         = n;%height(C);
% prob.rank      = r;

% prob.cost      = @(X) cost(X,L,invL);
% prob.egrad     = @(X) egrad(X,L,invL);      %euclidean gradient
% prob.ehess     = @(X, U) ehess(X,L,invL,U); %euclidean hessian
% prob.routine   = @routine;                  %power method routine



SDPinitialize = 1;
n = 100;
r = 3; % dimension
noise = 0.2; % [0,1]
lamreg = 0.1/n; % [0,1];
ax = 0.3; % axis of anchor;
ad = 15; % average degree
taus = (3*ad/(4*pi*n))^(1/3); % distance threshod between sensor
taua = sqrt(2)*ax*1.1; % [0,ax*sqrt(2)] distance threshod between sensor and anchor
[Na,Nx,An,Sen,dx2,da2] = SNL_data(n,r,ax,taus,taua,noise);
par.maxiter = 10000;
%par.verbose = verbose;
%par.tol = tol;
M         = SNLfactory(Sen,An,Nx,Na,dx2,da2); %Create a mainfold

disp('wait');


% [R,iter,fval,gradnorm,pfeas,ttime] = Rie_SNL(R0,Sen,An,Nx,Na,dx2,da2,lamreg,par);
% RMSD = norm(R-Sen,'fro')/sqrt(n);
% if (r == 3)&&(plotyes)
%     figure(1);
%     plot3(Sen(:,1),Sen(:,2),Sen(:,3),'r*');
%     hold on;
%     plot3(An(:,1),An(:,2),An(:,3),'kp');
%     plot3(R(:,1),R(:,2),R(:,3),'b+');
%     plot3([Sen(:,1),R(:,1)]',[Sen(:,2),R(:,2)]',[Sen(:,3),R(:,3)]','g-');
%     title(['Rie  ','RMSD=',num2str(RMSD)]);
% end
% if (useSQP==0)&&(useIPM==0)
%     ttime = ttime+ttime1;
% end
% fprintf('\n $n=%2d$& RDRSOM & %6.7e & %3.2e & %3.2e & %3.2e & %3.2e \\\\',n,fval,gradnorm,pfeas,RMSD,ttime);


%% Random initialization
R = randn(n,r);
AA = R'*R;
[P,V] = eig(full(AA));
d = diag(V);
d = sqrt(1./d);
V = spdiags(d,0,r,r);
X0 = R*(P*V*P');
para.X0 = X0;

tic;
opt.maxiter = 30000;
opt.tolgradnorm = para.epislon;

%[x, xcost, info, options] = trustregions(prob,[],opt); %manopt function
%[x, xcost, info, options] = steepestdescent(prob,[],opt);
Out = HSODM(prob,para);  %main function 
toc;



function y = routine(X,Xk,prob,para) %power method routine  %Xk is current iterate
    gk       = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    x1       = reshape(X(1:end-1),prob.n,prob.rank);
    x2       = X(end);
    y = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+x2*gk,[],1);
        gk(:).'*x1(:) - para.delta*x2];
end

function y = cost(X,L,invL)
    R = L*X;
    d = sum(X.*X,2);
%     y = 0.5*(X(:).'*R(:)) + (d'*invL*d)/4;
    y = 0.5*(X(:).'*R(:)) + (d'*(L\d))/4;
    %R = X'*X;
    %y = L(:).'*R(:);
end

function y = egrad(X,L,invL)
    d = sum(X.*X,2);
    %R = invL*d;
    R = L\d;
    y = L*X + R.*X;
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