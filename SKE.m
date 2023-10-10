clc;clear;
addpath("package\");
dataroot = "data\KSE\";

%saveroot = "result\hsodm\SKE\"; 
saveroot = "result\manopt\SKE\";

para.epislon   = 10^-6;%desired gradient accuracy
para.Maxiter   = 5000;%Maximum iterations
para.beta      = 0.5;%line search parameter: reduction
para.gamma     = 2;%line search parameter: a constant
para.Threshold = 2; %This is cap delta (trigangle) in the paper
para.nu        = 0.45;
para.delta     = 2; %the button right constant (control eigenvalue)
para.eta       = 1; %%initial line search step size
para.step      = 5;



for N = [1000,2000,5000,7000,10000]
    for j = [20,50]
        name = "n"+num2str(N)+"r"+num2str(j);
        load(dataroot+name+".mat");
        
        
        prob.n         = n;%height(C);
        prob.rank      = r;
        prob.M         = stiefelfactory(prob.n,prob.rank,1); %Create a mainfold
        invL           = inv(L);
        prob.cost      = @(X) cost(X,L,invL);
        prob.egrad     = @(X) egrad(X,L,invL);      %euclidean gradient
        prob.ehess     = @(X, U) ehess(X,L,invL,U); %euclidean hessian
        prob.routine   = @routine;                  %power method routine
        
        opt.maxiter = 30000;
        opt.tolgradnorm = para.epislon;
        [x, xcost, info, options] = trustregions(prob,[],opt); %manopt function
        
        %Out = HSODM(prob,para); 
        %save(saveroot+name+"-result.mat",'Out');
        save(saveroot+name+"-result.mat",'x', 'xcost', 'info', 'options');
    end
end



function y = routine(X,Xk,prob,para) %power method routine  %Xk is current iterate
    gk = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    x1 = reshape(X(1:end-1),prob.n,prob.rank);
    x2 = X(end);
    y  = [reshape(prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1)+para.delta*gk,[],1);
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