clc;clear;
  
set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
      'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};
saveroot = "result\hsodm\eta10\"; 
dataroot = "data\sdplib\";
% load("mcp250-3.mat");
% [At,b,c,K]=fromsdpa('mcp500-4.dat-s');
% n = width(At);
% C = reshape(c,n,n);



% At = zeros(n,n^2);
% idx = find(diag(ones(n,1)));
% for i = 1:n
%     At(i,idx(i)) = 1;
% end
% b = ones(n,1);
% c = reshape(C,[],1);
% K.s = n;
% [res,info]=SolveMosek(At',b,c,K);
% Truecost = res.sol.itr.pobjval;

rank =15;%Factorization
para.epislon = 10^-6;%desired gradient accuracy
para.Maxiter = 30000;%Maximum iterations

para.beta = 0.5;%line search parameter: reduction
para.gamma = 2;%line search parameter: a constant
para.Threshold =2; %This is cap delta (trigangle) in the paper
para.nu = 0.45;
para.delta = 2; %the button right constant (control eigenvalue)
para.eta  = 10; %initial line search step size
for i = 1:length(set)
    load(dataroot+set{i});
    name = split(set{i},'.');
    prob.n = height(C);
    prob.rank = rank;
    prob.C = C; %problem data
    prob.M = obliquefactory(prob.rank,prob.n,true); %Create a mainfold
    
    X0 = prob.M.rand(); %initial point
    
    obj = [];
    Grad = [];
    egrad = @(X) 2*C*X; 
    ehess = @(X, U) 2*U;

    stop ="";
    fprintf("iter  |   obj  |   opt   |  rgrad  |   eta   \n");
    Xk = X0;
    for iter = 1:para.Maxiter
        
        gk = prob.M.egrad2rgrad(Xk,egrad(Xk)); %euclidean gradient to Riemannian gradient
        Afun = @(x) F(x,Xk,prob,para);%Define a routine
        opts.v0 = [reshape(gk,[],1);rand(1)];%starting vector %This affects the convergence a lot!! %we start from the tangent space!
        [v,~] = eigs(Afun,prob.n*prob.rank+1,1,'smallestreal',opts);
        normv = norm(v);

        vk = reshape(v(1:end-1),prob.n,prob.rank);
        tk = v(end);%the last element is scalar t

        if abs(tk) >= para.nu
            dk    = vk/tk;
        else
            dk = sign(-prob.M.inner(Xk, gk, vk))*vk;
        end
        
        eta = lineserach(Xk,dk,prob,para,@cost);
        Xk = prob.M.retr(Xk,dk,eta);%retraction 
        normgk = norm(prob.M.egrad2rgrad(Xk,egrad(Xk)),'fro');%euclidean gradient to Riemannian gradient
        obj = [obj;cost(C,Xk)];
        Grad = [Grad;normgk];
        if normgk<para.epislon
            stop = "grad reach desired accuracy!";
            break;
        end
        if (mod(iter,10)== 0)
            fprintf("%d     %3.1f  %.1e  %+.3e   %3.4f \n",iter,obj(iter),(obj(iter)-Truecost)/abs(Truecost),normgk,eta);
        end
    end
    if iter == para.Maxiter 
        stop = "reach maximum iteration!";
    end
    Out.iter = iter;
    Out.obj  = obj;
    Out.grad = Grad;
    Out.stop = stop;
    Out.para = para;
    Out.X    = Xk;
    fprintf("Algorithm stop: " + stop + "\n");
    fprintf("iter: %d, ngrad: %+.3e \n",iter,Grad(end));
    save(saveroot+name{1}+"-result.mat",'Out');
end




function y = F(x,Xk,prob,para) %Eigenvector routine 
    egrad = @(X) 2*prob.C*X;
    ehess = @(X, U) 2*U;
    gk = prob.M.egrad2rgrad(Xk,egrad(Xk)); 
    x1 = reshape(x(1:end-1),prob.n,prob.rank);
    x2 = x(end);
    y = [reshape(prob.M.ehess2rhess(Xk,egrad(Xk),ehess(Xk,x1),x1)+para.delta*gk,[],1);
        gk(:).'*x1(:) - para.delta*x2];
end

function y = cost(C,X)
    R = X*X';
    y = C(:).'*R(:);
end

function y = lineserach(X,dk,prob,para,Cost)
    etak = para.eta;
    C    = prob.C;
    Costant = para.gamma/6*norm(dk,'fro')^3*etak^3;
    iter = 1;
    while iter <= 10
        Dk= Cost(C,X) - Cost(C, prob.M.retr(X,dk,etak));
        if Dk >= Costant*etak^3
            break;
        else
            etak = para.beta*etak;
        end
        iter = iter + 1;
    end
    y = etak;
end


