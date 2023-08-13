function Out = HSODM(prob,para)
    obj = [];
    Grad = [];
    X0 = prob.M.rand(); %initial point
    Xk = X0;
    stop ="";
    fprintf("iter  |   obj  |  rgrad  |   eta   \n");
    for iter = 1:para.Maxiter
        
        gk = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); %euclidean gradient to Riemannian gradient
        Afun = @(x) prob.routine(x,Xk,prob,para);%Define a routine
        opts.v0 = [reshape(gk,[],1);rand(1)];%starting vector %This affects the convergence a lot!! %we start from the tangent space!
        [v,~] = eigs(Afun,prob.n*prob.rank+1,1,'smallestreal',opts);

        vk = reshape(v(1:end-1),prob.n,prob.rank);%vector on the tangent space
        tk = v(end);%the last element is scalar t

        if abs(tk) >= para.nu
            dk    = vk/tk;
        else
            dk = sign(-prob.M.inner(Xk, gk, vk))*vk;
        end
        
        eta = lineserach(Xk,dk,prob,para);
        Xk = prob.M.retr(Xk,dk,eta);%retraction 
        normgk = norm(prob.M.egrad2rgrad(Xk,prob.egrad(Xk)),'fro');%euclidean gradient to Riemannian gradient
        obj = [obj;prob.cost(Xk)];
        Grad = [Grad;normgk];
        if normgk<para.epislon
            stop = "rgrad reaches desired accuracy!";
            break;
        end
        if (mod(iter,10)== 0)
            fprintf("%d     %3.1f  %+.3e   %3.4f \n",iter,obj(iter),normgk,eta);
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

end

function y = lineserach(X,dk,prob,para)
    etak = para.eta;
    %C    = prob.C;
    Costant = para.gamma/6*norm(dk,'fro')^3*etak^3;
    iter = 1;
    while iter <= 10
        Dk= prob.cost(X) - prob.cost(prob.M.retr(X,dk,etak));
        if Dk >= Costant*etak^3
            break;
        else
            etak = para.beta*etak;
        end
        iter = iter + 1;
    end
    y = etak;
end