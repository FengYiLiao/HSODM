function Out = HSODM(prob,para)
    %HSODM stops when Riemannian gradient reaches desired accuracy
    rng("default");
    
    obj  = [];
    Grad = [];
    if isfield(para,"X0")
        X0   = para.X0;
    else
        X0   = prob.M.rand(); %initial point
    end
    Xk   = X0;
    

    if ~isfield(para,"adp_delta")
        para.adp_delta = false;
    end

    if ~isfield(para,"delta_min")
        para.delta_min = 2; %default value
    end
    
    if ~isfield(para,"L")
        para.L = 2;%default value
    end

    if ~isfield(para,"para.")
        para.linesearch = true;
    end

    stop = "";
    fprintf("iter  |   obj  |  rgrad  |   eta   |   delta \n");
    rng(1);
    for iter = 1:para.Maxiter
        %Euclidean gradient to Riemannian gradient %As the initial guess for the power method
        gk      = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 

        %Define a routine for power method
        Afun    = @(x) prob.routine(x,Xk,prob,para);     

        %starting vector %This affects the convergence a lot!! %we start from the tangent space!
        opts.v0 = [reshape(gk,[],1);rand(1)];            
        [v,~]   = eigs(Afun,prob.n*prob.rank+1,1,'smallestreal',opts);

        %compensate computational error
        v       = real(v); 
        vk      = reshape(v(1:end-1),prob.n,prob.rank);%vector on the tangent space
        tk      = v(end);                              %the last element is scalar t

        %some critetion 
        if abs(tk) >= para.nu
            dk  = vk/tk;
        else
            dk  = sign(-prob.M.inner(Xk, gk, vk))*vk;
        end
        
        %line search for step size
        if para.linesearch
            eta    = lineserach(Xk,dk,prob,para);
        else
            eta    = 1;
        end
        
        
        %retraction 
        Xk     = prob.M.retr(Xk,dk,eta);
        normgk = norm(prob.M.egrad2rgrad(Xk,prob.egrad(Xk)),'fro');%euclidean gradient to Riemannian gradient
        obj    = [obj;prob.cost(Xk)];
        Grad   = [Grad;normgk];
        if normgk<para.epislon
            stop = "rgrad reaches desired accuracy!";
            break;
        end
        if (mod(iter,para.step)== 0)
            fprintf("%d     %3.1f  %+.3e   %3.4f   %3.4f\n",iter,obj(iter),normgk,eta,para.delta);
        end

        %adaptively change delta
        if para.adp_delta
%             if normgk <=10^-2
%                 para.delta = 2;
%             else
%                 para.delta = max(para.delta_min,normgk^(1/2));
%             end
             para.delta =normgk^(1/1.5);
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
    Costant = para.gamma/6*norm(dk,'fro')^3*etak^3;
    iter = 1;
    while iter <= 30
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
