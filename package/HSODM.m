function Out = HSODM(prob,para)
    %HSODM stops when Riemannian gradient reaches desired accuracy
    rng(2019);
    
    obj  = [];
    Grad = [];
    if isfield(prob,"X0")
        X0   = prob.X0;
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

    if ~isfield(para,"linesearch")
        para.linesearch = true;
    end

    if ~isfield(prob,"nummanifold")
        prob.nummanifold = 1;
    end

    if ~isfield(prob,'mani2vec')
        prob.mani2vec = @(gk,prob) reshape(gk,[],1);
    end

    stop = "";
    fprintf("iter  |   obj  |  rgrad  |   eta   |   delta \n");
    rng(1);
    prob.storedb = StoreDB(2);
    key = prob.storedb.getNewKey();

    for iter = 1:para.Maxiter
        %Euclidean gradient to Riemannian gradient %As the initial guess for the power method
        %gk      = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 

        %Create a store database and get a key for the current x

       
        
        gk = getGradient(prob, Xk,prob.storedb,key);
        %gk      = prob.M.egrad2rgrad(Xk,prob.egrad(Xk,storedb),storedb,key); 


        %Define a routine for power method
        Afun    = @(x) prob.routine(x,Xk,prob,para);     

        %starting vector %This affects the convergence a lot!! %we start from the tangent space!
        
 
        opts.v0 = prob.mani2vec(gk,prob);
        opts.v0 = [opts.v0;rand(1)];
        opts.tol = 1e-9;
        if iter > 1
            if normgk <= 1e-2
                para.eta   = 1;
            end
        end
        opts.maxit = 150;
        opts.p     = 50; %dimenion of the Kylov subspace 
        opts.fail  = 'keep';
%         if para.nummanifold  == 1
%             opts.v0 = [reshape(gk,[],1);rand(1)];
%         else
%             opts.v0 = prob.mani2vec(gk,prob);
%             %opts.v0 = zeros(prob.Tolvar+1);
%         end
        
        [v,~]   = eigs(Afun,prob.Tolvar+1,1,'smallestreal',opts);
        
        %compensate computational error
        %v       = real(v); 

        [vk,tk] = prob.vec2mani(v,prob);
        tk      = tk/para.ck;

%         vk      = reshape(v(1:end-1),prob.n,prob.rank);%vector on the tangent space
%         tk      = v(end);                              %the last element is scalar t

        %some critetion (Dead) %para.nu = 0
        if abs(tk) >= para.nu
            if prob.nummanifold == 1
                dk   = vk/tk;
            else
                dk   = vk;
                dk.U = vk.U/tk;
                dk.V = vk.V/tk;
            end
        else
            dk  = sign(-prob.M.inner(Xk, gk, vk))*vk;
        end
       
        
        %normdk = norm(dk,'fro');
        %normdk = prob.M.norm(Xk,dk);

        %line search for step size
        if para.linesearch
            eta = lineserach(Xk,dk,prob,para);
        else
            eta = para.eta;
        end
        
        
        %retraction 
        Xk     = prob.M.retr(Xk,dk,eta);
        %Xk     = prob.M.exp(Xk,dk,eta);

        %normgk = norm(prob.M.egrad2rgrad(Xk,prob.egrad(Xk)),'fro');%euclidean gradient to Riemannian gradient
        
%         %Create a store database and get a key for the current x
%         storedb = StoreDB(2);
%         key = storedb.getNewKey();

        newkey  = prob.storedb.getNewKey();
        [fk,gk] = getCostGrad(prob, Xk,prob.storedb,newkey);
        normgk  = prob.M.norm(Xk,gk);
        %norm(gk,'fro');


%         obj    = [obj;prob.cost(Xk)];
%         Grad   = [Grad;normgk];

        obj    = [obj;fk];
        Grad   = [Grad;normgk];

        if normgk<para.epislon
            stop = "rgrad reaches desired accuracy!";
            break;
        end
        if (mod(iter,para.step)== 0)
            fprintf("%d     %3.1f  %+.3e   %3.4f   %3.4f\n",iter,obj(iter),normgk,eta,para.delta);
        end

        %adaptively change delta
        if para.adp_delta && normgk <=1
%             if normgk <=10^-2
%                 para.delta = 2;ß
%             else
%                 para.delta = max(para.delta_min,normgk^(1/2));
%             end
             %para.delta =normgk^(1/1.5);
             %para.delta = sqrt(para.epislon);
             para.delta = 1e-3;
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
    Out.name = "HSODM";
    fprintf("Algorithm stop: " + stop + "\n");
    fprintf("iter: %d, ngrad: %+.3e \n",iter,Grad(end));
end

function y = lineserach(X,dk,prob,para)
    etak = para.eta;
    
    if (prob.nummanifold == 1)
        cartesian = false;
    else
        cartesian = true;
        dk_copy = dk;
    end
    if cartesian
        dk_copy.U = dk.U * etak;
        dk_copy.V = dk.V * etak;
        Costant = para.gamma/6*prob.M.norm(X,dk_copy)^3;
    else
        Costant = para.gamma/6*prob.M.norm(X,etak*dk)^3;
    end
    
    iter = 1;
    cost1 = getCost(prob, X);
    while iter <= 30
        %test = canGetCost(prob,Xk);
        
        cost2 = getCost(prob, prob.M.retr(X,dk,etak));
        %Dk= prob.cost(X) - prob.cost(prob.M.retr(X,dk,etak));
        Dk= cost1 - cost2;
        if Dk >= Costant
            break;
        else
            etak = para.beta*etak;
        end
        iter = iter + 1;
        if cartesian
            dk_copy.U = dk.U * etak;
            dk_copy.V = dk.V * etak;
            Costant = para.gamma/6*prob.M.norm(X,dk_copy)^3;
        else
            Costant = para.gamma/6*prob.M.norm(X,etak*dk)^3;
        end
       % Costant = para.gamma/6*prob.M.norm(X,etak*dk)^3;
    end
    y = etak;
end
