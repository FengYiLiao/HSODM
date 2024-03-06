function Out = HSODM_2(prob,option)
    %HSODM stops when Riemannian gradient reaches desired accuracy
    %This script considers the homogenous second-order descent framework
    %This code is not working!!!!!!!!!!!!
    rng("default");
    
    obj  = [];
    Grad = [];

    para = initialization(prob,option);
    Xk   = para.X0;
   
    stop = "";
    fprintf("iter  |   obj  |  rgrad  |   eta   |   delta \n");
    rng(1);
    pk = 0.5;%%%%%%%%%%%%
    for iter = 1:para.Maxiter

        %find a proper I_h
        if iter > 1 
            if pk > para.eta2 
                Ih_min = max(para.hmin,para.r4 * sqrt(hk));
                Ih_max = sqrt(hk);
            elseif (para.eta1 <= pk) && (pk<=para.eta2) 
                Ih_min = sqrt(hk)/para.r1;
                Ih_max = para.r2*sqrt(hk);                
            else
                Ih_min = para.r2*sqrt(hk);
                Ih_max = para.r3*sqrt(hk);
            end
        end


        %Euclidean gradient to Riemannian gradient %As the initial guess for the power method
        gk         = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 


        if iter > 1
        
        
            %Hessian    = prob.M.ehess2rhess(Xk,prob.egrad(Xk),prob.ehess(Xk,x1),x1);
            Hessian    = @(x) prob.Hessian_routine(x,Xk,prob);  
            optsh.v0   = reshape(gk,[],1);
            [v_l,e_l]      = eigs(Hessian,prob.n*prob.rank,1,'largestreal',optsh);
            [v_s,e_s]      = eigs(Hessian,prob.n*prob.rank,1,'smallestreal',optsh);
    
            %compensate computational error
            v_l          = real(v_l);
            e_l          = real(e_l);
            v_s          = real(v_s);
            e_s          = real(e_s);
    
            delta_low  = min(e_s,-sqrt(Ih_max));
            delta_high = max(((1+e_s)^2*(1+e_l)+abs(e_s))/Ih_min,norm(prob.M.egrad2rgrad(Xk,prob.egrad(Xk)),'fro')^2);
    
            delta_low  = -2;
            delta_high = -1;


            %Define a routine for power method
            Afun_low    = @(x) prob.routine_2(x,Xk,prob,delta_high); 
            Afun_high   = @(x) prob.routine_2(x,Xk,prob,delta_low); 
    
    
            opts.v0  = [reshape(gk,[],1);rand(1)];  
    
            %low
            [v1,e1]  = eigs(Afun_low,prob.n*prob.rank+1,1,'smallestreal',opts);
    
            %compensate computational error
            v1       = real(v1); 
            e1       = real(e1);
            vk1      = reshape(v1(1:end-1),prob.n,prob.rank);%vector on the tangent space
            tk1      = v1(end);                              %the last element is scalar t    
            dk1      = vk1/tk1;
            
            h_low    = h(-e1,dk1); %the value we need in the bisection 
    
    
            %high
            [v2,e2]  = eigs(Afun_high,prob.n*prob.rank+1,1,'smallestreal',opts);
    
            %compensate computational error
            v2       = real(v2); 
            e2       = real(e2);  
            vk2      = reshape(v2(1:end-1),prob.n,prob.rank);%vector on the tangent space
            tk2      = v2(end);                              %the last element is scalar t    
            dk2      = vk2/tk2;
       
            h_high   = h(-e2,dk2); %the value we need in the bisection
    
            h1       = h_low;
            h2       = h_high;
            delta_1  = delta_high;
            delta_2  = delta_low;
    
            while true
                
        
                %Define a routine for power method
                Afun    = @(x) prob.routine(x,Xk,prob,para);     
        
                %starting vector %This affects the convergence a lot!! %we start from the tangent space!
                opts.v0 = [reshape(gk,[],1);rand(1)];            
                [v,e]   = eigs(Afun,prob.n*prob.rank+1,1,'smallestreal',opts);
        
                %compensate computational error
                v       = real(v); 
                    
                vk      = reshape(v(1:end-1),prob.n,prob.rank);%vector on the tangent space
                tk      = v(end);                              %the last element is scalar t
                
                dk      = vk/tk;
                hk      = h(e,dk);
    
    
                if (Ih_min<= hk) && (hk <= Ih_max)
                    break;
                else
                    %bisection
                    midpoint_h        = (h1 + h2)/2;
                    midpoint_delta    = (delta_1+delta_2)/2;
                    if Ih_max<= midpoint_h
                       h2       = midpoint_h;
                       delta_1  = midpoint_delta;
                    else
                       h1 = midpoint_h;
                       delta_2  = midpoint_delta;
                    end
                    para.delta  = midpoint_delta;
                end
            end
    
            %some critetion 
    %         if abs(tk) >= para.nu
    %             dk  = vk/tk;
    %         else
    %             dk  = sign(-prob.M.inner(Xk, gk, vk))*vk;
    %         end
            
            %Compute 
            Xk_next   = prob.M.retr(Xk,dk,1);
            mk        = prob.cost(Xk) + trace(gk.'*dk) + 0.5*trace(gk.'*prob.M.ehess2rhess(Xk,prob.egrad(Xk)))+ sqrt(hk)/3*dk^3;
            pk        = (prob.cost(Xk_next) - prob.cost(Xk))/(mk-prob.cost(Xk));
    
            if pk>= para.eta1 
                Xk = Xk_next;
            end
        else
            %Euclidean gradient to Riemannian gradient %As the initial guess for the power method
            gk      = prob.M.egrad2rgrad(Xk,prob.egrad(Xk)); 
    
            %Define a routine for power method
            Afun    = @(x) prob.routine(x,Xk,prob,para);     
    
            %starting vector %This affects the convergence a lot!! %we start from the tangent space!
            opts.v0 = [reshape(gk,[],1);rand(1)];            
            [v,e]   = eigs(Afun,prob.n*prob.rank+1,1,'smallestreal',opts);
    
            %compensate computational error
            v       = real(v); 
            vk      = reshape(v(1:end-1),prob.n,prob.rank);%vector on the tangent space
            tk      = v(end);                              %the last element is scalar t
    
            dk  = vk/tk;
%             %some critetion 
%             if abs(tk) >= para.nu
%                 dk  = vk/tk;
%             else
%                 dk  = sign(-prob.M.inner(Xk, gk, vk))*vk;
%             end
            
            e = real(e);
            
            %line search for step size
%             eta    = lineserach(Xk,dk,prob,para);

            hk      = h(e,dk);
        end


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
            if normgk <=10^-2
                para.delta = 2;
            else
                para.delta = max(para.delta_min,normgk^(1/2));
            end
        end
    end
    if iter == para.Maxiter 
        stop = "reach maximum iteration!";
    end
    Out.iter = iter;
    Out.obj  = obj;
    Out.grad = Grad;
    Out.stop = stop;
    Out.para = option;
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


function y = h(theta,d)
    temp  = norm(d,'fro');
    temp1 = theta/norm(d,'fro');
    y     = (theta/norm(d,'fro'))^2;
end

function para = initialization(prob,option)
    para = option;
    
    if ~isfield(para,"X0")
%         X0   = para.X0;
%     else
        para.X0   = prob.M.rand(); %initial point
    end

    if ~isfield(para,"adp_delta")
        para.adp_delta  = false;
        para.adp_delta = false;
    end

    if ~isfield(para,"delta_min")
        para.delta_min   = 2;
        option.delta_min = 2; %default value
    end

    if isfield(para,"r1")
       para.r1 = option.r1;
    else
       para.r1 = 1.5;
    end

    if isfield(para,"r2")
       para.r2 = option.r2;
    else
       para.r2 = 2;
    end
        
    if isfield(para,"r3")
       para.r3 = option.r3;
    else
       para.r3 = 2.5;
    end

    if isfield(para,"r4")
       para.r4 = option.r4;
    else
       para.r4 = 0.5;
    end

    if isfield(para,"eta1")
       para.eta1 = option.eta1;
    else
       para.eta1 = 0.1;
    end

    if isfield(para,"eta2")
       para.eta2 = option.eta2;
    else
       para.eta2 = 0.9;
    end

    if isfield(para,"hmin")
       para.hmin = option.hmin;
    else
       para.hmin = 0.1;
    end    

end
