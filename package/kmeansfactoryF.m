function M = kmeansfactoryF(n,r)
   
    M.St = stiefelfactory(n, r, 1);
    

%     M.egrad2rgrad = @egrad2rgrad;
% 
%     function rgrad = egrad2rgrad(X,egrad)
%         
%     end
    M.inner = @(x, d1, d2) d1(:)'*d2(:);
    M.trnsp = @(X) X;
    M.proj  = @projection;

%     function up = projection(X,W) 
%         l     = ones(n,1);
%         x     = 1/n*(W*X'*l);
%         C     = x*l'+l*x';                %avoid duplicated computations
%         Sigma = 0.25*(W'*X+X'*W-2*X'*C*X);
%         up    = W - 2*X*Sigma - C*X;
%     end

    function up = projection(X,W)
        n     = height(W);
        q     = width(W);
        C     = X'*W; %avoid duplicated computations
        I     = eye(q);
        e     = ones(n,1);
        alpha = X'*e;
        alpha = alpha/norm(alpha);
        up = 0.5*X*(C-C')+(I-X*X')*W*(I-alpha*alpha');
    end


    M.egrad2rgrad = M.proj; 

    M.retr  = @retraction;


%     function Y = retraction(X,U,t)
%         A  = t*X'*U;
%         Ap = X*A*X';
%         B  = t*(U*X'-X*U')-2*Ap;
%         Y  = expm(B)*expm(Ap)*X;
%     end


    M.projFv  = @projectionFv;
    function Y = projectionFv(X)
        n     = height(X);
        q     = width(X);
        e     = ones(n,1);
        I     = eye(q);
        alpha = X'*e;
        alpha = alpha/norm(alpha);
        Y     = e*alpha'/(norm(e))+ X *(I - alpha*alpha');
    end


    function Y = retraction(X,U,t)
       temp = M.St.retr(X,t*U);
       Y    = projectionFv(temp);
    end


    M.norm    = @(x, d) norm(d(:));
    M.lincomb = @matrixlincomb;
    M.transp  = @(x1, x2, d) M.proj(x2, d);
    M.Hess    = @Hess;
    
    function Y = Hess(X,V)
        n           = height(X);
        q           = width(X);
        e           = ones(n,1);
        I           = eye(q);
        pj          = X;
        idx         = X <= 0;
        pj(idx)     = 0 ;
        gradpj      = zeros(n,q);
        gradpj(idx) = 1;
        temp1       = X'*(gradpj.*V);
        temp2       = V'*pj; 
        temp3       = e*e'*X;
        alpha       = X'*e;
        norm_alpha  = norm(alpha);
        alpha       = alpha/norm(alpha);
        A           = X*(temp1+temp2-temp1-temp2);
        B           = 2*(X*V'+V*X')*pj*(I - alpha*alpha');
        C           = 2*(I-X*X')*pj*((alpha*alpha'/sqrt(norm_alpha)*2)*temp3(:)'*V(:)-X*e*e'*V/norm_alpha/norm_alpha);
        Y           = A + B + C;
        Y           = projection(X,Y);
    end
    
    
    
end