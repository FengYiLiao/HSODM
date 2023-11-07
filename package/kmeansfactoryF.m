function M = kmeansfactory(n,r)
   
%     M.egrad2rgrad = @egrad2rgrad;
% 
%     function rgrad = egrad2rgrad(X,egrad)
%         
%     end
    M.inner = @(x, d1, d2) d1(:)'*d2(:);
    M.trnsp = @(X) X;
    M.proj  = @projection;

    function up = projection(X,W) 
        l     = ones(n,1);
        x     = 1/n*(W*X'*l);
        C     = x*l'+l*x';                %avoid duplicated computations
        Sigma = 0.25*(W'*X+X'*W-2*X'*C*X);
        up    = W - 2*X*Sigma - C*X;
    end

    M.egrad2rgrad = M.proj; 

    M.retr  = @retraction;


    function Y = retraction(X,U,t)
        A  = t*X'*U;
        Ap = X*A*X';
        B  = t*(U*X'-X*U')-2*Ap;
        Y  = expm(B)*expm(Ap)*X;
    end

    M.norm    = @(x, d) norm(d(:));
    M.lincomb = @matrixlincomb;
    M.transp = @(x1, x2, d) M.proj(x2, d);
end