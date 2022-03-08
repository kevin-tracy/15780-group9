function [x,s,z,y,cond_hist,iter,success] = solve_qp_ldl(qp, verbose, atol, easy_start, max_iters)

    qp.idx = gen_idx(qp);
    if easy_start
        [x,s,z,y] = naive_start(qp);
    else
        [x,s,z,y] = better_start(qp);
    end
    
    if verbose
        disp("iter     objv        gap       |Ax-b|    |Gx+s-h|    step")
        disp("---------------------------------------------------------")
    end
    
    cond_hist = zeros(max_iters,1);
    success = false;
    for iter = 1:max_iters
    
        if norm(kkt(qp,x,s,z,y))<atol
            cond_hist = cond_hist(1:iter-1);
            success = true;
            break
        end
    
        % get jacobian 
        K = kkt_jacobian(qp,s,z);
        cond_hist(iter) = condest(K);
        K = decomposition(K,'lu');
    
        % affine step 
        rhs_a = rhs_affine(qp,x,s,z,y);
        beta_a = K\rhs_a;
        ds_a = beta_a(qp.idx.s);
        dz_a = beta_a(qp.idx.z);
    
        % centering params 
        [sigma, mu] = centering_params(qp,ds_a, dz_a, s,z);
    
        % centering setp 
        rhs_c = rhs_centering(qp,ds_a, dz_a, sigma, mu);
        beta_c = K\rhs_c;
    
        % combine steps
        dbeta = beta_a + beta_c;
        ds = dbeta(qp.idx.s);
        dz = dbeta(qp.idx.z);
    
        % update primal and dual variables
        alpha = min(1,0.99*min(linesearch(s,ds), linesearch(z,dz)));
        x = x + alpha*dbeta(qp.idx.x);
        s = s + alpha*dbeta(qp.idx.s);
        z = z + alpha*dbeta(qp.idx.z);
        y = y + alpha*dbeta(qp.idx.y);
    
        if verbose
            logging(qp,s,x,z,y,alpha,iter)
        end
        
    
    end
end

function [K] = kkt_jacobian(qp,s,z)
idx = qp.idx;
[Q,q,A,b,G,h] = unpack_qp(qp);
K = sparse(idx.nt,idx.nt);
K(idx.x,idx.x) = Q;
K(idx.x,idx.z) = G';
K(idx.x,idx.y) = A';
K(idx.s,idx.s) = diag(z);
K(idx.s,idx.z) = diag(s);
K(idx.z,idx.x) = G;
K(idx.z,idx.s) = eye(idx.nz);
K(idx.y,idx.x) = A;
end

function rhs = rhs_affine(qp,x,s,z,y)
[Q,q,A,b,G,h] = unpack_qp(qp);
rhs = -[A'*y + G'*z + Q*x + q; 
       s .* z;               
       G*x + s - h;           
       A*x - b];
end

function rhs = rhs_centering(qp,ds_a, dz_a, sigma, mu)
rhs = zeros(qp.idx.nt,1);
rhs(qp.idx.s) = sigma*mu - (ds_a .* dz_a);
end



