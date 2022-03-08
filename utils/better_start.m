function [x,s,z,y] = better_start(qp)
    idx = qp.idx;
    Ni = idx.nx + idx.nz + idx.ny;
    idx_y = idx.y - idx.ns;
    A = sparse(Ni, Ni);
    A(idx.x,idx.x) = qp.Q;
    A(idx.x,idx.s) = qp.G';
    A(idx.x,idx_y) = qp.A';
    A(idx.s,idx.x) = qp.G;
    A(idx.s,idx.s) = -eye(idx.ns);
    A(idx_y,idx.x) = qp.A;
    init = A\[-qp.q;qp.h;qp.b];
    
    x = init(idx.x);
    z = init(idx.s);
    y = init(idx_y);
    
    alpha_p = -min(-z);
    if alpha_p < 0
        s = -z;
    else
        s = -z + (1 + alpha_p);
    end
    
    alpha_d = -min(z);
    if alpha_d >= 0
        z = z + (1 + alpha_d);
    end
end
