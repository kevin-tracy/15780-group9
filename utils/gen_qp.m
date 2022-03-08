
function [qp] = gen_qp(nx,ns,ny)
    % generate a random QP by starting with optimal primal and dual variables,
    % and using the KKT conditions to complete the problem definition
    nz = ns;
    fill_density = 0.2;
    x = randn(nx,1);
    Q = sprandn(nx,nx,fill_density/10);
    Q = Q'*Q + 1e-1*eye(nx);
    G = sprandn(ns,nx,fill_density);
    s = abs(randn(ns,1));
    z = abs(randn(ns,1));
    % this encodes which constraints are active
    for i = 1:ns
        if rand < 0.5
            s(i) = 0 ;
        else
            z(i) = 0 ;
        end
    end
    h = G*x + s;
    A = sprandn(ny,nx,fill_density);
    if sprank(A)<size(A,1)
        A(1:ny,1:ny) = eye(ny);
    end
    b = A*x;
    y = randn(ny,1);
    q = -Q*x - G'*z - A'*y;
    [qp] = pack_qp(Q,q,A,b,G,h);
    % make sure we setup the problem right
    assert(norm(kkt(qp,x,s,z,y))<1e-10)
end
