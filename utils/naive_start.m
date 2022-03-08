function [x,s,z,y] = naive_start(qp)
    x = zeros(qp.idx.nx,1);
    y = zeros(qp.idx.ny,1);
    s = ones(qp.idx.ns,1);
    z = ones(qp.idx.nz,1);
end
