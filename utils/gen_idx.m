
function idx = gen_idx(qp)
    nx = length(qp.q);
    ns = length(qp.h);
    ny = length(qp.b);
    % create indices for stacked [x;s;z;y] 
    nz = ns;
    idx.nx = nx; idx.ns = ns; idx.nz = nz; idx.ny = ny;
    idx.x = 1:nx;
    idx.s = (idx.x(end) + 1) : (idx.x(end) + idx.ns);
    idx.z = (idx.s(end) + 1) : (idx.s(end) + idx.nz);
    idx.y = (idx.z(end) + 1) : (idx.z(end) + idx.ny);
    idx.nt = idx.nx + idx.ns + idx.nz + idx.ny; % total # of variables 
end