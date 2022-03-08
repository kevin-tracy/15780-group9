function [sigma, mu] = centering_params(qp,ds_a, dz_a, s,z)
    mu = dot(s,z)/qp.idx.ns;

    alpha = min(linesearch(s,ds_a), linesearch(z,dz_a));

    sigma = (dot(s + alpha*ds_a, z + alpha*dz_a)/dot(s,z))^3;
end
