
function [qp] = pack_qp(Q,q,A,b,G,h)
    % packs the qp struct
    qp.Q = Q; qp.q = q; qp.A = A; qp.b = b; qp.G = G; qp.h = h;
end

