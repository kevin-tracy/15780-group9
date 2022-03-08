
function [Q,q,A,b,G,h] = unpack_qp(qp)
    % unpacks the qp struct
    Q = qp.Q; q = qp.q; A = qp.A; b = qp.b; G = qp.G; h = qp.h;
end
