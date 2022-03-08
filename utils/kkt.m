function [out] = kkt(qp,x,s,z,y)
    % KKT conditions for the QP 
    [Q,q,A,b,G,h] = unpack_qp(qp);
    out = [A'*y + G'*z + Q*x + q; % stationarity 
           s .* z;                % complimentarity
           G*x + s - h;           % primal feasibility
           A*x - b;               % primal feasability
           s - abs(s);            % dual feasability (s > 0) 
           z - abs(z)];           % dual feasibility (z > 0)
end