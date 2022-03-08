
function logging(qp,s,x,z,y,alpha,iter)
    [Q,q,A,b,G,h] = unpack_qp(qp);
    % rhs = -[A'*y + G'*z + Q*x + q; 
    J = 0.5*x'*Q*x + q'*x;
    gap = dot(s,z)/length(s);
    eq_res = norm(A*x -b);
    ineq_res = norm(G*x + s - h);
    disp(sprintf('%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n',...
                 [iter, J, gap, eq_res, ineq_res, alpha]))
end