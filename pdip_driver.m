clear

% run utils.m

% nx = 100 ;
% ns = 24;
% nz = ns ;
% ny = 8;
nx = 30 ;
ns = 14;
nz = ns ;
ny = 8;
% 
[qp] = gen_qp(nx,ns,ny);
% qp.Q = qp.Q + 1e2*eye(length(qp.q));
% qp = create_mpc(10);
% 
verbose = true;
atol = 1e-10;
easy_start = true;
max_iters = 10;
reg = 1e-6;
% 
% [x,s,z,y,cond_hist,iter,success] = solve_qp_lu(qp, verbose, atol, easy_start, max_iters);
[x,s,z,y,cond_hist,iter,success] = solve_qp_ldl(qp, verbose, atol, easy_start, max_iters,reg);

N = 100 
lu_cond = cell(N,1);
ldl_cond = cell(N,1)
for i = 1:N
    [qp] = gen_qp(nx,ns,ny);
    [x,s,z,y,lu_cond{i},iter,success] = solve_qp_lu(qp, verbose, atol, easy_start, max_iters);
    [x,s,z,y,ldl_cond{i},iter,success] = solve_qp_ldl(qp, verbose, atol, easy_start, max_iters,reg);

end
     


figure
hold on
% plot(cond_hist)
for i = 1:N
    plot(lu_cond{i},'b')
    plot(ldl_cond{i},'r')
    
end

p1 = plot(NaN,NaN,'b');
p2 = plot(NaN,NaN,'r')
legend([p1,p2],{'Original KKT', 'Symmetrized KKT'})
set(gca,'YScale','log')
xlabel('Iteration')
ylabel('Condition number')
hold off 
