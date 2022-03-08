function qp = create_MPC(N)
    nx = 6
    nu = 3
    xg = [-1;-2;-3;0;0;0]
    nz=(N*nx) + (N-1)*nu
    Ac = [zeros(3,3) eye(3);zeros(3,6)]
    Bc = [zeros(3,3);eye(3)]

    dt = 0.2
    H = expm(dt*[Ac Bc; zeros(3,9)])
    Ad = sparse(H(1:nx,1:nx))
    Bd = sparse(H(1:nx,nx+1:end));

%     xi = [(nu+nx)*(i-1) .+ (1:nx)  for i = 1:N]
% 
%     ui = [(nu+nx)*(i-1) .+ (nx .+ (1:nu))  for i = 1:N-1]
%     ci = [(nx)*(i-1) .+ (1:nx)  for i = 1:N]

    % cost function
    Qc = 10*eye(nx);
    Rc = eye(nu);
    Q = sparse(nz,nz);
    q = zeros(nz,1);
    for i = 1:N
        if i<N
            Q(xi(i,nu,nx),xi(i,nu,nx)) = Qc;
            q(xi(i,nu,nx)) = -Qc*xg;
            Q(ui(i,nu,nx),ui(i,nu,nx)) = Rc;
        else
            Q(xi(i,nu,nx),xi(i,nu,nx)) = Qc;
            q(xi(i,nu,nx)) = -Qc*xg;
        end
    end


    A = sparse(N*nx,nz);
    for i = 1:N-1
        % xkp1 - axk - buk
        A(ci(i,nu,nx),xi(i+1,nu,nx)) = eye(nx);
        A(ci(i,nu,nx),xi(i,nu,nx))= -Ad;
        A(ci(i,nu,nx),ui(i,nu,nx)) = -Bd;
    end
    A(ci(N,nu,nx),xi(1,nu,nx)) = eye(nx);
    b = zeros((N*nx),1);
    b(ci(N,nu,nx)) = [1;4;-8;.1;.2;-.3];

    G_up = sparse(nz,nz);
%     gi = [(nu)*(i-1) .+ (1:nu)  for i = 1:N-1]
    for i = 1:N-1
        idx = (nu)*(i-1) + (1:nu);
        G_up(idx,ui(i,nu,nx)) = eye(nu);
    end
    G = [G_up;-G_up];
    h = 2*ones(2*nz,1);

    qp = pack_qp(Q,q,A,b,G,h);
end


function idx = xi(i,nu,nx)
idx = (nu+nx)*(i-1) + (1:nx);
end
function idx = ui(i,nu,nx)
idx = (nu+nx)*(i-1) + (nx + (1:nu));
end
function idx = ci(i,nu,nx)
idx = (nx)*(i-1) + (1:nx) ;
end
