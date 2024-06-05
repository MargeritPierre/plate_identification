function U = mindlin(PLATE,w)
% MINDLIN computes a series of harmonic plate responses using a FEM
% implementation of the mindlin isotropic thick-plate model
% mesh is the FEM mesh with a list of nodes, elements & a quadrature rule
% includePlaneMotion = true ;
% nDOFsByNode = 3 + 2*includePlaneMotion ;

% Generate the mesh
    %MESH = fe.mesh.gridmesh(PLATE.L,PLATE.dx) ; % quad mesh
    MESH = PLATE.mesh ;
    MESH = fe.mesh.quad8(MESH) ; % serendipity element
    MESH = fe.mesh.quadrature(MESH,'quad2') ; % generate a quadrature rule
    
% Problem matrices
    % Degrees of freedom: u = [u1;u2;u3;phi1;phi2]
    % Build the interpolation matrices
        MESH = fe.mesh.interpolation(MESH) ; % generate interpolation matrices
        [N,Dx,Wq,O] = deal(MESH.N,MESH.Dx,MESH.Wq,MESH.O) ;
        NWN = N'*Wq*N ;
    % In-Plane motion
        % strains Em = [E11;E22;2E12] = Bm*u
            Bm = [  ...
                    Dx{1} O O O O ; .... E11 = u1,1
                    O Dx{2} O O O ; ... E22 = u2,2
                    Dx{2} Dx{1} O O O ; ... 2E12 = u1,2+u2,1
                 ] ;
        % stiffness
            Km = Bm'*kron(sparse(PLATE.A),Wq)*Bm ;
        % inertia
            Mm = PLATE.m*kron(sparse(diag([1 1 0 0 0])),NWN) ;
    % Bending motion
        % strains Eb = [K11;K22;2K12] = Bb*u
            Bb = [  ...
                    O O O Dx{1} O ; .... K11 = phi1,1
                    O O O O Dx{2} ; ... K22 = phi2,2
                    O O O Dx{2} Dx{1} ; ... KE12 = phi1,2+phi2,1
                 ] ;
        % stiffness
            Kb = Bb'*kron(sparse(PLATE.D),Wq)*Bb ;
        % inertia
            Mb = PLATE.J*kron(sparse(diag([0 0 0 1 1])),NWN) ;
    % Out-of-plane shear motion
        % strains Es = [d1;d2] = Bs*u
            Bs = [  ...
                    O O Dx{1} N O ; .... d1 = phi1+u3,1
                    O O Dx{2} O N ; ... d2 = phi2+u3,2
                 ] ;
        % stiffness F
            Ks = (Bs'*kron(sparse(PLATE.F),Wq)*Bs) ;
        % inertia
            Ms = PLATE.m*kron(sparse(diag([0 0 1 0 0])),NWN) ;
    % Total system
        K = Km+Kb+Ks ;
        M = Mm+Mb+Ms ; M = .5*(M+M') ;
% Boundary conditions on plate edges
    tol = sqrt(eps) ;
    boundaryNodes = any([MESH.X<tol MESH.X>PLATE.L-tol],2) ;
    switch PLATE.bc
        case 'clamped'
            fix = [1 1 1 1 1] ; % fix all dofs
        case 'supported'
            fix = [1 1 1 0 0] ; % free rotations
        case 'free'
            fix = [0 0 0 0 0] ; % all dofs are free
    end
    keepDOF = ~(boundaryNodes & fix) ;
    K = K(keepDOF,keepDOF) ;
    M = M(keepDOF,keepDOF) ;
% Applied load
    [~,iF] = min(sum((MESH.X-PLATE.xf).^2,2)) ;
    f = zeros(size(MESH.X,1),5) ;
    f(iF,3) = 1 ;
    f = f(keepDOF) ;
% Solve for all frequencies
    U = zeros(size(MESH.X,1)*5,numel(w)) ;
    for ww = 1:numel(w)
        U(keepDOF,ww) = (K-w(ww)^2*M)\f ;
    end
    U = reshape(U,size(MESH.X,1),5,numel(w)) ;
% Extract only the transverse displacement
    U = reshape(U(1:prod(MESH.nX+1),3,:),[],numel(w)) ;
end

