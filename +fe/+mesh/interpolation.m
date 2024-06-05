function MESH = interpolation(MESH)
% RETRUNS MESH INTERPOLATION MATRICES ASSOCIATED TO THE MESH QUADRATURE
% RULE
% f(x_qp) = N*f(x_n)
% df(x_qp)_dxi = Dx{i}*f(x_n)
% Wq contains the integration weights
% O is the zero matrix with the same size than N or Dx{i}

% Mesh info
    [nNodes,nCoord] = size(MESH.X) ; 
    [nElems,nNodesByElem] = size(MESH.Elems) ;
        
% Mesh interpolation
    % full quadrature
        nQPe = numel(MESH.wQ) ;
        nQP = nQPe*nElems ;
        XIq = repmat(MESH.xiQ,[nElems 1]) ;
    % Sparse matrices
        ii = repelem((1:nQP)',1,nNodesByElem) ;
        jj = repelem(MESH.Elems,nQPe,1) ;
    % f(XIq) = N*f(X)
        N = sparse(ii,jj,MESH.ShapeFcn(XIq),nQP,nNodes) ;
        %patch('vertices',N*X,'faces',(1:nQP)','marker','+') % plot quadrature points
    % df(XIq)_dxi = Dxi*f(X)
        Dxi = cellfun(@(dd)sparse(ii,jj,dd,nQP,nNodes),MESH.dShapeFcn_dxi(XIq),'uni',false) ;
    % df(XIq)_dX = Dx*f(X)
        dX_dxi = [Dxi{1}*MESH.X Dxi{2}*MESH.X] ; % [dX1_dxi1 dX2_dxi1 dX1_dxi2 dX2_dxi2]
        detJ = dX_dxi(:,1).*dX_dxi(:,4)-dX_dxi(:,2).*dX_dxi(:,3) ;
        iJ = (1./detJ).*[dX_dxi(:,4) -dX_dxi(:,3) -dX_dxi(:,2) dX_dxi(:,1)] ; % [dxi1_dX1 dxi2_dX1 dxi1_dX2 dxi2_dX2]
        Dx = cell(nCoord,1) ;
        Dx{1} = diag(sparse(iJ(:,1)))*Dxi{1} + diag(sparse(iJ(:,2)))*Dxi{2} ;
        Dx{2} = diag(sparse(iJ(:,3)))*Dxi{1} + diag(sparse(iJ(:,4)))*Dxi{2} ;
        %norm(Dx{1}*X-[1 0]),norm(Dx{2}*X-[0 1]) % check derivatives
    % Quadrature weights
        Wq = diag(sparse(detJ.*repmat(MESH.wQ(:),[nElems 1]))) ;
        %abs(sum(diag(Wq))-prod(L)) % check integration
    % Zero matrix
        O = sparse(nQP,nNodes) ;
        
        
% PUT IN THE MESH
    MESH.N = N ;
    MESH.Dx = Dx ;
    MESH.Wq = Wq ;
    MESH.O = O ;
        
end

