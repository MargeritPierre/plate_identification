function MESH = gridmesh(L,dx)
% Construct a uniform grid mesh of quad elements
% L: grid length in each dimensions
% dx: approximate element size

% BUILD THE QUAD MESH
    L = L + 0*dx ; dx = dx + 0*L ;
    nX = floor(L./dx) ; % number of elements in each direction
    dx = L./nX ;
    X = arrayfun(@linspace,L*0,L,nX+1,'uni',false) ;
    [X{:}] = ndgrid(X{:}) ;
    X = cat(3,X{:}) ;
    X = reshape(X,[],2) ;
    iX = reshape(1:prod(nX+1),nX+1) ;
    Elems = cat(3,iX(1:end-1,1:end-1),iX(1:end-1,2:end),iX(2:end,2:end),iX(2:end,1:end-1)) ;
    Elems = reshape(Elems,[],4) ;
    
% Mesh structure
    MESH = [] ;
    MESH.L = L ;
    MESH.dx = dx ; 
    MESH.nX = nX ;
    MESH.X = X ; 
    MESH.Elems = Elems ;
    
% Element interpolation 
    % shape functions, xi in [-1;1]x[-1;1]
        MESH.ShapeFcn = @(xi) [ ...
                           1/4*(1-xi(:,1)).*(1-xi(:,2)) ...
                           1/4*(1+xi(:,1)).*(1-xi(:,2)) ...
                           1/4*(1+xi(:,1)).*(1+xi(:,2)) ...
                           1/4*(1-xi(:,1)).*(1+xi(:,2)) ...
                         ] ;
    % shape function derivatives
        MESH.dShapeFcn_dxi = @(xi) {[ ...
                                   -1/4.*(1-xi(:,2)) ...
                                   1/4.*(1-xi(:,2)) ...
                                   1/4.*(1+xi(:,2)) ...
                                   -1/4.*(1+xi(:,2)) ...
                              ];[ ...
                                   -1/4*(1-xi(:,1)) ...
                                   -1/4*(1+xi(:,1)) ...
                                   1/4*(1+xi(:,1)) ...
                                   1/4*(1-xi(:,1)) ...
                              ]} ;
    
end

