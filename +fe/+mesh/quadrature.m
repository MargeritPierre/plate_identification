function MESH = quadrature(MESH,scheme)
% Library of quadrature rules for function integration on meshes

if nargin<2 ; scheme = 'quad3' ; end
MESH.quadrature = scheme ;
    
switch MESH.quadrature
    case 'quad1' % second-order scheme
        MESH.xiQ = [0 0] ;
        MESH.wQ = 4 ;   
    case 'quad2' % second-order scheme
        MESH.xiQ = [-1 -1 ; 1 -1 ; 1 1 ; -1 1]*1/sqrt(3) ;
        MESH.wQ = [1 1 1 1] ;   
    case 'quad3' % third-order scheme
        MESH.xiQ = [-1 -1 ; 1 -1 ; 1 1 ; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0]*0.774596669241483 ;
        MESH.wQ = [25 ; 25 ; 25 ; 25 ; 40 ; 40 ; 40 ; 40 ; 64]/81 ;
end