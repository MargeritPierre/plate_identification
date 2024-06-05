function p = plotmesh(mesh,u,varargin)
% plot a mesh a a function if needed

% Initialize the mesh patch
    p = patch(...
            'vertices',mesh.X...
            ,'faces',mesh.Elems...
            ,'FaceColor','w'...
            ,'marker','.'...
            ,varargin{:}) ;

% Plot the function
    if nargin>=2
    % Function interpolation
        switch size(u,1)
            case 0 % empty field
            case size(mesh.X,1) % 1 value by node
                p.FaceColor = 'interp' ;
            case size(mesh.Elems,1) % 1 value by element
                p.FaceColor = 'flat' ;
            otherwise
                error('wrong shape for the mesh function u') ;
        end
    % Function components
        switch size(u,2) % number of function components ?
            case 0 % empty data
                u = mesh.X*0 ;
            case 1 % assume u3 was given
                u = padarray(u,[0 2],0,'pre') ; % add [u1 u2]
            case 2 % assume [u1 u2] was given
                u = padarray(u,[0 1],0,'post') ;  % add u3
            case 3 % assume [u1 u2 u3] was given, do nothing
            otherwise % [u1 u2 u3 phi1 phi2] may have been given
                u = u(:,1:3) ; % remove [phi1 phi2]
        end
    % Deformed shape
        p.Vertices = padarray(mesh.X,[0 size(u,2)-size(mesh.X,2)],0,'post') + u ;
    % Color data
        p.FaceVertexCData = sqrt(sum(abs(u).^2,2)) ;
    end


end