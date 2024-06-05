%% SYNTHESIZE A MEASUREMENT ON A RIBBON
clc
clear all

%% RECTANGULAR, ISOTROPIC, HOMOGENEOUS PLATE
clc,clearvars ;
PLATE = [] ;
% Plate parameters
    PLATE.L = [100 80] ; % plate dimensions (mm)
    PLATE.dx = 2*[1 1] ; % spatial step (mm)
    PLATE.h = .1 ; % plate thickness (mm)
    PLATE.E = 70e3.*(1+.05i) ; % Young Modulus (MPa)
    PLATE.nu = .33 ; % Poisson ratio
    PLATE.rho = 2600e-12 ; % Density (tons/mm3)
    PLATE.bc = 'clamped' ; % kind of boundary conditions
    PLATE.xf = PLATE.L.*[3/10 4/10] ; % location of the point load application
    w = 2*pi*15e3 ; % Frequency (rad/s) 
% Basic discretization
    PLATE.mesh = fe.mesh.gridmesh(PLATE.L,PLATE.dx) ;
% Plate behavior
    PLATE = plate.behavior(PLATE) ;
% Solve using FEM
    u = fe.mindlin(PLATE,w) ;
    U = reshape(u,[PLATE.mesh.nX+1 numel(w)]) ;
    X = reshape(PLATE.mesh.X,[PLATE.mesh.nX+1 2]) ;
% Spatial response
    clf ; axis equal tight
    %fe.mesh.plotmesh(PLATE.mesh,real(u(:,end))) ;
    surf(X(:,:,1),X(:,:,2),real(U),'edgecolor','k','facecolor','interp') ;
    title("Harmonic response at f="+string(w/2/pi)+"Hz")
    xlabel 'X(mm)' ; ylabel 'Y(mm)' ;
%% Fourier response
    [Fu,K] = spectrum.fourier(U,PLATE.dx) ;
    clf ; axis equal tight ;
    surf(K{1},K{2},0*abs(Fu),abs(Fu),'edgecolor','none','facecolor','interp') ;
    title("Wavevector spectrum at f="+string(w/2/pi)+"Hz")
    xlabel 'Wavenumber kx (rad/mm)' ; ylabel 'Wavenumber ky (rad/mm)' ;
% theoretical wavenumbers 
    theta = linspace(0,2*pi,100)' ;
    Kp = plate.wavenumbers(PLATE,w,theta) ;
    plot(real(Kp.kb).*cos(theta),real(Kp.kb).*sin(theta),':r') ;
    
%% ESPRIT

out = ESPRIT_fcn(U,'DIMS_K',[1 2],'R0',5) ;
    clf ; axis equal tight ;
    surf(K{1},K{2},0*abs(Fu),abs(Fu),'edgecolor','none','facecolor','interp') ;
    title("Wavevector spectrum at f="+string(w/2/pi)+"Hz")
    xlabel 'Wavenumber kx (rad/mm)' ; ylabel 'Wavenumber ky (rad/mm)' ;
    ke = out.K.'./PLATE.dx ;
    plot(real(ke(:,1)),real(ke(:,2)),'.w','markersize',30) ;        
         