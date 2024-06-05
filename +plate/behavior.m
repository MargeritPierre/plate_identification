function PLATE = behavior(PLATE)
% Computes the stiffness & interia matrices of a given plate
PLATE.xi = pi^2./12 ; % shear correction factor

% Plate stiffnesses
    PLATE.Q = (PLATE.E./(1-PLATE.nu.^2)).*[1 PLATE.nu 0 ; PLATE.nu 1 0 ; 0 0 .5*(1-PLATE.nu)] ; % plane-stress stiffness
    PLATE.A = PLATE.h.*PLATE.Q ; % membrane stifness
    PLATE.D = PLATE.h.^3/12.*PLATE.Q ; % bending stiffness
    PLATE.F = PLATE.xi*PLATE.h.*PLATE.E./2./(1+PLATE.nu).*[1 0 ; 0 1] ; % out-of-plane shear stiffness
% Plate Inertia
    PLATE.m = PLATE.rho.*PLATE.h ; % section interia
    PLATE.J = PLATE.rho.*PLATE.h.^3./12 ; % rotary inertia
    
end

