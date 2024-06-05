function K = wavenumbers(PLATE,w,theta)
% Computes the wavenumber spectrums of a thick plate over a range of
% frequencies w and propagation angles theta

% Thin plate bending
K.kb = sqrt(w).*(PLATE.m./PLATE.D(1,1)).^.25 ;

% Membrane waves
K.kl = w.*sqrt(PLATE.m./PLATE.A(1,1)) ;
K.kt = w.*sqrt(PLATE.m./PLATE.A(3,3)) ;



end

