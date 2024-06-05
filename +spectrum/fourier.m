function [Fu,K] = fourier(U,dx)
% 2D Wavevector spectrum given by the Fourier Transform
nX = size(U,1:2) ;

K = arrayfun(@linspace,-pi./dx,pi./dx,nX,'uni',false) ;
[K{:}] = ndgrid(K{:}) ;

Fu = fftshift(fft(fft(U,[],1),[],2)) ;
end

