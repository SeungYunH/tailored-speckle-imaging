function F_E = FFT2(E)
    F_E = fftshift(fftshift(fft2(ifftshift(ifftshift(E,1),2)),1),2) / sqrt(length(E(:,1,1))*length(E(1,:,1)));
end