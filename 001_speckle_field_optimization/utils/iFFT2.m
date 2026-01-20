function F_E = iFFT2(E)
    F_E = fftshift(fftshift(ifft2(ifftshift(ifftshift(E,1),2)),1),2) * sqrt(length(E(:,1,1))*length(E(1,:,1)));
end