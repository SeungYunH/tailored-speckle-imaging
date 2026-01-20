% example_pdf_generation.m
addpath('../utils');

SLMsize = [800 600];

% Super-Rayleigh (order 3)
sup3 = generateSuperRayleighPDF(3, SLMsize);
% Rayleigh
R    = generateSuperRayleighPDF(1, SLMsize);
% Binary
Bin  = generateBinaryPDF(SLMsize, 1);
% Delta (sub-Rayleigh)
Dt   = generateDeltaPDF(SLMsize);

% Plot them
figure; hold on;
hist_err(sup3.pdf, sup3.name);
hist_err(R.pdf, 'R');
hist_err(Bin.pdf, 'Bin1');
hist_err(Dt.pdf, 'Dt');
legend;
title('Example Target PDFs');