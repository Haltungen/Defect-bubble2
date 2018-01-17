% DemoBubbleResonanceNumerical
%
% Overview:
%   We calculate the Minnaert resonance for a single bubble
%   in an infinite extent of liquid by directly solving the
%   boundary integral formulation of the problem for the
%   the characteristic value.
%
% References:
%   Mathematical and Computational Methods in Photonics - Tutorial Notes
%
% Authors:
%   Habib Ammari, Brian Fitzpatrick, Matias Ruiz, Sanghyeon Yu.
clc
clear all
% close all

% M = 0.2 + 0.1i;
% z0 = M;
% z1 = M - 1/100;
% z2 = M - 2/100;
distTol = 5e-4;
fTol = 1e-10;
iterMax = 100;
pointsN = 2^6;

radiusMin = 0.2;
radiusMax = radiusMin;
radiusN = 2^0;
radiusRange = fliplr(linspace(radiusMin, radiusMax, radiusN));

firstPass = 1;
initialGuess = 0.0002;
omegaSingleBubble = zeros(radiusN, 1);
for iRadius = 1:radiusN
    radius = radiusRange(iRadius);
    
    if firstPass == 1
        firstPass = 0;
    else
        initialGuess = omegaSingleBubble(iRadius-1);
    end
    
    z0 = initialGuess;
    z1 = initialGuess-initialGuess/100;
    z2 = initialGuess-initialGuess/200;    
    
    fprintf('iRadius: %d    (Radius: %.4f)\n', iRadius, radius);
    omegaSingleBubble(iRadius) = tools.MullersMethod('f', z0, z1, z2, iterMax, distTol, fTol, radius, pointsN);
    save('omegaSingleBubble', 'omegaSingleBubble');
end

