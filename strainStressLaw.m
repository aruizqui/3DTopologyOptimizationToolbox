function C = strainStressLaw(E, v)

% Generate stress-Strain law matrix from Hooke's law in 3d
% Inputs:
%   - E: Young's Modulus [MPa].
%   - v: Poisson ratio (non-dimensional).
% v1: Feb 2021 - Vicente Cholvi Gil
% v2: Dec 2024 - A. Ruiz: Change matrix C sintaxis to improve legibility.

ci = E/((1+v)*(1-2*v));

C = ci * [
    1-v, v, v, 0, 0, 0;
    v, 1-v, v, 0, 0, 0;
    v, v, 1-v, 0, 0, 0;
    0, 0, 0, (1-2*v)/2, 0, 0;
    0, 0, 0, 0, (1-2*v)/2, 0;
    0, 0, 0, 0, 0, (1-2*v)/2
];

