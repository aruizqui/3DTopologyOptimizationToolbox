function [Ki,B_accum, VOL] = hex3DShapeFunction(nc,C)
    % Calculate the B matrix, Jacobian, and volume for a HEX8 element
    % Adrian Ruiz Quiñones Dec 2024
    %
    % Inputs:
    %   nc - 8x3 matrix with global coordinates of the 8 nodes of the element
    %
    % Outputs:
    %   B       - Strain-displacement matrix (6x24)
    %   B_avg   - B matrix averaged in all Gauss Points (6x24)
    %   J       - Jacobian matrix (3x3)
    %   VOL     - Volume of the HEX8 element

    [gausepoint, weights] = gaussPointsHex8();
    numGP = size(gausepoint, 1);
    betaLength = size(nc,1);
    Ki = zeros (24,24);
    B_accum = zeros(6,24);
    VOL = 0;
    
    for k = 1:numGP
        xi = gausepoint(k, 1); eta = gausepoint(k, 2); zeta = gausepoint(k, 3);
        w = weights(k);
        RSTderivation = 1/8 * [-(1-eta)*(1-zeta),  (1-eta)*(1-zeta),...
            (1+eta)*(1-zeta), -(1+eta)*(1-zeta), ...
            -(1-eta)*(1+zeta),  (1-eta)*(1+zeta),...
            (1+eta)*(1+zeta), -(1+eta)*(1+zeta);

            -(1-xi)*(1-zeta), -(1+xi)*(1-zeta),...
            (1+xi)*(1-zeta),  (1-xi)*(1-zeta), ...
            -(1-xi)*(1+zeta), -(1+xi)*(1+zeta),...
            (1+xi)*(1+zeta),  (1-xi)*(1+zeta);

            -(1-xi)*(1-eta), -(1+xi)*(1-eta),...
            -(1+xi)*(1+eta), -(1-xi)*(1+eta), ...
             (1-xi)*(1-eta),  (1+xi)*(1-eta),...
             (1+xi)*(1+eta),  (1-xi)*(1+eta)
        ];


        J = nc'*RSTderivation';
        invJ = inv(J');
        XYZderivation = invJ * RSTderivation;
        detJ = det(J');

        if detJ <= 0
            error('Negative or zero Jacobian determinant at Gauss point %d', k);
        end

        B = zeros(6, 3*betaLength); % B matrix (accumulates contributions from all Gauss points)
        
        % Epsilon x y z 
        B(1,1:betaLength)= XYZderivation(1,:); % Derivadas respecto de x
        B(2,betaLength+1:2*betaLength)=XYZderivation(2,:); % Derivadas respecto a y
        B(3, 2*betaLength+1:3*betaLength) = XYZderivation(3, :);  % Derivadas respecto a z
        % Gamma xy
        B(4, 1:betaLength) = XYZderivation(2, :);  % Derivadas respecto a y
        B(4, betaLength+1:2*betaLength) = XYZderivation(1, :);  % Derivadas respecto a x
        % Gamma yz
        B(5, betaLength+1:2*betaLength)= XYZderivation(3,:);
        B(5,2*betaLength+1:3*betaLength) = XYZderivation(2,:);
        % Gamma xz
        B(6, 1:betaLength)= XYZderivation(3,:);
        B(6,2*betaLength+1:3*betaLength) = XYZderivation(1,:);
        
        B_accum = B_accum + B * w * detJ/1e3;
        VOL = VOL + w * det(J)/1e3;
        Ki = Ki +  B' * C * B* w * detJ;
        
    end
    
end

function [gp,w] = gaussPointsHex8()
        % gaussPointsHex8: Provides Gauss points and weights for HEX8 elements
        % Outputs:
        %   gp - [8 x 3] array, each row contains (ξ, η, ζ) coordinates of Gauss points
        %   w  - [8 x 1] array, weights for each Gauss point
        % Gauss points in 1D and weights
        g1D = [-sqrt(1/3), sqrt(1/3)]; % Gauss points for 2-point quadrature in 1D
        w1D = [1, 1];                  % Corresponding weights in 1D

        % Generate 3D Gauss points and weights
        [xi,eta,zeta] =ndgrid(g1D,g1D,g1D);
        gp= [xi(:),eta(:),zeta(:)];
        [wx, wy, wz] = ndgrid(w1D, w1D, w1D);
        w = wx(:) .* wy(:) .* wz(:);
end

