clear; clc; close all;

% Define Nodes and Grid Spacing
n = 51; % Number of circumferential nodes
m = 52; % Number of axial nodes
dx = 2 * pi / (n - 1);
dy = 1 / (m - 1);
X = linspace(0, 2*pi, n); % Angular grid
Y = linspace(0, 1, m);    % Axial grid

% Define Geometry Parameters
lambdaValues = [0, 1, 2, 4]; % Length-to-width ratios
eccentricityValues = linspace(0, 1.0, 21); % Eccentricities from 0 to 1 (21 values)

% Plot Styling
lineStyles = {'b-', 'r--', 'k-.', 'g:'};

% Store Results
Results(length(lambdaValues)) = struct( ...
    'BearingNumber', [], ...
    'MinFilmThickness', [], ...
    'AttitudeAngle', [], ...
    'FrictionCoefficient', [], ...
    'FlowRate', [], ...
    'LeakageRatio', [], ...
    'MaxPressureRatio', [], ...
    'MaxPressureAngle', [], ...
    'TerminatingPressureAngle', []);

% Pressure and Film Thickness Initialization
P = zeros(n, m); % Pressure grid
H = zeros(n, m); % Film thickness grid

% Solver Parameters
relFactor = 1.8; % Relaxation factor for SOR
tol = 1e-6;      % Convergence tolerance
maxIter = 5000;  % Maximum number of iterations

fprintf('Starting Reynolds equation solver...\n');

% Loop Over Lambda Values
for lambdaIdx = 1:length(lambdaValues)
    lambda = lambdaValues(lambdaIdx);
    geomFactor = (lambda / 2)^2;

    % Loop Over Eccentricity Values
    for eccIdx = 1:length(eccentricityValues)
        ecc = eccentricityValues(eccIdx);

        % Update Film Thickness
        for i = 1:n
            for j = 1:m
                H(i, j) = 1 + ecc * cos(X(i));
                if H(i, j) < 0 % Boundary condition
                    H(i, j) = 0;
                end
            end
        end

        % Initialize error and iteration count
        ERR = 1.0;
        iter = 0;

        % Iterate Using SOR
        while ERR > tol && iter < maxIter
            ERR = 0;
            iter = iter + 1;

            for i = 2:n-1
                for j = 2:m-1
                    hR = (H(i, j) + H(i+1, j)) / 2;
                    hL = (H(i, j) + H(i-1, j)) / 2;
                    hUp = geomFactor * H(i, j)^3 * (dx / dy)^2;

                    sourceTerm = -6 * ecc * sin(X(i));
                    sorUpdate = (hR^3 * P(i+1, j) + hL^3 * P(i-1, j) + ...
                                 hUp * P(i, j+1) + hUp * P(i, j-1) - sourceTerm) ...
                                 / (hR^3 + hL^3 + 2 * hUp);

                    newPressure = (1 - relFactor) * P(i, j) + relFactor * sorUpdate;
                    P(i, j) = max(0, newPressure); % Clamp pressure
                    ERR = ERR + abs(P(i, j) - newPressure);
                end
            end
        end

        % Compute Derived Quantities
        Wx = sum(sum(P .* sin(X') * dx * dy));
        Wz = sum(sum(-P .* cos(X') * dx * dy));
        resultantForce = sqrt(Wx^2 + Wz^2);

        if resultantForce < 1e-9; resultantForce = 1e-9; end

        % Bearing Number
        Bj = 1 / (pi * resultantForce);

        % Attitude Angle
        phiAtt = atan2d(Wx, Wz);

        % Friction Coefficient
        frictionForce = 0;
        for i = 1:n-1
            dP_dTheta = (P(i+1, :) - P(i, :)) / dx;
            shearStress = (H(i, :) / 2) .* dP_dTheta + (1 ./ H(i, :));
            frictionForce = frictionForce + sum(shearStress) * dx * dy;
        end
        frictionCoeff = frictionForce / resultantForce;

        % Flow Rate
        dP_dx = (P(2, :) - P(1, :)) / dx;
        fluxIn = (H(1, :) / 2) - (H(1, :).^3 / 12 .* dP_dx);
        qCircum = sum(fluxIn) * dy;
        QVar = 2 * pi * qCircum;

        % Leakage Ratio
        leakageSum = 0;
        for i = 1:n
            dP_dy_1 = (P(i, 2) - P(i, 1)) / dy;
            dP_dy_2 = (P(i, m) - P(i, m-1)) / dy;
            leakageSum = leakageSum + ((H(i, 1)^3 / 12) * dP_dy_1 + abs(H(i, m)^3 / 12 * dP_dy_2)) * dx;
        end
        leakageRatio = geomFactor * leakageSum / qCircum;

        % Pressure Metrics
        [maxPressure, pressureIdx] = max(P(:));
        [rowMax, ~] = ind2sub(size(P), pressureIdx);
        maxPressureAngle = X(rowMax) * 180 / pi;

        Results(lambdaIdx).BearingNumber(end+1) = Bj;
        Results(lambdaIdx).MinFilmThickness(end+1) = 1 - ecc;
        Results(lambdaIdx).AttitudeAngle(end+1) = phiAtt;
        Results(lambdaIdx).FrictionCoefficient(end+1) = frictionCoeff;
        Results(lambdaIdx).FlowRate(end+1) = QVar;
        Results(lambdaIdx).LeakageRatio(end+1) = leakageRatio;
        Results(lambdaIdx).MaxPressureRatio(end+1) = resultantForce / (2 * maxPressure);
        Results(lambdaIdx).MaxPressureAngle(end+1) = maxPressureAngle - phiAtt;
        Results(lambdaIdx).TerminatingPressureAngle(end+1) = maxPressureAngle + phiAtt - 180;
    end
end

fprintf('Main calculations completed.\n');

% Generate Figures (1–8)
figure(2); hold on;
for k = 1:length(lambdaValues)
    semilogx(Results(k).BearingNumber, Results(k).MinFilmThickness, lineStyles{k}, 'LineWidth', 1.5);
end
title('Figure 2: Minimum Film Thickness vs Bearing Number');
xlabel('Bearing Number, B_j');
ylabel('Dimensionless Minimum Film Thickness, h_{min}');
xlim([0.01 100]);
ylim([0 1]);
grid on;
legend('\lambda=0', '\lambda=1', '\lambda=2', '\lambda=4', 'Location', 'SouthEast');

figure(3); hold on;
for k = 1:length(lambdaValues)
    semilogx(Results(k).BearingNumber, Results(k).AttitudeAngle, lineStyles{k}, 'LineWidth', 1.5);
end
title('Figure 3: Attitude Angle vs Bearing Number');
xlabel('Bearing Number, B_j');
ylabel('Attitude Angle, \phi (deg)');
xlim([0.01 100]);
ylim([0 100]);
grid on;
legend('\lambda=0', '\lambda=1', '\lambda=2', '\lambda=4', 'Location', 'SouthEast');

% Repeat similar plotting commands for Figures 4–8...