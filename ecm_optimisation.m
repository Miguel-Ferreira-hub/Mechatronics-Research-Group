% This script fits parameters to the equivalent circuit model (ECM) of the
% battery by solving the optimisation problem, and estimates the state of
% charge (SoC) through the predicted open circuit voltage (OCV).
close all;

% Plot ECM fit graphs if set to True
display=true;

%% Optimise for different cells by changing "cell" to "12","10","1" or "old"
% old cell is the original cell that is 100% SoH (labelling of cells comes
% from original position on test rig

cell = "12";

%% !! The code does not automatically change the director, so for each cell the MATLAB directory
%% must be changed to where the data is stored for each cell !!

switch cell

    case "12"
        % Cell 12
        data_list = [85,80,75,70,65,60,55,50,40,20,15,10,05];
        labels = [96.00,90.35,84.71,79.06,73.41,67.76,62.12,56.47,45.17,22.59,16.94,11.29,5.65];
        SoC = [5.65,11.29,16.94,22.59,45.17,56.47,62.12,67.76,73.41,79.06,84.71,90.35,96.00];
        Vt = [3.361,3.453,3.489,3.528,3.665,3.727,3.770,3.821,3.864,3.905,3.946,3.989,4.033];
        V_array = [3.387,3.473,3.508,3.547,3.685,3.748,3.791,3.841,3.882,3.923,3.965,4.008,4.052];
        I = 210e-3;

    case "10"
        % Cell 10
        data_list = [85,70,65,60,55,45,40,35,30,20,10,05];
        labels = flip([6.00,12.00,24.00,36.01,42.01,48.01,54.01,66.01,72.01,78.01,84.01,99.72]);
        SoC = flip([6.00,12.00,24.00,36.01,42.01,48.01,54.01,66.01,72.01,78.01,84.01,99.72]);
        Vt = [3.403,3.465,3.559,3.625,3.643,3.661,3.684,3.796,3.843,3.884,3.925,4.069];
        V_array = [3.450,3.503,3.589,3.656,3.673,3.693,3.716,3.829,3.874,3.915,3.958,4.100];
        I = 210e-3;

    case "1"
        % Cell 1
        data_list = flip([85,80,75,70,60,55,50,40,35,30,25,20,15,05]);
        labels = flip([98.91,93.09,87.27,81.45,69.81,64.00,58.18,46.54,40.72,34.91,29.09,23.27,17.45,5.82]);
        SoC = [5.82,17.45,23.27,29.09,34.91,40.72,46.54,58.18,64.00,69.81,81.45,87.27,93.09,98.91];
        Vt = [3.380,3.522,3.566,3.599,3.617,3.633,3.648,3.700,3.736,3.796,3.913,3.952,3.992,4.035];
        V_array = [3.423,3.550,3.593,3.626,3.643,3.660,3.676,3.731,3.766,3.824,3.940,3.979,4.020,4.064];
        I = 210e-3;

    case "old"
        % 100% SoH cell
        data_list = [100, 90, 80, 70, 60, 50, 40, 30, 20, 00];
        labels = flip([100.0, 90.7, 81.4, 72.0, 62.7, 53.4, 44.1, 34.8, 25.4, 0.00]);
        Vt = [2.75,3.394,3.480,3.535,3.566,3.599,3.654,3.769,3.849,4.008];
        I = 419e-3;

end

OCV_original = V_array;

%% Fit both 2RC and 5RC models

n_data = length(data_list);

% Number of parameters:
x_opts_2RC = zeros(n_data,5);
x_opts_5RC = zeros(n_data,13);

% Store impedance data and fitted impedance
Z_original = cell(n_data,1);
Z_2RC = cell(n_data,1);
Z_5RC = cell(n_data,1);

% Store frequency
frequency_store = cell(n_data,1);

% Store RMSE
RMSE_2RC = zeros(n_data,2);
RMSE_5RC = zeros(n_data,2);

for i = 1:n_data

    %% Read data

    filename = sprintf('%02d.txt', data_list(i));

    data = readtable(filename);

    freq = data.freq__Hz;
    ZI = data.Z1_ohm;
    ZII = data.Z2_ohm;

    % Remove positive imaginary impedance
    mask = ZII < 0;

    freq = freq(mask);
    ZI = ZI(mask);
    ZII = ZII(mask);

    % Smoothed impedance used for fitting
    Z = sgolayfilt(ZI + 1i*ZII, 3, 11);

    label = sprintf('%.1f%%', labels(i));

    %% Store original data

    Z_original{i} = Z;
    frequency_store{i} = freq;

    %% 2RC FIT

    [x2, rmse2I, rmse2II, Zfit2] = ...
        fit_ECM(Z, freq, '2RC');

    x_opts_2RC(i,:) = x2;

    Z_2RC{i} = Zfit2;

    RMSE_2RC(i,:) = [rmse2I rmse2II];

    %% 5RC FIT

    [x5, rmse5I, rmse5II, Zfit5] = ...
        fit_ECM(Z, freq, '5RC');

    x_opts_5RC(i,:) = x5;

    Z_5RC{i} = Zfit5;

    RMSE_5RC(i,:) = [rmse5I rmse5II];

    %% Display results

    fprintf('\nSoC = %s\n', label);

    fprintf('2RC RMSE real = %.6e\n', rmse2I);
    fprintf('2RC RMSE imag = %.6e\n', rmse2II);

    fprintf('5RC RMSE real = %.6e\n', rmse5I);
    fprintf('5RC RMSE imag = %.6e\n', rmse5II);

end

%% OCV estimation

OCV_2RC = nan(size(Vt));
OCV_5RC = nan(size(Vt));

for i = 1:length(Vt)

    %% 2RC OCV
    OCV_2RC(i) = Vt(i) + I * sum(x_opts_2RC(i,1:3));


    %% 5RC OCV
    OCV_5RC(i) = Vt(i) + I * ...
        (sum(x_opts_5RC(i,1:6)) + x_opts_5RC(i,12));

end

%% Interpolate
SoC_interp_2RC = interp1(V_array, SoC, OCV_2RC, 'spline'); 
SoC_interp_5RC = interp1(V_array, SoC, OCV_5RC, 'spline'); 

x_interp = linspace(0,100,500);

% Original OCV curve
OCV_interp_original = spline(SoC, OCV_original, x_interp);
OCV_interp_2RC = spline(SoC_interp_2RC, OCV_2RC, x_interp);
OCV_interp_5RC = spline(SoC_interp_5RC, OCV_5RC, x_interp);

%% Plot OCV curves
figure('Name','OCV Comparison');

plot(SoC, OCV_original, ...
    'o', ...
    'DisplayName','Original OCV', ...
    'LineWidth',1);

hold on

plot(x_interp, OCV_interp_2RC, ...
    'DisplayName','2RC OCV', ...
    'LineWidth',1);

plot(x_interp, OCV_interp_5RC, ...
    'DisplayName','5RC OCV', ...
    'LineWidth',1);

hold off

xlabel('$SoC\ (\%)$', ...
    'Interpreter','latex', ...
    'FontSize',15);

ylabel('$OCV\ (V)$', ...
    'Interpreter','latex', ...
    'FontSize',15);

legend('Location','northwest', ...
    'FontSize',10, ...
    'Box','off');

xlim([0 100]);

box on

%% Nyquist plots - individual SoCs
for i = 1:n_data

    figure('Name', sprintf('Nyquist - SoC %.1f%%', labels(i)));

    plot(real(Z_original{i}), -imag(Z_original{i}), ...
        'o', ...
        'DisplayName','Data', ...
        'LineWidth',1);

    hold on

    plot(real(Z_2RC{i}), -imag(Z_2RC{i}), ...
        'DisplayName','2RC', ...
        'LineWidth',1);

    plot(real(Z_5RC{i}), -imag(Z_5RC{i}), ...
        'DisplayName','5RC', ...
        'LineWidth',1);

    hold off

    title(sprintf('Nyquist Plot - SoC = %.1f%%', labels(i)));

    xlabel('$Z''\ (\Omega)$', ...
        'Interpreter','latex');

    ylabel('$-Z''''\ (\Omega)$', ...
        'Interpreter','latex');

    legend('Location','best', ...
        'FontSize',8, ...
        'Box','off');

    box on

end

%% Bode magnitude plots - individual SoCs

for i = 1:n_data

    figure('Name', sprintf('Bode Magnitude - SoC %.1f%%', labels(i)));

    semilogx(frequency_store{i}, ...
        M(real(Z_original{i}),imag(Z_original{i})), ...
        'o', ...
        'DisplayName','Data', ...
        'LineWidth',1);

    hold on

    semilogx(frequency_store{i}, ...
        M(real(Z_2RC{i}),imag(Z_5RC{i})), ...
        'DisplayName','2RC', ...
        'LineWidth',1);

    semilogx(frequency_store{i}, ...
        M(real(Z_5RC{i}),imag(Z_5RC{i})), ...
        'DisplayName','5RC', ...
        'LineWidth',1);

    hold off

    title(sprintf('Bode Magnitude - SoC = %.1f%%', labels(i)));

    xlabel('$Frequency\ (Hz)$', ...
        'Interpreter','latex');

    ylabel('$|Z|\ (dB)$', ...
        'Interpreter','latex');

    legend('Location','best', ...
        'FontSize',8, ...
        'Box','off');

    box on

end

%% Bode phase plots - individual SoCs

for i = 1:n_data

    figure('Name', sprintf('Bode Phase - SoC %.1f%%', labels(i)));

    semilogx(frequency_store{i}, ...
        P(real(Z_original{i}),imag(Z_original{i})), ...
        'o', ...
        'DisplayName','Data', ...
        'LineWidth',1);

    hold on

    semilogx(frequency_store{i}, ...
        P(real(Z_2RC{i}),imag(Z_2RC{i})), ...
        'DisplayName','2RC', ...
        'LineWidth',1);

    semilogx(frequency_store{i}, ...
        P(real(Z_5RC{i}),imag(Z_5RC{i})), ...
        'DisplayName','5RC', ...
        'LineWidth',1);

    hold off

    title(sprintf('Bode Phase - SoC = %.1f%%', labels(i)));

    xlabel('$Frequency\ (Hz)$', ...
        'Interpreter','latex');

    ylabel('$Phase\ (deg)$', ...
        'Interpreter','latex');

    legend('Location','best', ...
        'FontSize',8, ...
        'Box','off');

    box on

end

% Function to solve optimisation problem and fit ECM
function [x_opt, RMSEI, RMSEII, Zfit] = ...
    fit_ECM(Zdata, frequencydata, model_type)

    %% Frequency

    frequency = frequencydata(:);
    omega = frequency * 2*pi;

    %% Select model

    switch model_type

        case '2RC'

            % -------------------------------------------------------------
            % 2RC MODEL
            %
            % x = [R0 R1 R2 C1 C2]
            % -------------------------------------------------------------

            x0 = [ ...
                0.01, ...      % R0
                0.001, ...     % R1
                0.001, ...     % R2
                0.01, ...      % C1
                0.01];        % C2

            lb = [ ...
                0, 0, 0, 0, 0];

            ub = [ ...
                1, 1, 1, 1000, 1000];

        case '5RC'
            % 5RC MODEL

            x0 = [0.01,0.001,0.001,0.02,0.002,0.002,0.01,0.01,0.01,1,7,0.002,0.1]; 
            lb = [0,0,0,0,0,0,0,0,0,0,0,0,0]; 
            ub = [1,1,1,1,1,1,10,10,10,10,2000,200,1000];

        otherwise
            error('Unknown ECM model Specified.');
    end

    %% Objective function

    objective = @(x) [ ...
        real(ECM_impedance(x,omega,model_type)) - real(Zdata);
        imag(ECM_impedance(x,omega,model_type)) - imag(Zdata)
        ];

    %% Optimisation options

    opts = optimoptions('lsqnonlin', ...
        'Display','iter', ...
        'FunctionTolerance',1e-12, ...
        'OptimalityTolerance',1e-12);

    problem = createOptimProblem('lsqnonlin', ...
        'x0',x0, ...
        'objective',objective, ...
        'lb',lb, ...
        'ub',ub, ...
        'options',opts);

    %% MultiStart for best solution (initialised from random strating points (initial guesses))

    ms = MultiStart('UseParallel',true);

    [x_opt,~] = run(ms,problem,200);

    %% Fitted impedance

    Zfit = ECM_impedance(x_opt,omega,model_type);

    %% RMSE

    RMSEI = sqrt(mean( ...
        (real(Zdata)-real(Zfit)).^2));

    RMSEII = sqrt(mean( ...
        (imag(Zdata)-imag(Zfit)).^2));

end

% Function for impedance model to be optimised
function Z = ECM_impedance(x,omega,model_type)

    %% Ohmic resistance

    R0 = x(1);

    %% 2RC

    if strcmp(model_type,'2RC')

        R1 = x(2);
        R2 = x(3);

        C1 = x(4);
        C2 = x(5);

        Z1 = R1 ./ ...
            (1 + 1i*omega*R1*C1);

        Z2 = R2 ./ ...
            (1 + 1i*omega*R2*C2);

        Z = R0 + Z1 + Z2;

    %% 5RC

    elseif strcmp(model_type,'5RC')
        R0 = x(1);  
        R1 = x(2);  
        R2 = x(3);
        R3 = x(4);  
        R4 = x(5);  
        R5 = x(6);
        C1 = x(7);  
        C2 = x(8);
        C3 = x(9);  
        C4 = x(10);  
        C5 = x(11);
        RD   = x(12);
        tauD = x(13);
    
        % RC elements in parallel
        Z1 = R1 ./ (1 + 1i * omega * R1 * C1);
        Z2 = R2 ./ (1 + 1i * omega * R2 * C2);
        Z3 = R3 ./ (1 + 1i * omega * R3 * C3);
        Z4 = R4 ./ (1 + 1i * omega * R4 * C4);
    
        % (Warburg + R5) in parallel with C5
        ZC5 = -1i./(omega*C5);
        ZW = RD .* tanh(sqrt(1i*omega*tauD)) ./ sqrt(1i*omega*tauD);
        Z5 = 1 ./ (1 ./ ZC5 + 1 ./ (R5 + ZW));
    
        % Total Impedance
        Z = R0 + Z1 + Z2 + Z3 + Z4 + Z5;

    end

end

% Magnitude of impedance
function magnitude=M(ZI, ZII) 
    magnitude=20*log10(((ZI.^2+ZII.^2).^(1/2)));
end

% Phase of impedance
function phase=P(ZI, ZII) 
    phase=atand(ZII./ZI);
end