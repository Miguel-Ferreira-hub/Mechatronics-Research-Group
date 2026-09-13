%% DRT for four different states of health

close all;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell / SoH information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cells(1).name = '100.0% SoH';
cells(1).directory = ...
    "C:\Users\migue\Desktop\Part IV\Quant\Up to Date ML for EIS\Old Data";
cells(1).data_list = [100, 90, 80, 70, 60, 50, 40, 30, 20, 00];
cells(1).labels = [100.0, 90.7, 81.4, 72.0, 62.7, 53.4, ...
                   44.1, 34.8, 25.4, 0.00];

cells(2).name = '85.9% SoH';
cells(2).directory = ...
    "C:\Users\migue\Desktop\Part IV\Quant\Up to Date ML for EIS\Cell 1";
cells(2).data_list = [85,80,75,70,60,55,50,40,35,30,25,20,15,05];
cells(2).labels = [98.91,93.09,87.27,81.45,69.81,64.00, ...
                   58.18,46.54,40.72,34.91,29.09,23.27, ...
                   17.45,5.82];

cells(3).name = '88.6% SoH';
cells(3).directory = ...
    "C:\Users\migue\Desktop\Part IV\Quant\Up to Date ML for EIS\Cell 12";
cells(3).data_list = [85,80,75,70,65,60,55,50,40,20,15,10,5,0];
cells(3).labels = flip([0.00,5.65,11.29,16.94,22.59,45.17, ...
                        56.47,62.12,67.76,73.41,79.06,84.71, ...
                        90.35,96.00]);

cells(4).name = '83.3% SoH';
cells(4).directory = ...
    "C:\Users\migue\Desktop\Part IV\Quant\Up to Date ML for EIS\Cell 10";
cells(4).data_list = [85,70,65,60,55,45,40,35,30,20,10,05,00];
cells(4).labels = flip([0.00,6.00,12.00,24.00,36.01,42.01, ...
                        48.01,54.01,66.01,72.01,78.01,84.01,99.72]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target SoCs for Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_soc = [25, 50, 75, 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRT figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','DRT - Four States of Health');

tiledlayout(2,2);

% Same colours for the same target SoC in every subplot
colors = winter(length(target_soc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage for DRT peak and time constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peak = nan(length(cells), length(target_soc));
tau  = nan(length(cells), length(target_soc));

% Store the actual selected SoCs
selected_soc = nan(length(cells), length(target_soc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through each state of health
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c = 1:length(cells)

    % Create subplot
    ax_drt = nexttile;
    hold(ax_drt, 'on');

    % Find the four SoCs closest to 25, 50, 75 and 100%
    selected_indices = zeros(1,length(target_soc));

    for j = 1:length(target_soc)

        [~, selected_indices(j)] = ...
            min(abs(cells(c).labels - target_soc(j)));

        selected_soc(c,j) = cells(c).labels(selected_indices(j));

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate DRT for selected SoCs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j = 1:length(target_soc)

        idx = selected_indices(j);

        % Data filename
        filename = sprintf('%02d.txt', cells(c).data_list(idx));

        % Full path
        path = fullfile(cells(c).directory, filename);

        % Read data
        data = readtable(path);

        freq = data.freq__Hz;
        ZI   = data.Z1_ohm;
        ZII  = data.Z2_ohm;

        % Remove positive imaginary impedance points
        mask = ZII < 0;

        ZI   = ZI(mask);
        ZII  = ZII(mask);
        freq = freq(mask);

        % Smooth impedance data
        Z = sgolayfilt(ZI + 1i*ZII, 3, 11);

        % Actual SoC used
        actual_soc = cells(c).labels(idx);

        label = sprintf('%.1f%%', actual_soc);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DRT calculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [~, ~, freq_fine, gamma_ridge_fine] = ...
            DRTCalculation( ...
                real(Z), ...
                imag(Z), ...
                freq, ...
                'Gaussian', ...
                label, ...
                29.8e-6, ...
                0.5, ...
                ax_drt, ...
                colors(j,:) ...
            );

        % Find maximum DRT peak
        [max_gamma, peak_idx] = max(gamma_ridge_fine);

        peak(c,j) = max_gamma;

        % Corresponding relaxation time
        tau(c,j) = 1 ./ freq_fine(peak_idx);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Formatting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    title(ax_drt, cells(c).name);

    xlabel(ax_drt, '$\tau\ (s)$', ...
        'Interpreter','latex', ...
        'FontSize',15);

    ylabel(ax_drt, '$\gamma\ (\Omega)$', ...
        'Interpreter','latex', ...
        'FontSize',15);

    set(ax_drt, ...
        'XScale','log', ...
        'FontSize',12, ...
        'XTick',10.^(-10:2:10), ...
        'TickLabelInterpreter','latex');

    ylim(ax_drt,[0 0.035]);

    legend(ax_drt, ...
        'Location','northwest', ...
        'FontSize',9, ...
        'Box','off');

    box(ax_drt,'on');

    % Formatting of DRT plot

    if c == 1
        % Top-left DRT: all five labels
        text(ax_drt, 0.00096, 0.01, '$1$', ...
            'Interpreter', 'Latex', 'FontSize', 12);
    
        text(ax_drt, 0.015, 0.0027, '$2$', ...
            'Interpreter', 'Latex', 'FontSize', 12);
    
        text(ax_drt, 0.17, 0.0025, '$3$', ...
            'Interpreter', 'Latex', 'FontSize', 12);
    
        text(ax_drt, 3, 0.0040, '$4$', ...
            'Interpreter', 'Latex', 'FontSize', 12);
    
        text(ax_drt, 45, 0.003, '$5$', ...
            'Interpreter', 'Latex', 'FontSize', 12);
    
    else
        % Remaining three DRTs: labels 4 and 5 only
        text(ax_drt, 3, 0.0043, '$4$', ...
            'Interpreter', 'Latex', 'FontSize', 12);

        text(ax_drt, 45, 0.0030, '$5$', ...
            'Interpreter', 'Latex', 'FontSize', 12);
    end
    
    box(ax_drt, 'on');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display which SoCs were actually selected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nSelected SoCs:\n');
fprintf('------------------------------------------------------\n');

for c = 1:length(cells)

    fprintf('%s:\n', cells(c).name);

    for j = 1:length(target_soc)

        fprintf('  Target %.0f%% -> Selected %.2f%%\n', ...
            target_soc(j), selected_soc(c,j));

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display peak and tau values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nDRT peak values:\n');
disp(peak);

fprintf('\nDRT time constants:\n');
disp(tau);

% Formatting of Nyquist plot
xlabel(ax, '$ZI (\Omega)$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel(ax, '$-ZII (\Omega)$', 'Interpreter', 'Latex', 'FontSize', 15);
legend(ax, 'show', 'location', 'northwest', 'FontSize', 10, 'Box', 'off');
box(ax, 'on');

% Formatting of magnitude plot
xlabel(ax_mag, '$Frequency (Hz)$', 'Interpreter', 'Latex', 'FontSize', 15);
set(ax_mag, 'XScale', 'log')
ylabel(ax_mag, '$Magnitude (dB)$', 'Interpreter', 'Latex', 'FontSize', 15);
legend(ax_mag, 'show', 'location', 'northeast', 'FontSize', 10, 'Box', 'off');
box(ax_mag, 'on');

% Formatting of phase plot
xlabel(ax_phase, '$Frequency (Hz)$', 'Interpreter', 'Latex', 'FontSize', 15);
set(ax_phase, 'XScale', 'log')
ylabel(ax_phase, '$Phase (rad)$', 'Interpreter', 'Latex', 'FontSize', 15);
legend(ax_phase, 'show', 'location', 'east', 'FontSize', 10, 'Box', 'off');
box(ax_phase, 'on');

% Data for following plots (taken from previous plots) (old)
% SoC=[0, 25.4, 34.8, 44.1, 53.4, 62.7, 72, 81.4, 90.7, 100];
% tau=[62.0484, 54.1741, 52.7235, 52, 49.9377, 47.2991, 54.1741, 51.3117, 46.0325, 54.1741];
% peak=[0.19684, 0.02533, 0.023674, 0.019866, 0.01753, 0.017611, 0.035274, 0.029108, 0.020588, 0.023204];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell 1 - 100%
SoCo1 = [100.0, 90.7, 81.4, 72.0, 62.7, 53.4, 44.1, 34.8, 25.4, 0.00];
tauo1 = [54.1741,47.2991,52.7235,55.6647,47.2991,49.4753,52.7235,52.7235,54.1741,62.0484];
peako1 = [0.0232,0.0206,0.0291,0.0353,0.0176,0.0175,0.0199,0.0237,0.0253,0.1968];

% Cell 1 - 85.94%
SoC1 = [98.91,93.09,87.27,81.45,69.81,64.00,58.18,46.54,40.72,34.91,29.09,23.27,17.45,5.82];
tau1 = [50.8568,51.5573,51.5573,49.0316,50.8568,56.7503,54.7136,53.5859,50.8568,56.8117,55.6267,55.6942,53.8463,58.0228];
peak1 = [0.0196,0.0195,0.0191,0.0179,0.0267,0.0392,0.0323,0.0230,0.0200,0.0211,0.0191,0.0234,0.0248,0.0552];

% Cell 12 - 88.55%
SoC12 = flip([0.00,5.65,11.29,16.94,22.59,45.17,56.47,62.12,67.76,73.41,79.06,84.71,90.35,96.00]);
tau12 = [49.0316,49.0316,49.0316,50.8568,50.8568,52.0092,55.7482,54.7136,50.8568,52.7500,52.7500,52.0092,58.6061,65.8836];
peak12 = [0.0211,0.0212,0.0219,0.0244,0.0258,0.0298,0.0401,0.0341,0.0230,0.0258,0.0283,0.0302,0.0605,0.2439];

% Cell 10 - 83.30%
SoC10 = flip([0.00,12.00,24.00,30.00,36.01,42.01,48.01,54.01,66.01,72.01,78.01,84.01,99.72]);
tau10 = [49.6989,51.8337,51.8337,51.8337,54.0601,54.0601,51.8337,49.6989,49.6989,51.8337,56.3821,62.6717,66.8808];
peak10 = [0.0208,0.0242,0.0252,0.0269,0.0353,0.0293,0.0253,0.0225,0.0209,0.0250,0.0331,0.0736,0.1437];

colors = winter(4);

% Diffusion peak vs SoC plot
figure('Name', 'Diffusion peak vs SoC');
plot(SoCo1, peako1, 'linewidth', 1.5, 'DisplayName','100.0%','Color',colors(1,:));
hold on
plot(SoC12, peak12, 'linewidth', 1.5, 'DisplayName','88.55%','Color',colors(2,:));
plot(SoC1, peak1, 'linewidth', 1.5, 'DisplayName','85.94%','Color',colors(3,:));
plot(SoC10, peak10, 'linewidth', 1.5, 'DisplayName','83.30%','Color',colors(4,:));
xlabel('$SoC (\%)$', 'Interpreter', 'Latex','Fontsize',15);
ylabel('$\gamma (\Omega)$', 'Interpreter', 'Latex','Fontsize',15);
legend('location', 'northwest', 'FontSize', 10, 'Box', 'off');
box on

% Diffusion time constant vs SoC plot
figure('Name', 'Tau vs SoC');
plot(SoCo1, tauo1, 'linewidth', 1.5,'DisplayName','100.0%','Color',colors(1,:));
hold on
plot(SoC12, tau12, 'linewidth', 1.5, 'DisplayName','88.55%','Color',colors(2,:));
plot(SoC1, tau1, 'linewidth', 1.5, 'DisplayName','85.94%','Color',colors(3,:));
plot(SoC10, tau10, 'linewidth', 1.5, 'DisplayName','83.30%','Color',colors(4,:));
xlabel('$SoC (\%)$', 'Interpreter', 'Latex','Fontsize',15);
ylabel('$\tau (s)$', 'Interpreter', 'Latex','Fontsize',15);
% legend('location', 'northwest', 'FontSize', 10, 'Box', 'off');
box on

% Magnitude of impedance data
function magnitude=M(ZI, ZII)
    magnitude=20*log10(((ZI.^2+ZII.^2).^(1/2)));
end

% Phase of impedance data
function phase=P(ZI, ZII)
    phase=atand(ZII./ZI);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This portion of the code was built using the MATLAB script developed by
% 1. Ciucci et al. "DRTtools" with the following references:
% Wan, T. H., Saccoccio, M., Chen, C., & Ciucci, F. (2015). Influence of 
% the discretization methods on the distribution of relaxation times 
% deconvolution: implementing radial basis functions with DRTtools. 
% Electrochimica Acta, 184, 483-499.
% 2. M. Saccoccio, T. H. Wan, C. Chen, and F. Ciucci, "Optimal Regularization 
% in Distribution of Relaxation Times applied to Electrochemical Impedance 
% Spectroscopy: Ridge and Lasso Regression Methods - 
% A Theoretical and Experimental Study,” Electrochimica Acta, vol. 147, 
% pp. 470–482, Nov. 2014, doi:
% 3. Liu, J., Wan, T. H., & Ciucci, F. (2020).A Bayesian view on the 
% Hilbert transform and the Kramers-Kronig transform of electrochemical 
% impedance data: Probabilistic estimates and quality scores. 
% Electrochimica Acta, 357, 136864.**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute DRT
function [x_ridge, epsilon, freq_fine, gamma_ridge_fine] = ...
    DRTCalculation(ZI, ZII, freq, rbf_type, SoC, lambda, coeff, ax, colour)
    hold on
    taumax = ceil(max(log10(1./freq)))+0.5;    
    taumin = floor(min(log10(1./freq)))-0.5;
    freq_fine = logspace(-taumin, -taumax, 10*numel(freq));
%   bounds ridge regression
    lb = zeros(numel(freq)+2,1);
    ub = Inf*ones(numel(freq)+2,1);
    x_0 = ones(size(lb));
    options = optimset('algorithm','interior-point-convex','Display','off','TolFun',1e-15,'TolX',1e-10,'MaxFunEvals', 1E5);
    Z_exp=ZI+1i*ZII;
    b_re = real(Z_exp);
    b_im = imag(Z_exp);
%   compute epsilon
    shape_control=0.5;
    epsilon = compute_epsilon(freq, coeff, rbf_type, shape_control);
%   calculate the A_matrix
    A_re = assemble_A_re(freq, epsilon, rbf_type);
    A_im = assemble_A_im(freq, epsilon, rbf_type);

%   adding the resistence column to the A_re_matrix
    A_re(:,2) = 10e-3;
    
%   calculate the M_matrix
    %M = assemble_M_1(freq, epsilon, rbf_type);
    M = assemble_M_2(freq, epsilon, rbf_type);
    %lambda=10e-5;
    [H_combined,f_combined] = quad_format_combined(A_re, A_im, b_re, b_im, M, lambda);
    x_ridge = quadprog(H_combined, f_combined, [], [], [], [], lb, ub, x_0, options);
    
    % map x to gamma
    [gamma_ridge_fine,freq_fine] = map_array_to_gamma(freq_fine, freq, x_ridge(3:end), epsilon, rbf_type);
    freq_fine = freq_fine';
    
    plot(1./freq_fine, gamma_ridge_fine, 'DisplayName', SoC, 'LineWidth', 1);
    y_min = 0; 
    y_max = 0.035;
    xlabel('$\tau (s)$', 'Interpreter', 'Latex','Fontsize',15)
    ylabel('$\gamma (\Omega)$', 'Interpreter', 'Latex','Fontsize',15);
    set(gca,'xscale','log','xlim',[min(1./freq_fine), max(1./freq_fine)],'ylim',[y_min, 1.1*y_max],'Fontsize',15,'xtick',10.^(-10:2:10),'TickLabelInterpreter','latex');
    hold off
    legend('location', 'northwest', 'FontSize', 10, 'Box', 'off');
end

function out_A_re = assemble_A_re(freq, epsilon, rbf_type)

%   this function assembles the A_re matrix
%   we will assume that the tau vector is identical to the freq vector

%   first get number of frequencies
    N_freq = numel(freq);

%   the define the A_re output matrix
    out_A_re_temp = zeros(N_freq);
    out_A_re = zeros(N_freq, N_freq+2);

%   we compute if the frequencies are sufficiently log spaced
    std_diff_freq = std(diff(log(1./freq)));
    mean_diff_freq = mean(diff(log(1./freq)));

%   if they are, we apply the toeplitz trick
    toeplitz_trick = std_diff_freq/mean_diff_freq<0.01;

%   if terms are evenly distributed & do not use PWL
%   then we do compute only N terms
%   else we compute all terms with brute force
    if toeplitz_trick && ~strcmp(rbf_type,'Piecewise linear')

        % define vectors R and C
        R = zeros(1,N_freq);
        C = zeros(N_freq,1);

        % for clarity the C and R computations are separated
        for iter_freq_n = 1: N_freq

            freq_n = freq(iter_freq_n);
            freq_m = freq(1);
            C(iter_freq_n, 1) = g_i(freq_n, freq_m, epsilon, rbf_type);

        end

        for iter_freq_m = 1: N_freq

            freq_n = freq(1);
            freq_m = freq(iter_freq_m);
            R(1, iter_freq_m) = g_i(freq_n, freq_m, epsilon, rbf_type);

        end

        out_A_re_temp= toeplitz(C,R);

    else
    % compute using brute force

        for iter_freq_n = 1: N_freq

            for iter_freq_m = 1: N_freq

                freq_n = freq(iter_freq_n);
                freq_m = freq(iter_freq_m);

                % this is the usual PWL approximation
                if strcmp(rbf_type,'Piecewise linear')

                    if iter_freq_m == 1

                        freq_m_plus_1 = freq(iter_freq_m+1);
                        out_A_re_temp(iter_freq_n, iter_freq_m) = 0.5/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m_plus_1)/(1/freq_m));

                    elseif iter_freq_m == N_freq

                        freq_m_minus_1 = freq(iter_freq_m-1);
                        out_A_re_temp(iter_freq_n, iter_freq_m) = 0.5/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m)/((1/freq_m_minus_1)));

                    else

                        freq_m_plus_1 = freq(iter_freq_m+1);
                        freq_m_minus_1 = freq(iter_freq_m-1);
                        out_A_re_temp(iter_freq_n, iter_freq_m) = 0.5/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m_plus_1)/(1/freq_m_minus_1));

                    end

                else

                    % compute all RBF terms
                    out_A_re_temp(iter_freq_n, iter_freq_m) = g_i(freq_n, freq_m, epsilon, rbf_type);

                end

            end

        end

    end

    % the first and second row are reserved for L and R respectively
    out_A_re(:, 3:end) = out_A_re_temp;

end

function out_A_im = assemble_A_im(freq, epsilon, rbf_type)

%   this function assembles the A_im matrix
%   we will assume that the tau vector is identical to the freq vector

%   first get number of frequencies
    N_freq = numel(freq);

%   define the A_im output matrix
    out_A_im_temp = zeros(N_freq);
    out_A_im = zeros(N_freq, N_freq+2); 

%   we compute if the frequencies are sufficiently log spaced
    std_diff_freq = std(diff(log(1./freq)));
    mean_diff_freq = mean(diff(log(1./freq)));

%   if they are, we apply the toeplitz trick
    toeplitz_trick = std_diff_freq/mean_diff_freq<0.01;

%   if terms are evenly distributed & do not use PWL
%   then we do compute only N terms
%   else we compute all terms with brute force
    if toeplitz_trick && ~strcmp(rbf_type,'Piecewise linear')

        % define vectors R and C
        R = zeros(1,N_freq);
        C = zeros(N_freq,1);

        % for clarity the C and R computations are separated
        for iter_freq_n = 1: N_freq

            freq_n = freq(iter_freq_n);
            freq_m = freq(1);
            C(iter_freq_n, 1) = -g_ii(freq_n, freq_m, epsilon, rbf_type);

        end  

        for iter_freq_m = 1: N_freq

            freq_n = freq(1);
            freq_m = freq(iter_freq_m);
            R(1, iter_freq_m) = -g_ii(freq_n, freq_m, epsilon, rbf_type);

        end

        out_A_im_temp = toeplitz(C,R);

    else 
    % compute using brute force

        for iter_freq_n = 1: N_freq

            for iter_freq_m = 1: N_freq

                freq_n = freq(iter_freq_n);
                freq_m = freq(iter_freq_m);

                % this is the usual PWL approximation
                if strcmp(rbf_type,'Piecewise linear')

                    if iter_freq_m == 1

                        freq_m_plus_1 = freq(iter_freq_m+1);
                        out_A_im_temp(iter_freq_n, iter_freq_m) = -0.5*(2*pi*freq_n/freq_m)/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m_plus_1)/(1/freq_m));

                    elseif iter_freq_m == N_freq

                        freq_m_minus_1 = freq(iter_freq_m-1);
                        out_A_im_temp(iter_freq_n, iter_freq_m) = -0.5*(2*pi*freq_n/freq_m)/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m)/((1/freq_m_minus_1)));

                    else

                        freq_m_plus_1 = freq(iter_freq_m+1);
                        freq_m_minus_1 = freq(iter_freq_m-1);
                        out_A_im_temp(iter_freq_n, iter_freq_m) = -0.5*(2*pi*freq_n/freq_m)/(1+((2*pi*freq_n/freq_m))^2)*log((1/freq_m_plus_1)/(1/freq_m_minus_1));

                    end

                else

                   % compute all RBF terms
                   out_A_im_temp(iter_freq_n, iter_freq_m) = -g_ii(freq_n, freq_m, epsilon, rbf_type);

                end

            end
            
        end
    
    end

%   the first and second row are reserved for L and R respectively
    out_A_im(:, 3:end) = out_A_im_temp;

end

function out_M = assemble_M_1(freq, epsilon, rbf_type)

%   this function assembles the M matrix which contains the 
%   the inner products of 1st-derivative of the discretization rbfs
%   size of M matrix depends on the number of collocation points, i.e. tau vector
%   we assume that the tau vector is the inverse of the freq vector

%   first get number of frequencies
    N_freq = numel(freq);
    
%   define the M output matrix
    out_M = zeros(N_freq+2, N_freq+2); 
    out_M_temp = zeros(N_freq);
    
%   compute if the frequencies are sufficienly log spaced
    std_diff_freq = std(diff(log(1./freq)));
    mean_diff_freq = mean(diff(log(1./freq)));

%   if they are, then we apply the toeplitz trick
    toeplitz_trick = std_diff_freq/mean_diff_freq<0.01;

%   if terms are evenly distributed & do not use PWL,
%   then we do compute only N terms
%   else we compute all terms with brute force
    if toeplitz_trick && ~strcmp(rbf_type,'Piecewise linear')
        
        % define vectors R and C
        R = zeros(1,N_freq);
        C = zeros(N_freq,1);
        
        % for clarity the C and R computations are separated
        for iter_freq_n = 1: N_freq
            
            freq_n = freq(iter_freq_n);
            freq_m = freq(1);
            C(iter_freq_n, 1) = inner_prod_rbf_1(freq_n, freq_m, epsilon, rbf_type);
                    
        end  

        for iter_freq_m = 1: N_freq

            freq_n = freq(1);
            freq_m = freq(iter_freq_m);
            R(1, iter_freq_m) = inner_prod_rbf_1(freq_n, freq_m, epsilon, rbf_type);

        end

            out_M_temp = toeplitz(C,R);

    % if piecewise linear discretization
    elseif strcmp(rbf_type,'Piecewise linear') 

            out_L_temp = zeros(N_freq-1, N_freq);

                for iter_freq_n = 1:N_freq-1

                    delta_loc = log((1/freq(iter_freq_n+1))/(1/freq(iter_freq_n)));
                    out_L_temp(iter_freq_n,iter_freq_n) = -1/delta_loc;
                    out_L_temp(iter_freq_n,iter_freq_n+1) = 1/delta_loc;

                end

            out_M_temp = out_L_temp'*out_L_temp;
            
    % compute rbf with brute force
    else 
        for iter_freq_n = 1: N_freq

            for iter_freq_m = 1: N_freq

                freq_n = freq(iter_freq_n);
                freq_m = freq(iter_freq_m);
                out_M_temp(iter_freq_n, iter_freq_m) = inner_prod_rbf_1(freq_n, freq_m, epsilon, rbf_type);

            end
        end
    %             out_L_temp = chol(out_M_temp);
    end

    out_M(3:end, 3:end) = out_M_temp;

end

function out_M = assemble_M_2(freq, epsilon, rbf_type)

%   this function assembles the M matrix which contains the 
%   the inner products of 2nd-derivative of the discretization rbfs
%   size of M matrix depends on the number of collocation points, i.e. tau vector
%   we assume that the tau vector is the inverse of the freq vector

%   first get number of frequencies
    N_freq = numel(freq);
    
%   define the M output matrix    
    out_M = zeros(N_freq+2, N_freq+2);
    out_M_temp = zeros(N_freq);    
    
%   we compute if the frequencies are sufficienly log spaced    
    std_diff_freq = std(diff(log(1./freq)));
    mean_diff_freq = mean(diff(log(1./freq)));

%   if they are we apply the toeplitz trick
    toeplitz_trick = std_diff_freq/mean_diff_freq<0.01;
    
%   if terms are evenly distributed & do not use PWL,
%   then we do compute only N terms
%   else we compute all terms with brute force
    if toeplitz_trick && ~strcmp(rbf_type,'Piecewise linear')
        
        % define vectors R and C
        R = zeros(1,N_freq);
        C = zeros(N_freq,1);
        
        % for clarity the C and R computations are separated
        for iter_freq_n = 1: N_freq
                
            freq_n = freq(iter_freq_n);
            freq_m = freq(1);
            C(iter_freq_n, 1) = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type);
            
        end  

        for iter_freq_m = 1: N_freq

            freq_n = freq(1);
            freq_m = freq(iter_freq_m);
            R(1, iter_freq_m) = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type);

        end

        out_M_temp = toeplitz(C,R);

    % if piecewise linear discretization
    elseif strcmp(rbf_type,'Piecewise linear') 

        out_L_temp = zeros((N_freq-2), N_freq);

            for iter_freq_n = 1:(N_freq-2)

                delta_loc = log((1/freq(iter_freq_n+1))/(1/freq(iter_freq_n)));

                 if iter_freq_n == 1 || iter_freq_n == N_freq-2
                     
                     out_L_temp(iter_freq_n,iter_freq_n) = 2./(delta_loc^2);
                     out_L_temp(iter_freq_n,iter_freq_n+1) = -4./(delta_loc^2);
                     out_L_temp(iter_freq_n,iter_freq_n+2) = 2./(delta_loc^2);

                 else
                     
                     out_L_temp(iter_freq_n,iter_freq_n) = 1./(delta_loc^2);
                     out_L_temp(iter_freq_n,iter_freq_n+1) = -2./(delta_loc^2);
                     out_L_temp(iter_freq_n,iter_freq_n+2) = 1./(delta_loc^2);

                 end

            end

        out_M_temp = out_L_temp'*out_L_temp;
        
    % compute rbf with brute force
    else 

        for iter_freq_n = 1: N_freq
            
            for iter_freq_m = 1: N_freq

                freq_n = freq(iter_freq_n);
                freq_m = freq(iter_freq_m);
                out_M_temp(iter_freq_n, iter_freq_m) = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type);

            end
        end
    %             out_L_temp = chol(out_M_temp);
    end    

    out_M(3:end, 3:end) = out_M_temp;

end

function out_val = compute_epsilon(freq, coeff, rbf_type, shape_control)
    
%   this function is used to compute epsilon, i.e., the shape factor of
%   the rbf used for discretization. user can directly set the shape factor
%   by selecting 'Shape Factor' for the shape_control. alternatively, 
%   when 'FWHM Coefficient' is selected, the shape factor is such that 
%   the full width half maximum (FWHM) of the rbf equals to the average 
%   relaxation time spacing in log space over coeff, i.e., FWHM = delta(ln tau)/coeff

    switch rbf_type
        case 'Gaussian'
            rbf_gaussian_4_FWHM = @(x) exp(-(x).^2)-1/2;
            FWHM_coeff = 2*fzero(@(x) rbf_gaussian_4_FWHM(x), 1);
        case 'C2 Matern'
            rbf_C2_matern_4_FWHM = @(x) exp(-abs(x)).*(1+abs(x))-1/2;
            FWHM_coeff = 2*fzero(@(x) rbf_C2_matern_4_FWHM(x), 1);
        case 'C4 Matern'
            rbf_C4_matern_4_FWHM = @(x) 1/3*exp(-abs(x)).*(3+3*abs(x)+abs(x).^2)-1/2;
            FWHM_coeff = 2*fzero(@(x) rbf_C4_matern_4_FWHM(x), 1);
        case 'C6 Matern'
            rbf_C6_matern_4_FWHM = @(x) 1/15*exp(-abs(x)).*(15+15*abs(x)+6*abs(x).^2+abs(x).^3)-1/2;
            FWHM_coeff = 2*fzero(@(x) rbf_C6_matern_4_FWHM(x), 1);
        case 'Inverse quadratic'
            rbf_inverse_quadratic_4_FWHM = @(x)  1./(1+(x).^2)-1/2;
            FWHM_coeff =  2*fzero(@(x) rbf_inverse_quadratic_4_FWHM(x), 1);
        case 'Inverse quadric'
            rbf_inverse_quadric_4_FWHM = @(x)  1./sqrt(1+(x).^2)-1/2;
            FWHM_coeff = 2*fzero(@(x) rbf_inverse_quadric_4_FWHM(x), 1);
        case 'Cauchy'
            rbf_cauchy_4_FWHM = @(x)  1./(1+abs(x))-1/2;
            FWHM_coeff = 2*fzero(@(x) rbf_cauchy_4_FWHM(x) ,1);
        case 'Piecewise linear'
            FWHM_coeff = 0 ;
    end
            delta = mean(diff(log(1./freq)));
            out_val  = coeff*FWHM_coeff/delta;

end

function out_val = g_i(freq_n, freq_m, epsilon, rbf_type)

%   this function generate the elements of A_re

    alpha = 2*pi*freq_n/freq_m;

%   choose among positive definite RBFs
%   choose a function from a switch
    switch rbf_type
        case 'Gaussian'
            rbf = @(x) exp(-(epsilon*x).^2);
        case 'C0 Matern'
            rbf = @(x) exp(-abs(epsilon*x));
        case 'C2 Matern'
            rbf = @(x) exp(-abs(epsilon*x)).*(1+abs(epsilon*x));
        case 'C4 Matern'
            rbf = @(x) 1/3*exp(-abs(epsilon*x)).*(3+3*abs(epsilon*x)+abs(epsilon*x).^2);
        case 'C6 Matern'
            rbf = @(x) 1/15*exp(-abs(epsilon*x)).*(15+15*abs(epsilon*x)+6*abs(epsilon*x).^2+abs(epsilon*x).^3);
        case 'Inverse quadratic'
            rbf = @(x) 1./(1+(epsilon*x).^2);
        case 'Inverse quadric'
            rbf = @(x) 1./sqrt(1+(epsilon*x).^2);
        case 'Cauchy'
            rbf = @(x) 1./(1+abs(epsilon*x));
        otherwise
            warning('Unexpected RBF input');
    end
%   end of switch

%   integrate the following function
    integrand_g_i = @(x) 1./(1+alpha^2*exp(2*x)).*rbf(x);

%   integrate from -infty to +infty with relatively tight tolerances
    out_val = integral(integrand_g_i, -inf, inf,'RelTol',1E-9,'AbsTol',1E-9);

end

function out_val = g_ii(freq_n, freq_m, epsilon, rbf_type)

%   this function generate the elements of A_im

    alpha = 2*pi*freq_n/freq_m;

%   choose among positive definite RBFs
%   choose a function from a switch
    switch rbf_type
        case 'Gaussian'
            rbf = @(x) exp(-(epsilon*x).^2);
        case 'C0 Matern'
            rbf = @(x) exp(-abs(epsilon*x));
        case 'C2 Matern'
            rbf = @(x) exp(-abs(epsilon*x)).*(1+abs(epsilon*x));
        case 'C4 Matern'
            rbf = @(x) 1/3*exp(-abs(epsilon*x)).*(3+3*abs(epsilon*x)+abs(epsilon*x).^2);
        case 'C6 Matern'
            rbf = @(x) 1/15*exp(-abs(epsilon*x)).*(15+15*abs(epsilon*x)+6*abs(epsilon*x).^2+abs(epsilon*x).^3);
        case 'Inverse quadratic'
            rbf = @(x) 1./(1+(epsilon*x).^2);
        case 'Inverse quadric'
            rbf = @(x) 1./sqrt(1+(epsilon*x).^2);
        case 'Cauchy'
            rbf = @(x) 1./(1+abs(epsilon*x));
        otherwise
            warning('Unexpected RBF input');
    end
%   end of switch

%   integrate the following function
    integrand_g_ii = @(x) alpha./(1./exp(x)+alpha^2*exp(x)).*rbf(x);

%   integrate from -infty to +infty with relatively tight tolerances
    out_val = integral(integrand_g_ii, -inf, inf,'RelTol',1E-9,'AbsTol',1E-9);

end

function [H,c] = quad_format_combined(A_re,A_im,b_re,b_im,M,lambda) 

% this function reformats the DRT regression 
% as a quadratic program - this uses both re and im

    H = 2*((A_re'*A_re+A_im'*A_im)+lambda*M);
    H = (H'+H)/2;
    c = -2*(b_im'*A_im+b_re'*A_re);

end

function out_IP = inner_prod_rbf_1(freq_n, freq_m, epsilon, rbf_type)

%   this function output the inner product of the first derivative of the
%   rbf with respect to log(1/freq_n) and log(1/freq_m)

    a = epsilon*log(freq_n/freq_m);

%   choose among positive definite RBFs
%   choose a function from a switch
    switch rbf_type
        case 'Gaussian'
            out_IP = -epsilon*(-1+a^2)*exp(-(a^2/2))*sqrt(pi/2);
        case 'C0 Matern'
            out_IP = epsilon*(1-abs(a))*exp(-abs(a));
        case 'C2 Matern'
            out_IP = epsilon/6*(3+3*abs(a)-abs(a)^3)*exp(-abs(a));
        case 'C4 Matern'
            out_IP = epsilon/30*(105+105*abs(a)+30*abs(a)^2-5*abs(a)^3-5*abs(a)^4-abs(a)^5)*exp(-abs(a));
        case 'C6 Matern'
            out_IP = epsilon/140*(10395 +10395*abs(a)+3780*abs(a)^2+315*abs(a)^3-210*abs(a)^4-84*abs(a)^5-14*abs(a)^6-abs(a)^7)*exp(-abs(a));
        case 'Inverse quadratic'
            out_IP = 4*epsilon*(4-3*a^2)*pi/((4+a^2)^3);
        case 'Inverse quadric'
            y_n = -log(freq_n);
            y_m = -log(freq_m);
            % could only find numerical version
            rbf_n = @(y) 1./sqrt(1+(epsilon*(y-y_n)).^2);
            rbf_m = @(y) 1./sqrt(1+(epsilon*(y-y_m)).^2);
            % compute derivative
            delta = 1E-8;
            sqr_drbf_dy = @(y) 1/(2*delta).*(rbf_n(y+delta)-rbf_n(y-delta)).*1/(2*delta).*(rbf_m(y+delta)-rbf_m(y-delta));
            out_IP = integral(@(y) sqr_drbf_dy(y),-Inf,Inf);       
        case 'Cauchy'
            if a == 0
                out_IP = 2/3*epsilon;
            else
                num = abs(a)*(2+abs(a))*(4+3*abs(a)*(2+abs(a)))-2*(1+abs(a))^2*(4+abs(a)*(2+abs(a)))*log(1+abs(a));
                den = abs(a)^3*(1+abs(a))*(2+abs(a))^3;
                out_IP = 4*epsilon*num/den;
            end

        otherwise
            warning('Unexpected RBF input.');
    end
%   end switch
    
end

function [out_gamma, freq_fine] = map_array_to_gamma(freq_map, freq_coll, x, epsilon, rbf_type)

%   this function map x, i.e., the magnitude at the specific frequency, to 
%   gamma, i.e., the DRT profile. rbf refer to the discretization function. 
%   on multiplying rbf function with x, one can obtain the gamma profile.

%   choose among positive definite RBFs
%   choose a function from a switch
    switch rbf_type
        case 'Gaussian'
            rbf = @(y, y0) exp(-(epsilon*(y-y0)).^2);
        case 'C0 Matern'
            rbf = @(y, y0)  exp(-abs(epsilon*(y-y0)));
        case 'C2 Matern'
            rbf = @(y, y0)  exp(-abs(epsilon*(y-y0))).*(1+abs(epsilon*(y-y0)));
        case 'C4 Matern'
            rbf = @(y, y0)  1/3*exp(-abs(epsilon*(y-y0))).*(3+3*abs(epsilon*(y-y0))+abs(epsilon*(y-y0)).^2);
        case 'C6 Matern'
            rbf = @(y, y0)  1/15*exp(-abs(epsilon*(y-y0))).*(15+15*abs(epsilon*(y-y0))+6*abs(epsilon*(y-y0)).^2+abs(epsilon*(y-y0)).^3);
        case 'Inverse quadratic'
            rbf = @(y, y0)  1./(1+(epsilon*(y-y0)).^2);
        case 'Inverse quadric'
            rbf = @(y, y0)  1./sqrt(1+(epsilon*(y-y0)).^2);
        case 'Cauchy'
            rbf = @(y, y0)  1./(1+abs(epsilon*(y-y0))); 
        case 'Piecewise linear'
            out_gamma = x;
            freq_fine = freq_coll';
            return

        otherwise
            warning('Unexpected RBF input');
    end
    % end of switch

    freq_fine = freq_map;
    y0 = -log(freq_coll);
    out_gamma = zeros(size(freq_map))';

    for iter_freq_map = 1: numel(freq_map)

        freq_map_loc = freq_map(iter_freq_map);
        y = -log(freq_map_loc);
        out_gamma(iter_freq_map) = x'*rbf(y, y0);

    end

end

function out_IP = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)

%   this function output the inner product of the second derivative of the
%   rbf with respect to log(1/freq_n) and log(1/freq_m)

    a = epsilon*log(freq_n/freq_m);

%   choose among positive definite RBFs
%   choose a function from a switch
    switch rbf_type
        case 'Gaussian'
            out_IP = epsilon^3*(3-6*a^2+a^4)*exp(-(a^2/2))*sqrt(pi/2);
        case 'C0 Matern'
            out_IP = epsilon^3*(1+abs(a))*exp(-abs(a));
        case 'C2 Matern'
            out_IP = epsilon^3/6*(3+3*abs(a)-6*abs(a)^2+abs(a)^3)*exp(-abs(a));
        case 'C4 Matern'
            out_IP = epsilon^3/30*(45+45*abs(a)-15*abs(a)^3-5*abs(a)^4+abs(a)^5)*exp(-abs(a));
        case 'C6 Matern'
            out_IP = epsilon^3/140*(2835+2835*abs(a)+630*abs(a)^2-315*abs(a)^3-210*abs(a)^4-42*abs(a)^5+abs(a)^7)*exp(-abs(a));
        case 'Inverse quadratic'
            out_IP = 48*(16+5*a^2*(-8 + a^2))*pi*epsilon^3/((4 + a^2)^5);
        case 'Inverse quadric'
            y_n = -log(freq_n);
            y_m = -log(freq_m);
            % could only find numerical version
            rbf_n = @(y) 1./sqrt(1+(epsilon*(y-y_n)).^2);
            rbf_m = @(y) 1./sqrt(1+(epsilon*(y-y_m)).^2);
            % compute derivative
            delta = 1E-4;
            sqr_drbf_dy = @(y) 1/(delta^2).*(rbf_n(y+delta)-2*rbf_n(y)+rbf_n(y-delta)).*1/(delta^2).*(rbf_m(y+delta)-2*rbf_m(y)+rbf_m(y-delta));
            out_IP = integral(@(y) sqr_drbf_dy(y),-Inf,Inf);       
        case 'Cauchy'
            if a == 0
                out_IP = 8/5*epsilon^3;
            else
                num = abs(a)*(2+abs(a))*(-96 +abs(a)*(2+abs(a))*(-30 +abs(a)*(2+abs(a)))*(4+abs(a)*(2+abs(a))))+... 
                        12*(1+abs(a))^2*(16+abs(a)*(2+abs(a))*(12+abs(a)*(2+abs(a))))*log(1+abs(a));
                den = abs(a)^5*(1+abs(a))*(2+abs(a))^5;
                out_IP = 8*epsilon^3*num/den;
            end    
        otherwise
            warning('Unexpected RBF input.');
    end
%   end switch

end

function out_dict = HT_single_est(theta, Z_exp, A, A_H, M, N_freqs, N_taus)

%   step 1 - identify the vector of hyperparameters
    options = optimset('Display','iter','TolX', 1E-10, 'TolFun', 1E-10,'Algorithm','quasi-newton');
    opt_theta = fminunc(@(c) NMLL_fct(c, Z_exp, A, M, N_freqs, N_taus), log(theta), options);
%     options = optimset('Display','iter','TolX', 1E-10, 'TolFun', 1E-10,'Algorithm','sqp');
%     lower = [-10,-10,-10];
%     upper = [10,10,10];
%     opt_theta = fmincon(@(c) NMLL_fct(c, Z_exp, A, M, N_freqs, N_taus), log(theta),[],[],[],[],lower,upper,[],options); 
%   collect the optimized theta'
    opt_theta = exp(opt_theta);
    sigma_n = opt_theta(1);
    sigma_beta = opt_theta(2);
    sigma_lambda = opt_theta(3);

%   step 2 - compute the pdf's of data regression
%   $K_agm = A.T A +\lambda L.T L$
    W = 1/(sigma_beta^2)*eye(N_taus+1) + 1/(sigma_lambda^2)*M;
    K_agm = 1/(sigma_n^2)*(A'*A) + W;

%   Cholesky factorization
    L_agm = chol(K_agm,'lower'); % the default matrix for matlab and python are npt the same
    inv_L_agm = inv(L_agm);
    inv_K_agm = inv_L_agm'*inv_L_agm;

%   compute the gamma ~ N(mu_gamma, Sigma_gamma)
    Sigma_gamma = inv_K_agm;
    mu_gamma = 1/(sigma_n^2)*(Sigma_gamma*A')*Z_exp;

%   compute Z ~ N(mu_Z, Sigma_Z) from gamma
    mu_Z = A*mu_gamma;
    Sigma_Z = A*(Sigma_gamma*A') + sigma_n^2*eye(N_freqs);
    
%   compute Z_DRT ~ N(mu_Z_DRT, Sigma_Z_DRT) from gamma
    A_DRT = A(:,2:end);
    mu_gamma_DRT = mu_gamma(2:end);
    Sigma_gamma_DRT = Sigma_gamma(2:end,2:end);
    mu_Z_DRT = A_DRT*mu_gamma_DRT;
    Sigma_Z_DRT = A_DRT*(Sigma_gamma_DRT*A_DRT');
    
%   compute Z_H_conj ~ N(mu_Z_H_conj, Sigma_Z_H_conj) from gamma   
    mu_Z_H = A_H*mu_gamma(2:end);
    Sigma_Z_H = A_H*(Sigma_gamma(2:end,2:end)*A_H');

    out_dict = struct('mu_gamma', mu_gamma,...
                      'Sigma_gamma', Sigma_gamma,...
                      'mu_Z', mu_Z,...
                      'Sigma_Z', Sigma_Z,...
                      'mu_Z_DRT', mu_Z_DRT,...
                      'Sigma_Z_DRT', Sigma_Z_DRT,...
                      'mu_Z_H', mu_Z_H,...
                      'Sigma_Z_H', Sigma_Z_H,...
                      'theta', opt_theta);

end

function [ZIKKRe, ZIIKKIm]=KK(Z, freq, coeff, shape_control, rbf_type)  
    
    omega_vec = 2*pi*freq;
    N_freqs = numel(freq);
    N_taus = numel(freq);
    
%   Step 1 construct the A matrix, A_H_matrix
    epsilon = compute_epsilon(freq, coeff, rbf_type, shape_control);

    A_re_0 = assemble_A_re(freq, epsilon, rbf_type);
    A_im_0 = assemble_A_im(freq, epsilon, rbf_type);

    A_H_re = A_re_0(:,3:end);
    A_H_im = A_im_0(:,3:end);

%   add resistence column and inductance column to A_re and A_im, remove
%   the unused zero column from A_re_0 and A_im_0
    A_re = [ones(N_freqs,1),A_re_0(:,3:end)];
    A_im = [omega_vec,A_im_0(:,3:end)];
    
    b_re = real(Z);
    b_im = imag(Z);

%   Step 2 construct L_matrix
    M = assemble_M_2(freq, epsilon, rbf_type);
            

    M = M(2:end,2:end);
    
%   Step 3 testing HT_single_est 
%   try until no error occur for the HT_single_est
    while true
        try
            % Randomly select three inital points between 10^4 to 10^-4 for optimization
            theta_0 = 10.^(8*rand(3,1)-4);
            out_dict_real = HT_single_est(theta_0, b_re, A_re, A_H_im, M, N_freqs, N_taus); % for v2 it input the  
            out_dict_imag = HT_single_est(out_dict_real.theta, b_im, A_im, A_H_re, M, N_freqs, N_taus);
            break
            
        catch
            disp('Error Occur, Try Another Inital Condition')
            
        end
    end
    
%   Step 4 testing EIS scoring
    N_MC_samples = 50000;

    handles.out_scores = EIS_score(theta_0, freq, Z, out_dict_real, out_dict_imag, N_MC_samples);

%   Step 5 print out scores
    fprintf('The EIS scores are as follow:\n');
    fprintf('s_res_re = %f %f %f\n', out_scores.s_res_re);
    fprintf('s_res_im = %f %f %f\n', out_scores.s_res_im);
    fprintf('s_mu_re = %f \n', out_scores.s_mu_re);
    fprintf('s_mu_im = %f \n', out_scores.s_mu_im);
    fprintf('s_HD_re = %f \n', out_scores.s_HD_re);
    fprintf('s_HD_im = %f \n', out_scores.s_HD_im);
    fprintf('s_JSD_re = %f \n', out_scores.s_JSD_re);
    fprintf('s_JSD_im = %f \n', out_scores.s_JSD_im);
    fprintf('opt_theta_real = %f %f %f\n', out_dict_real.theta);
    fprintf('opt_theta_imag = %f %f %f\n', out_dict_imag.theta);

%   Step 6 shows the band and the hilbert fitting in the real part and the imag part
%   1. real data
%   1.1 Bayesian regression
    handles.mu_Z_re = out_dict_real.mu_Z;
    ZIKKRe=handles.mu_Z_re(:);
    disp(ZIKKRe);

%   2. imaginary data
%   2.1 Bayesian regression
    handles.mu_Z_im = out_dict_imag.mu_Z;
    ZIIKKIm=handles.mu_Z_im(:);
    disp(ZIIKKIm);
end