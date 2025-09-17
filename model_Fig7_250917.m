% setting parameters
dt = 0.01;  
T_total = 1000;  
t = 0:dt:T_total;

% Intrinsic angular velocities
omega_A = 2 * pi / 22.7; % A: AVP  % change here (22.7 -> 23.7) for Sup. Fig. 11d
omega_B = 2 * pi / 24.3; % B: non-AVP (nA)

% Parameter search ranges
K_BtoA_values = linspace(0, 0.2, 41);          % = K_nA_to_A
K_VIP_values = linspace(0, 0.4, 41);           % = K_V
K_AtoB_ratio_list = [0, 0.2, 1, 5, 25];        % = K_A_to_nA

figure('Position', [100, 100, 1200, 900]);

for r = 1:length(K_AtoB_ratio_list)
    ratio = K_AtoB_ratio_list(r);

    per_A_map = NaN(length(K_BtoA_values), length(K_VIP_values));
    per_B_map = NaN(length(K_BtoA_values), length(K_VIP_values));
    sync_map   = NaN(length(K_BtoA_values), length(K_VIP_values));

    for i = 1:length(K_BtoA_values)
        K_BtoA = K_BtoA_values(i);
        K_AtoB = ratio * K_BtoA;

        for j = 1:length(K_VIP_values)
            K_VIP = K_VIP_values(j);

            theta_A = zeros(1, length(t));
            theta_B = zeros(1, length(t));
            theta_A(1) = 0;
            theta_B(1) = 0.1;

            for k = 1:length(t)-1
                dtheta_A = omega_A + K_BtoA * sin(theta_B(k) - theta_A(k)) - K_VIP*max(sin(theta_B(k)), 0)*max(sin(theta_A(k)- 3*pi/2), 0); 
                % change above term (1*pi/2) for Sup.Fig.11a (0*pi/2), Sup.Fig.11b (2*pi/2), and Sup.Fig.11c (3*pi/2) 
                dtheta_B = omega_B + K_AtoB * sin(theta_A(k) - theta_B(k));
                theta_A(k+1) = theta_A(k) + dtheta_A * dt;
                theta_B(k+1) = theta_B(k) + dtheta_B * dt;
            end

            n_stable = floor(length(t)/2);
            delta_t = t(end) - t(n_stable);
            delta_theta_A = theta_A(end) - theta_A(n_stable);
            delta_theta_B = theta_B(end) - theta_B(n_stable);

            per_A = 2*pi * delta_t / delta_theta_A;
            per_B = 2*pi * delta_t / delta_theta_B;

            per_A_map(i,j) = per_A;
            per_B_map(i,j) = per_B;
            per_A_map(isinf(per_A_map)) = NaN;
            per_B_map(isinf(per_B_map)) = NaN;

            if abs(per_A - per_B) < 0.2
                sync_map(i,j) = per_A;
            elseif isfinite(per_A) && isfinite(per_B)
                sync_map(i,j) = -1; % A period and B period is different -> black
            else
                sync_map(i,j) = NaN; % period is not determined -> gray
            end
        end
    end

    % hetmap presentation

    % A period
    subplot(3, 5, r);
    h1 = imagesc(K_VIP_values, K_BtoA_values, per_A_map, 'AlphaData', ~isnan(per_A_map));
    set(gca, 'YDir', 'normal');
    colormap(gca, 'jet');
    clim([22, 26]);
    colorbar;
    xlabel('K_{VIP}');
    ylabel('K_{nA→A}');
    title(['AVP period (K_{A→nA} = ', num2str(ratio), '×K_{nA→A})']);
    set(gca,'Color', [0.5 0.5 0.5]); % NaN -> gray

    % B (nA) period
    subplot(3, 5, r+5);
    h2 = imagesc(K_VIP_values, K_BtoA_values, per_B_map, 'AlphaData', ~isnan(per_B_map));
    set(gca, 'YDir', 'normal');
    colormap(gca, 'jet');
    clim([22, 26]);
    colorbar;
    xlabel('K_{VIP}');
    ylabel('K_{nA→A}');
    title(['non-AVP period (K_{A→nA} = ', num2str(ratio), '×K_{nA→A})']);
    set(gca,'Color', [0.5 0.5 0.5]); % NaN -> gray

    % Ensemble period
    subplot(3, 5, r+10);
    imagesc(K_VIP_values, K_BtoA_values, sync_map, 'AlphaData', ~isnan(sync_map));
    set(gca, 'YDir', 'normal');
    colormap(gca, [0 0 0; jet(256)]); % black (not synced), jet(synced)
    clim([22, 26]);
    colorbar;
    xlabel('K_{VIP}');
    ylabel('K_{nA→A}');
    title(['ensemble period (K_{A→nA} = ', num2str(ratio), '×K_{nA→A})']);
    set(gca,'Color', [0.5 0.5 0.5]); % NaN -> gray
end
