%% half_car_master.m
% Modular half-car master script
% Sections: eigenvalue analysis, base sim, speed sweep, q1..q4 plotting
clc; clear; close all;

%% ====================== 1. FLAGS (toggle sections) ======================
run_eigen       = true;   % eig analysis of linearised 4-DOF model
run_baseSim     = true;   % single simulation at base speed (50 km/h)
run_speedSweep  = true;   % multi-speed comparison (overlay)
run_qVariables  = true;   % q1..q4 formatted plots

%% ====================== 2. COMMON PARAMETERS ============================
% Mass / inertia
M   = 600;        % sprung mass (half-car) [kg]
J   = 300;        % pitch inertia [kg*m^2]
mf  = 40;         % front unsprung mass [kg]
mr  = 40;         % rear unsprung mass [kg]

% Stiffness & damping (tuned / as requested)
kf  = 2.0e4;      % front suspension stiffness [N/m]
kr  = 2.0e4;      % rear suspension stiffness [N/m]
cf  = 1000;       % front damping [N*s/m]
cr  = 1000;       % rear damping [N*s/m]

ktf = 3.0e5;      % front tyre stiffness [N/m]
ktr = 3.0e5;      % rear tyre stiffness [N/m]

% Geometry & gravity
lf  = 1.3;        % CoM to front axle [m]
lr  = 1.3;        % CoM to rear axle [m]
g   = 9.81;       % gravity [m/s^2]

% Road bump (user confirmed)
D0  = 0;          % start of bump [m]
W   = 20;         % bump width [m]
H   = 1;          % bump height [m]  <-- user choice

% Time discretization
dt    = 0.001;    % time step [s]
t_end = 5;        % simulation duration [s]
t     = 0:dt:t_end;
N     = length(t);

% Base speed (km/h -> m/s)
v_base_kmph = 50;
v_base = v_base_kmph / 3.6;

% Speeds for sweep (km/h)
speeds_kmph = [50, 100, 120, 160];

%% ====================== 3. EIGENVALUE ANALYSIS ==========================
if run_eigen
    fprintf('\n--- Eigenvalue analysis (4-DOF linearised) ---\n');
    % Map parameters into the simplified linear stiffness mapping used earlier
    K = (kf + kr)/2;        % representative suspension stiffness (assume symmetric)
    Kt = (ktf + ktr)/2;     % representative tyre stiffness
    L_half = lf;            % mapping used previously (lf=lr=1.3)
    
    % Mass matrix (4 DOF: [zb, theta, zaxf, zaxr] mapping to [M, J, m, m])
    M_mat = [M     0     0   0;
             0     J     0   0;
             0     0     mf  0;
             0     0     0   mr];
    
    % Stiffness mapping consistent with earlier snippet
    % (the mapping uses K and geometric coupling via L_half)
    K_mat = [2*K       2*K*L_half   -K        -K;
             2*K*L_half 2*K*L_half^2 -K*L_half -K*L_half;
             -K       -K*L_half    K + Kt    0;
             -K       -K*L_half    0         K + Kt];
    
    % Solve generalised eigenvalue problem
    [phi, lambda] = eig(K_mat, M_mat);
    lambda_vec = diag(lambda);
    
    % Remove (or keep) tiny/negative numerical noise
    lambda_vec(lambda_vec < 0) = 0;
    
    omega = sqrt(lambda_vec);       % rad/s
    freq_hz = omega / (2*pi);       % Hz
    
    disp('Natural Frequencies (Hz):');
    disp(freq_hz);
    disp('Mode shapes (columns correspond to eigenvectors):');
    disp(phi);
    
    % Quick plot of mode participation (scaled)
    figure('Name','Eigen Modes','NumberTitle','off');
    for k = 1:size(phi,2)
        subplot(2,2,k);
        bar(abs(phi(:,k)));
        xticklabels({'zb','\theta','z_{axf}','z_{axr}'});
        ylabel('|mode amplitude|');
        title(sprintf('Mode %d: f = %.3f Hz', k, freq_hz(k)));
        grid on;
    end
end

%% ====================== 4. BASE SIMULATION (single speed) ==============
if run_baseSim
    fprintf('\n--- Base simulation at %d km/h ---\n', v_base_kmph);
    
    % Initialize state vector x = [zb; zb_dot; theta; theta_dot; zf; zf_dot; zr; zr_dot]
    x = zeros(8, N);
    
    % Road inputs for base speed
    x_vehicle = v_base * t;
    zrf = roadInput(x_vehicle, D0, W, H);
    zrr = roadInput(x_vehicle - (lf + lr), D0, W, H);
    
    % Time stepping (forward Euler)
    for i = 1:N-1
        dx = halfcarEoM(x(:,i), M, J, mf, mr, kf, kr, cf, cr, ktf, ktr, lf, lr, zrf(i), zrr(i));
        x(:,i+1) = x(:,i) + dx * dt;
    end
    
    % Extract signals for plotting
    zb   = x(1,:);    theta = x(3,:);
    zaxf = x(5,:);    zaxr = x(7,:);
    rel_front = zaxf - zrf;
    rel_rear  = zaxr - zrr;
    
    % Plot results (compact)
    figure('Name','Base Simulation (single speed)','NumberTitle','off');
    subplot(4,1,1);
    plot(t, zrf, 'r--', t, zrr, 'b--'); hold on;
    ylabel('Road elevation [m]');
    legend('Front road input','Rear road input');
    title(sprintf('Road Inputs at %d km/h (H=%.2f m)', v_base_kmph, H));
    grid on;
    
    subplot(4,1,2);
    plot(t, zb); ylabel('Body heave z_b [m]'); grid on;
    
    subplot(4,1,3);
    plot(t, theta*180/pi); ylabel('Pitch \theta [deg]'); grid on;
    
    subplot(4,1,4);
    plot(t, rel_front, 'r', t, rel_rear, 'b'); ylabel('Suspension deflection [m]');
    legend('Front deflection','Rear deflection');
    xlabel('Time [s]'); grid on;
end

%% ====================== 5. SPEED SWEEP (overlay comparisons) ============
if run_speedSweep
    fprintf('\n--- Speed sweep: %s km/h ---\n', num2str(speeds_kmph));
    % Prepare figure
    figure('Name','Speed sweep (overlay)','NumberTitle','off','Units','normalized','Position',[0.1 0.1 0.7 0.7]);
    
    % We'll overlay body heave and front deflection as examples
    subplot(2,1,1); hold on; grid on; title('Body heave: speed comparison'); ylabel('z_b [m]');
    subplot(2,1,2); hold on; grid on; title('Front suspension deflection: speed comparison'); ylabel('Front deflection [m]');
    
    legends_heave = cell(1,length(speeds_kmph));
    legends_front = cell(1,length(speeds_kmph));
    for idx = 1:length(speeds_kmph)
        vk = speeds_kmph(idx);
        v_ms = vk/3.6;
        
        % simulate
        x = zeros(8,N);
        x_vehicle = v_ms * t;
        zrf = roadInput(x_vehicle, D0, W, H);
        zrr = roadInput(x_vehicle - (lf + lr), D0, W, H);
        for i = 1:N-1
            dx = halfcarEoM(x(:,i), M, J, mf, mr, kf, kr, cf, cr, ktf, ktr, lf, lr, zrf(i), zrr(i));
            x(:,i+1) = x(:,i) + dx * dt;
        end
        
        zb = x(1,:);
        zaxf = x(5,:);
        rel_front = zaxf - zrf;
        
        % plot
        subplot(2,1,1);
        plot(t, zb);
        subplot(2,1,2);
        plot(t, rel_front);
        
        legends_heave{idx} = sprintf('v = %d km/h', vk);
        legends_front{idx} = sprintf('v = %d km/h', vk);
    end
    
    subplot(2,1,1);
    legend(legends_heave,'Location','best'); xlabel('Time [s]');
    subplot(2,1,2);
    legend(legends_front,'Location','best'); xlabel('Time [s]');
end

%% ====================== 6. q1..q4 formatted visualization ==============
if run_qVariables
    fprintf('\n--- q1..q4 formatted run (base speed) ---\n');
    % run at base speed
    x = zeros(8,N);
    x_vehicle = v_base * t;
    zrf = roadInput(x_vehicle, D0, W, H);
    zrr = roadInput(x_vehicle - (lf + lr), D0, W, H);
    for i = 1:N-1
        dx = halfcarEoM(x(:,i), M, J, mf, mr, kf, kr, cf, cr, ktf, ktr, lf, lr, zrf(i), zrr(i));
        x(:,i+1) = x(:,i) + dx * dt;
    end
    
    % Map q variables
    q1 = x(1,:); q1d = x(2,:);
    q2 = x(3,:); q2d = x(4,:);
    q3 = x(5,:); q3d = x(6,:);
    q4 = x(7,:); q4d = x(8,:);
    rel_front = q3 - zrf; rel_rear = q4 - zrr;
    
    % Plot grid layout 4x2
    figure('Name','q1..q4 layout','NumberTitle','off','Units','normalized','Position',[0.05 0.05 0.85 0.85]);
    sgtitle('Half-Car: q1..q4 displacements and velocities (base speed)');
    subplot(4,2,1); plot(t, q1, 'b'); ylabel('q1 (zb) [m]'); grid on;
    subplot(4,2,2); plot(t, q1d, 'b'); ylabel('q1d [m/s]'); grid on;
    subplot(4,2,3); plot(t, q2*180/pi, 'r'); ylabel('q2 (\theta) [deg]'); grid on;
    subplot(4,2,4); plot(t, q2d*180/pi, 'r'); ylabel('\theta_{dot} [deg/s]'); grid on;
    subplot(4,2,5); plot(t, q3, 'y'); ylabel('q3 (z_{axf}) [m]'); grid on;
    subplot(4,2,6); plot(t, q3d, 'y'); ylabel('q3d [m/s]'); grid on;
    subplot(4,2,7); plot(t, q4, 'm'); ylabel('q4 (z_{axr}) [m]'); grid on;
    subplot(4,2,8); plot(t, q4d, 'm'); ylabel('q4d [m/s]'); xlabel('Time [s]'); grid on;
    
    % Optional figure for road & relative defs
    figure('Name','Road and relative deflections','NumberTitle','off');
    subplot(3,1,1);
    plot(t, zrf, 'r--', t, zrr, 'b--'); ylabel('Road input [m]'); legend('Front input','Rear input'); grid on;
    subplot(3,1,2); plot(t, rel_front); ylabel('Front suspension deflection [m]'); grid on;
    subplot(3,1,3); plot(t, rel_rear); ylabel('Rear suspension deflection [m]'); xlabel('Time [s]'); grid on;
end

%% ====================== 7. HELPER FUNCTIONS ============================
function dx = halfcarEoM(x, M, J, mf, mr, kf, kr, cf, cr, ktf, ktr, lf, lr, zrf, zrr)
    % State vector: x = [zb; zb_dot; th; th_dot; zf; zf_dot; zr; zr_dot]
    zb = x(1); zb_dot = x(2);
    th = x(3); th_dot = x(4);
    zf = x(5); zf_dot = x(6);
    zr = x(7); zr_dot = x(8);
    % Relative motions
    delta_f = (zb + lf*th - zf);
    delta_r = (zb - lr*th - zr);
    delta_f_dot = (zb_dot + lf*th_dot - zf_dot);
    delta_r_dot = (zb_dot - lr*th_dot - zr_dot);
    % Forces (spring + damper)
    Fs_f = kf*delta_f + cf*delta_f_dot;
    Fs_r = kr*delta_r + cr*delta_r_dot;
    % Tyre forces (no tensile tyre force / simple linear)
    Ft_f = ktf * (zf - zrf);
    Ft_r = ktr * (zr - zrr);
    % Equations of Motion
    zb_ddot = (-Fs_f - Fs_r) / M;
    th_ddot = (-lf*Fs_f + lr*Fs_r) / J;
    zf_ddot = (Fs_f - Ft_f) / mf;
    zr_ddot = (Fs_r - Ft_r) / mr;
    dx = [zb_dot; zb_ddot; th_dot; th_ddot; zf_dot; zf_ddot; zr_dot; zr_ddot];
end

function z = STEP(d, d1, z1, d2, z2)
    % Smooth cubic step function between (d1,z1) and (d2,z2)
    if d <= d1
        z = z1;
    elseif d >= d2
        z = z2;
    else
        alpha = (d - d1) / (d2 - d1);
        z = z1 + (z2 - z1) * (alpha^2 * (3 - 2*alpha));
    end
end

function q = roadInput(d, D0, W, H)
    % Generates smooth bump: rises then falls using two STEP segments
    q = zeros(size(d));
    for ii = 1:length(d)
        q(ii) = STEP(d(ii), D0, 0, D0 + W/2, H) + STEP(d(ii), D0 + W/2, H, D0 + W, 0);
    end
end
