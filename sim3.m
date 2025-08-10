%%General Gate simulation 
%Lena Aly - NYU Abu Dhabi 
%lwa2019@nyu.edu

clear; clc; close all;

global dt start rho_t rho amplitude_mw Lz Lx Ly Sz Sx Sy Fx Fy Fz



% Parameters
omega0_L = 11.0 * 2 * pi;
omega0_S = 11.0 * 2 *pi;
omegaL_mw = 11.0 * 2 * pi;
omegaS_mw = 15.0 * 2 * pi;


deltaL_mw = omegaL_mw - omega0_L;
deltaS_mw = omegaS_mw -omega0_S;
amplitude_mw = 2 * pi;
phase_0 = 90;  % not used directly
J = 20 * pi; %coupling term


% Pauli Matrices
Iz = [1 0; 0 -1] / 2;
Ix = [0 1; 1 0] / 2;
Iy = [0 -1i; 1i 0] / 2;

%Two-spin Pauli Matrices
Lz=kron(Iz, eye(2));
Sz=kron(eye(2),Iz);
Lx=kron(Ix, eye(2));
Sx=kron(eye(2),Ix);
Ly=kron(Iy, eye(2));
Sy=kron(eye(2),Iy);
Fx=Lx+Sx;
Fy=Ly+Sy;
Fz=Lz+Sz;

% Time parameters
dt = (pi/(1*J))*0.01;


start =0;
CNOT_gate = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];

% Initial pure state density matrix 
psi1 = [0; 1];
psi0 = [1; 0];
r =  [0;0;1;0];
rho = r*r';
rho_t = zeros(4, 4, 0);  % ensemble tracking
disp("After step definition:");
disp(rho);

%__________________________________________________

function pulse(theta, axis, direction)

global dt start rho_t rho Lx Ly Sx Sy Fx Fy amplitude_mw;

if axis == "12x"
    H = amplitude_mw * Fx; 

elseif axis == "12y"
    H = amplitude_mw * Fy;

elseif axis == "1x"
    H = amplitude_mw * Lx;

elseif axis == "2x"
    H = amplitude_mw * Sx;

elseif axis == "1y"
    H = amplitude_mw * Ly;

elseif axis == "2y"
    H = amplitude_mw * Sy;

end 

NoP_segment = round( theta *(pi/180)/(dt* amplitude_mw));
U = expm(-1i * H * dt);

    for i= start +1 : start +NoP_segment
        rho_t(:, :, i) = rho;
        if direction == 1
           rho = U*rho*U'; 

        elseif direction == -1
           rho = U'*rho*U; 
        end

    end 
start = start + NoP_segment;
rho_t(:,:,start) = rho;
end 
%__________________________________________________________________________

function evolution(H, value)

global dt rho_t rho start;

NoP_evolution = round ((value)/dt);

U = expm(-1i * H * dt);
NoP_segment = NoP_evolution;

 for i= start +1 : start +NoP_segment
        rho_t(:, :, i) = rho;
           rho = U*rho*U'; 
 end 

 start = start + NoP_segment;
 rho_t(:,:,start) = rho;
end

%________________________________________________________________

function L(J)

global Lz Sz Lx Sx Ly Sy amplitude_mw Fz

H_evol = amplitude_mw*Fz+ J*(Lz*Sz);
value = pi/(2*J);

evolution( H_evol,value);
pulse(180, "12x", 1);
evolution( H_evol, value);
pulse(180, "12x", -1);

end 

% Pulse Sequence 

%SWAP
%{
pulse(180, "12x", -1);
L(J)
pulse(180, "12x", 1);
pulse(180, "12y",-1);
L (J)
pulse(180, "12y", 1);
L(J)
%}

%_______________________________________
%CNOT



value = pi/(2*J);
H_evol = amplitude_mw*Fz+ J*(Lz*Sz);

pulse( 90, "2y", -1);
L(J);
pulse(90, "2y", 1);
pulse( 90, "1x", -1);
pulse( 90, "1x", 1);
pulse( 90, "1y", -1);
pulse(90, "1x", -1);

%}

%pulse( 90, "2y", 1);
%L(J);
%pulse( 90, "2x", 1);
 


%}


%Spin-Spin Interaction Hamiltonain
pulse( 90, "2y", 1);
evolution(H_evol, 4*J);
pulse( 90, "1x", -1);

pulse(180, "1y", 1);

pulse( 90, "1x", 1);
evolution(H_evol, 4*J);
pulse( 90, "1x", -1);
pulse( 180, "1y", -1);
pulse( 90, "1x", 1);
pulse( 180, "12x", 1);
evolution(H_evol, 4*J);
pulse( 90, "1x", -1);
pulse( 180, "1y", 1);
pulse( 90, "1x", 1);
evolution(H_evol, 4*J);
pulse( 90, "1x", -1);
pulse( 180, "1y", -1);
pulse( 90, "1x", 1);
pulse( 180, "12x", -1);
pulse( 90, "2y", -1);
pulse(90, "1x", -1);
pulse( 90, "1y", 1);
pulse( 180, "12x", 1);

%________________________________________

% X-gate

% pulse(180, "1x", 1)

total_steps = start; 

% Extract probabilities
P0 = zeros(1, total_steps);
P1 = zeros(1, total_steps);
P2 = zeros(1, total_steps);
P3= zeros(1, total_steps);

 
for i = 1:total_steps

    P0(i) = real(rho_t(1,1,i)); 
    P1(i) = real(rho_t(4,4,i));
    P2(i) = real(rho_t(2,2,i));
    P3(i) = real(rho_t(3,3,i));

end

t = (0 : total_steps - 1) * dt;



 
%trace(rho^2)
%time = (0:total_steps-1) * dt;

%% Plot probabilities
figure(1);

plot(t, P0, 'r');hold on;plot(t, P2, 'b');...
hold on;plot(t, P3, 'g');hold on;plot(t, P1, 'black');
title('Population evolution under CNOT');
xlabel('Time (a.u.)');
ylabel('Probability');
legend('00', '01', '10', '11');
grid on;
 
% Plot phase evolution of amplitudes
rho_offdiag = squeeze(rho_t(1,4,:) + rho_t(1,3,:)+ rho_t(1,2)+rho_t(2,3)+rho_t(2,4)+rho_t(3,4));  % ρ₀₁ = α*β*
figure(2);
plot(t, real(rho_offdiag), 'r'); hold on;
plot(t, imag(rho_offdiag), 'b');hold on; plot(t, P0, '');
title('Phase evolution from Re/Im[ρ₀₁]');
xlabel('Time (a.u.)');
ylabel('Amplitude');
legend('Re(ρ_{01})', 'Im(ρ_{01})');
grid on;


disp(total_steps)

%% ================================================================
%%  Bloch‑sphere visualisation for the final state  (ADDED LINES)
%% ================================================================
% Reduced density matrices of spin‑1 and spin‑2
rho_final = rho;                                          % rho after the sequence
rho1 = [rho_final(1,1)+rho_final(2,2),  rho_final(1,3)+rho_final(2,4); ...
        rho_final(3,1)+rho_final(4,2),  rho_final(3,3)+rho_final(4,4)];
rho2 = [rho_final(1,1)+rho_final(3,3),  rho_final(1,2)+rho_final(3,4); ...
        rho_final(2,1)+rho_final(4,3),  rho_final(2,2)+rho_final(4,4)];

% Bloch‑vector components  ⟨σx⟩,⟨σy⟩,⟨σz⟩  (Ix = σx/2 etc.)
bloch1 = [2*real(trace(rho1*Ix)),  2*real(trace(rho1*Iy)),  2*real(trace(rho1*Iz))];
bloch2 = [2*real(trace(rho2*Ix)),  2*real(trace(rho2*Iy)),  2*real(trace(rho2*Iz))];
%{
% ----------------------------------------------------------------
%   Plot two Bloch spheres side‑by‑side
% ----------------------------------------------------------------
figure(3); clf;

% ---------- Spin 1 ----------
subplot(1,2,1);
[Xs,Ys,Zs] = sphere(60);
surf(Xs,Ys,Zs, 'FaceAlpha',0.08, 'EdgeColor','none');  colormap gray;  hold on;
plot3([0 bloch1(1)], [0 bloch1(2)], [0 bloch1(3)], 'r', 'LineWidth',2);
axis equal vis3d; grid on; xlabel('X'); ylabel('Y'); zlabel('Z');
title('Bloch sphere – spin 1'); view([135 30]);

% ---------- Spin 2 ----------
subplot(1,2,2);
surf(Xs,Ys,Zs, 'FaceAlpha',0.08, 'EdgeColor','none');  colormap gray;  hold on;
plot3([0 bloch2(1)], [0 bloch2(2)], [0 bloch2(3)], 'b', 'LineWidth',2);
axis equal vis3d; grid on; xlabel('X'); ylabel('Y'); zlabel('Z');
title('Bloch sphere – spin 2'); view([135 30]);
%}