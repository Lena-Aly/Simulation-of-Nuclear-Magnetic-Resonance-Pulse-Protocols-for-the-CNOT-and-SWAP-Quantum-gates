%%SWAP Gate simulation based on Volkov and Salikhov 2011
%Lena Aly - NYU Abu Dhabi 
%lwa2019@nyu.edu

clear; clc; close all;

global dt start rho_t NoP_90 NoP_180 NoP_evolution rho error 


addpath(genpath('/path/to/spinach'));

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
error = 19; %in degrees 

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
Ix=Lx+Sx;
Iy=Ly+Sy;
Iz=Lz+Sz;

% Time parameters
dt = (1/J )*0.1;
flip_1 = 0.5* pi;   % 90 Rotation angle
flip_2 = 1*pi;      % 180

t_90 = flip_1 / amplitude_mw;
t_180 = flip_2 / amplitude_mw;          


NoP_90 = round(t_90 / dt);
NoP_180  = round(t_180 / dt);
NoP_evolution = round ((pi/ (2*J))/dt);
start =0;
CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];

% Initial pure state density matrix 
psi1 = [0; 1];
psi0 = [1; 0];
r =  [0;1;0;0];
rho = r*r';
rho_t = zeros(4, 4, 0);  % ensemble tracking
disp("After step definition:");
disp(rho);


function pulse( theta, axis, H, direction)


global dt start rho_t NoP_90 NoP_180 rho;

    U = expm(-1i * H * dt);

    axis_list =["1y", "2y", "1x", "2x", "12x", "12y"];
    if theta == 90 && ismember(axis, axis_list)
  
        NoP_segment = NoP_90;


    elseif theta == 180 &&  ismember(axis, axis_list)
        NoP_segment = NoP_180;

    else
        print ("try again; you have not entered a valid axis")

    end 

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

function evolution(H)
global dt NoP_evolution start rho_t rho;

U = expm(-1i * H * dt);
NoP_segment = NoP_evolution;

 for i= start +1 : start +NoP_segment
        rho_t(:, :, i) = rho;
           rho = U*rho*U'; 

 end 

 start = start + NoP_segment;
 rho_t(:,:,start) = rho;
end

disp("After step funtions:");
disp(rho);

% Pulse Sequence 


H_evol = J*(Lz*Sz+ Lx*Sx + Ly*Sy);

%180 -12x
H = amplitude_mw* Ix;
pulse(180, "12x", H, -1);
disp("After step 2:");
disp(rho);


%L 
evolution( H_evol);

H = amplitude_mw* Ix;
pulse(180, "12x", H, 1);
disp("After step 2:");
disp(rho);

evolution( H_evol);

H = amplitude_mw* Ix;
pulse(180, "12x", H, -1);
disp("After step 2:");
disp(rho);

%_________________________

%180 12x
H = amplitude_mw* Ix;
pulse(180, "12x", H, 1);
disp("After step 2:");
disp(rho);

%180 -12y
H = amplitude_mw* Iy;
pulse(180, "12y", H, -1);
disp("After step 2:");
disp(rho);

%L 
evolution( H_evol);

H = amplitude_mw* Ix;
pulse(180, "12x", H, 1);
disp("After step 2:");
disp(rho);

evolution( H_evol);

H = amplitude_mw* Ix;
pulse(180, "12x", H, -1);
disp("After step 2:");
disp(rho);

%_________________________

%180 12y
H = amplitude_mw* Iy;
pulse(180, "12x", H, 1);
disp("After step 2:");
disp(rho);


%L 
evolution( H_evol);

H = amplitude_mw* Ix;
pulse(180, "12x", H, 1);
disp("After step 2:");
disp(rho);

evolution( H_evol);

H = amplitude_mw* Ix;
pulse(180, "12x", H, -1);
disp("After step 2:");
disp(rho);

%_________________________

disp("After step end:");
disp(rho);

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

trace(rho^2)
time = (0:total_steps-1) * dt;

%% Plot probabilities
figure(1);
plot(time, P0, 'r');hold on;plot(time, P2, 'b');...
hold on;plot(time, P3, 'g');hold on;plot(time, P1, 'black');
title('Population evolution under SWAP');
xlabel('Time (a.u.)');
ylabel('Probability');
legend('00', '01', '10', '11');
grid on;


%Fidility Calculation
psi_target = SWAP * r;            
rho_target = psi_target * psi_target'; 
Fid = real(trace(rho * rho_target));
Fid