clear all
close all
clc

% Cf = 155494.663 ;  % N / Rad
% Cr = 155494.663 ;  % N / rad 
% Iz = 1436.24 ; %m kgm^2
% Lf = 1.165 ;  % m
% Lr = 1.165 ; % m 
% m = 1140 ; % kg
% Vx = 18 ;  



%%% Rajamani Parameters
Cf = 80000 ;  % N / Rad
Cr = 80000 ;  % N / rad 
Iz = 2873 ; %m kgm^2
Lf = 1.1 ;  % m
Lr = 1.58 ; % m 
m = 1573 ; % kg
Vx = 30 ;  

% 100kph= 27.77778m/s
% 30kph= 8.333333m/s
 
Np = 12; % Predict Horizon
Nc = 3;  % Control Horizon
 
 %%%% Dynamic Models in terms of error with respect to road 
 % states : distance of the line center, distance rate of line center,
 % orientation error, orientation error rate
 % inputs : desired yaw rate(Vx/R), steering angle
 
 
A = [0 1 0 0;
     0, -(Cf + Cr)/(m * Vx), ( Cf + Cr) / m , (- Cf * Lf + Cr * Lr) / (m * Vx) ;
     0 0 0 1;
     0, -( Cf * Lf -  Cr * Lr)/(Iz * Vx), ( Cf * Lf - Cr * Lr)/Iz, -( Cf * Lf^2 +  Cr * Lr^2)/(Iz * Vx)];
 
 
 B = [ 0 0 ;
      (Cf/m), (-( Cf * Lf -  Cr * Lr)/(m * Vx)) - Vx ;
       0 0 ;
       (Cf * Lf) / Iz, -( Cf * Lf^2 + Cr * Lr^2)/(Iz * Vx)];

C = [1,0,0,0;
     0,0,1,0] ;
 
D= [0 0;0 0] ;

%%%% State Feedback
P = [-5-3i; -5+3i ; -7 ; -10] ; %% desired poles
B1 = B(:,1) ; %% (A-B1K)x + B2desYawRate
K = place(A,B1,P) ;

A = (A - B1*K) ;


%%%%%% Parameters for sim
rw = 0.1 ;
ts = 0.01 ;
K = 1/1000 ;

N_sim = 421 ; % # of sample

xm = [1.5; 0 ; deg2rad(-30); 0]; %initial states 
Xf = [0;0;0;0;1.5 ; deg2rad(-30) ] ;  %% delta states

u = [0];
y = [0;0];

ref_signal = zeros(2,N_sim) ;
ref_signal(1,:) = 0 * ones(N_sim,1) ;  % referenece signal
ref_signal(2,:) = 0 * ones(N_sim,1) ;  % referenece signal
file = matfile("frenet.mat") ;
u2 = (Vx .* K) * ones(1,N_sim);
%  u2 = (Vx .* file.rk);

%%% for constraints
u0(1) = [0];

amplitude_constraint = 0.7853981634; %45 derece 
% amplitude_constraint = 1.0471975512; % 60 derece
rate_constraint = 0.2 ;

%%%%
[system.Ad, system.Bd, system.Cd, system.Dd] = c2dm(A,B,C,D,ts); % Model Discritization
system.Bd1 = system.Bd(:,1) ;
system.Bd2 = system.Bd(:,2) ;

[ aug_system.A_e, aug_system.B_e, aug_system.Bw_e, aug_system.C_e ] = generate_aug_model(system); % Augmented Model (adding integral action)

[F, phi, omega, phiT_Rs, phiT_phi, phiT_F, phiT_omega] = calculate_mpc_gains( aug_system, Nc, Np ) ;

%%% simulation
% [y,u, deltaUs] = unconstrained_mpc_simulation(system, aug_system, xm, Xf, F, phi, phiT_Rs, Nc, Np, ref_signal,rw, N_sim,u,y, u2) ;
[y,u, deltaUs] = constrained_mpc_simulation(system, aug_system, xm, Xf, F, phi, phiT_Rs, Nc, Np, ref_signal,rw, N_sim,u,y, u2, u0, rate_constraint, amplitude_constraint, omega) ;

%%%%%%%% PLOTTING 
k = 0 : (N_sim -1);
aaa(1,k+1) = rad2deg(amplitude_constraint);
rrr(1,k+1) = rad2deg(rate_constraint) ;

figure
hold on
plot(k,y(1,:))
plot(k,ref_signal(1,:))
ylabel(' distance of the vehicle from center line[m] ')
legend('system response ',' ref signal ')
grid on

figure
hold on
plot(k,rad2deg(y(2,:)))
plot(k,ref_signal(2,:))
ylabel(' orientation error respect to road[deg] ')
legend('system response ',' ref signal ')
grid on

figure
hold on
plot(k,rad2deg(u(1,:)))
ylabel(' steering angle[deg]')
plot(k,aaa , "r --" )
plot(k,-aaa , "r --" )
legend('steering angle ',' constraint ')
grid on

figure
hold on
plot(k,rad2deg(deltaUs(1,:)))
ylabel(' steering angle rate ')
plot(k,rrr , "r --" )
plot(k,-rrr , "r --" )
legend('steering angle rate',' constraint ')
grid on

figure
plot(k,rad2deg(u2(1,:)))
ylabel(' desired yaw rate, determined from R (deg / sec) ')
grid on

%%%% Frenet
figure
plot(file.x, file.y, 'xb', 'DisplayName', 'input');
hold on
plot(file.rx, file.ry, '-r', 'DisplayName', 'spline');
grid on;
%axis equal;
xlabel('x[m]');
ylabel('y[m]');
legend;

figure
plot(file.s, rad2deg(file.ryaw), '-r', 'DisplayName', 'yaw');
grid on;
legend;
xlabel('line length[m]');
ylabel('yaw angle[deg]');

figure
plot(file.s, file.rk, '-r', 'DisplayName', 'curvature');
grid on;
legend;
xlabel('line length[m]');
ylabel('curvature [1/m]');








