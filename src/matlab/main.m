clear all
close all
clc

parameters.Cf = 155494.663 ;  % N / Rad
parameters.Cr = 155494.663 ;  % N / rad 
parameters.Iz = 1436.24 ; %m kgm^2
parameters.Lf = 1.165 ;  % m
parameters.Lr = 1.165 ; % m 
parameters.m = 1140 ; % kg
parameters.Vx = 30 ;  



%%% Rajamani Parameters
% parameters.Cf = 2 * 80000 ;  % N / Rad  % cornering stiffness of front tires
% parameters.Cr = 2 * 80000 ;  % N / rad  % cornering stiffness of rear tires
% parameters.Iz = 2873 ; %m kgm^2     % aw-plane rotational inertia 
% parameters.Lf = 1.1 ;  % m          % distance from C.G. to front axle 
% parameters.Lr = 1.58 ; % m          % distance from C.G. to rear axle 
% parameters.m = 1573 ; % kg          % vehicle mass 
% parameters.Vx = 18 ;  % m/sec       % longitudinal velocity

% 100kph= 27.77778m/s
% 30kph= 8.333333m/s
 
Np = 12; % Predict Horizon
Nc = 3;  % Control Horizon
 
 %%%% Dynamic Models in terms of error with respect to road 
 % states : distance of the line center, distance rate of line center,
 % orientation error, orientation error rate
 % inputs : desired yaw rate(Vx/R), steering angle
 
 
A = [0 1 0 0;
     0, -(parameters.Cf + parameters.Cr)/(parameters.m * parameters.Vx),...
     ( parameters.Cf + parameters.Cr) / parameters.m , (- parameters.Cf * parameters.Lf + parameters.Cr * parameters.Lr) / (parameters.m * parameters.Vx) ;
     0 0 0 1;
     0, -( parameters.Cf * parameters.Lf -  parameters.Cr * parameters.Lr)/(parameters.Iz * parameters.Vx), ...
     ( parameters.Cf * parameters.Lf - parameters.Cr * parameters.Lr)/parameters.Iz, -( parameters.Cf * parameters.Lf^2 +  parameters.Cr * parameters.Lr^2)/(parameters.Iz * parameters.Vx)];
 
 
 B = [ 0 0 ;
      (parameters.Cf/parameters.m), (-( parameters.Cf * parameters.Lf -  parameters.Cr * parameters.Lr)/(parameters.m * parameters.Vx)) - parameters.Vx ;
       0 0 ;
       (parameters.Cf * parameters.Lf) / parameters.Iz, -( parameters.Cf * parameters.Lf^2 + parameters.Cr * parameters.Lr^2)/(parameters.Iz * parameters.Vx)];

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
% parameters.curv = 1 / 250 ; %meter



xm = [1.5; 0 ; deg2rad(-30); 0]; %initial states 
Xf = [0;0;0;0;1.5 ; deg2rad(-30) ] ;  %% delta states

% xm = [0; 0 ; 0; 0]; %initial states 
% Xf = [0;0;0;0;0 ; 0 ] ;  %% augmented states




file = matfile("frenet.mat") ;
N_sim = length(file.rk) ; % # of sample
% N_sim = 300 ;
u = [0];
y = [0;0];

ref_signal = zeros(2,N_sim) ;
ref_signal(1,:) = 0 * ones(N_sim,1) ;  % referenece signal
ref_signal(2,:) = 0 * ones(N_sim,1) ;  % referenece signal


% u2 = (parameters.Vx .* parameters.curv) * ones(1,N_sim); % constant curvature
u2 = (parameters.Vx .* file.rk);
curvature = file.rk ;

%%% for constraints
u0(1) = [0];

amplitude_constraint = deg2rad(35); %30 derece 
% amplitude_constraint = 0.7853981634; %45 derece 
% amplitude_constraint = 1.0471975512; % 60 derece

rate_constraint = deg2rad(15) ;

%%%%
[system.Ad, system.Bd, system.Cd, system.Dd] = c2dm(A,B,C,D,ts); % Model Discritization
system.Bd1 = system.Bd(:,1) ;
system.Bd2 = system.Bd(:,2) ;

[ aug_system.A_e, aug_system.B_e, aug_system.Bw_e, aug_system.C_e ] = generate_aug_model(system); % Augmented Model (adding integral action)

[F, phi, omega, phiT_Rs, phiT_phi, phiT_F, phiT_omega] = calculate_mpc_gains( aug_system, Nc, Np ) ;

%%% simulation
% [y,u, deltaUs] = unconstrained_mpc_simulation(system, aug_system, xm, Xf, F, phi, phiT_Rs, Nc, Np, ref_signal,rw, N_sim,u,y, u2) ;
[y,u, deltaUs] = constrained_mpc_simulation(system, aug_system, xm, Xf, F, phi, phiT_Rs, Nc, Np, ref_signal,rw, N_sim,u,y, u2, u0, rate_constraint, amplitude_constraint, omega, K, parameters) ;


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
plot(k, rad2deg(u2(1,:)))
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
% 
% figure
% plot(file.s, rad2deg(file.ryaw), '-r', 'DisplayName', 'yaw');
% grid on;
% legend;
% xlabel('line length[m]');
% ylabel('yaw angle[deg]');
% 
figure
plot(file.s, file.rk, '-r', 'DisplayName', 'curvature');
grid on;
legend;
xlabel('line length[m]');
ylabel('curvature [1/m]');








