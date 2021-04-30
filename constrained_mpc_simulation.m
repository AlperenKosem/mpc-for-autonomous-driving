function [y,u0, deltaUs] = constrained_mpc_simulation(system, aug_system, xm, Xf, F, phi, phiTRs, Nc, Np, ref_signal,rw, N_sim,u,y, u2, u0, rate_constraint,amplitude_constraint, omega)
    
    A_aug = aug_system.A_e ;
    B_aug = aug_system.B_e ;
    C_aug = aug_system.C_e ;
    
   
    Ad = system.Ad ;
    Bd1 = system.Bd1 ; %use first column of Bd matrix
    Bd2 = system.Bd2 ; %use first column of Bd matrix
    Cd = system.Cd ;
    [row_phi,X ] = size(phi);
    u2(end + 1 : (numel(u2) + Nc)) = 0 ;  
     
    deltaW = zeros(Nc,1);
              
    
    for kk = 2 : N_sim 
        
        
        deltaW(1,1) = u2(kk) - u2(kk - 1) ;
        deltaW(2,1) = u2(kk + 1) - u2(kk) ;
        deltaW(3,1) = u2(kk + 2) * u2(kk + 1) ;
       %%% constraints amplitude and rate of change
        M = [1 0 0  ;
            -1 0 0 ;
            0 1 0 ;
            0 -1 0 ;
            0 0 1  ;
            0 0 -1 ;
            1 0 0  ;
            1 1 0 ;
            1 1 1 ;
            -1 0 0 ;
            -1 -1 0  ;
            -1 -1 -1];
%         M(13 : 12+row_phi,:) =  phi(:,:) ;
%         output_constraint = 5;
%         M(13 + row_phi : 12 + row_phi + row_phi ,:) = - phi(:,:) ;
        
    %     gamma = [3 ; 1.5; 3 ; 1.5 ; 3 ; 1.5; 1 - u1(kk-1) ; 1 - u1(kk-1) ; 1 - u1(kk-1) ; 1.5 + u1(kk-1) ; 1.5 + u1(kk-1)  ; 1.5 + u1(kk-1) ] ;

        gamma = [rate_constraint ; rate_constraint ; rate_constraint ; rate_constraint ; rate_constraint ; rate_constraint ; 
                 amplitude_constraint - u0(kk-1) ; amplitude_constraint - u0(kk-1) ; amplitude_constraint - u0(kk-1) ; 
                 amplitude_constraint + u0(kk-1) ; amplitude_constraint + u0(kk-1)  ; amplitude_constraint + u0(kk-1) ] ;
%         gamma(13 : 12 + row_phi) =  output_constraint - F * Xf ;
%         gamma(13 + row_phi : 12 + row_phi + row_phi) = - (output_constraint - F * Xf);
    %       Constraints on amplitude of Control
    %     M = [1 0 0  ;
    %     1 1 0 ;
    %     1 1 1 ;
    %     -1 0 0 ;
    %     -1 -1 0  ;
    %     -1 -1 -1 ];
    % 
    %     gamma = [1 - u1(kk-1) ; 1 - u1(kk-1) ; 1 - u1(kk-1) ; 1.5 + u1(kk-1) ; 1.5 + u1(kk-1)  ; 1.5 + u1(kk-1) ] ;

        E =  2 * (phi' * phi + rw * eye(Nc)) ;
        Fconstraint = - 2 *  (phiTRs * ref_signal(:,kk)  - phi' * F * Xf - phi' * omega * deltaW)  ; 
        deltaU = Qphild(E,Fconstraint,M,gamma) ;

            
         deltaUs(kk) = deltaU(1);
         u = u + deltaU(1) ;
         u0(kk) = u ;
         xm_old = xm ;
         xm = Ad * xm + Bd1 * u  +  Bd2 * u2(:,kk);
         y(:,kk) = Cd * xm ;

         Xf(1:4) = xm - xm_old ;
         Xf(5:6) = y(1:2,kk);
    end

end