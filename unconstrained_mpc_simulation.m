function [y,u1, deltaUs] = unconstrained_mpc_simulation(system, aug_system,xm,Xf, F ,phi,phiTRs,Nc,Np,ref_signal,rw, N_sim,u, y, u2)
%     ref_signal
    
    A_aug = aug_system.A_e ;
    B_aug = aug_system.B_e ;
    C_aug = aug_system.C_e ;
   
    Ad = system.Ad ;
    Bd1 = system.Bd1 ; %use first column of Bd matrix
    Bd2 = system.Bd2 ; %use first column of Bd matrix
    Cd = system.Cd ;
    
    [output_size, X ] = size(C_aug);
    [X,input_size] = size(Bd1) ;
    
    for kk = 1 : 1 : N_sim
        
         deltaU =  (phi' * phi + rw * eye(Nc * input_size) )^-1  * (phiTRs * ref_signal(:,kk)  - phi' * F * Xf) ;
         deltaUs(:,kk) = deltaU(1);
         u = u + deltaU(1) ;
         u1(:,kk) = u ;

         xm_old = xm ;
         xm = Ad * xm + Bd1 * u ; %  +  Bd2 * u2(:,kk) ;
         y(:,kk) = Cd * xm ;
        %      
         Xf(1:4) = xm - xm_old ;
         Xf(5:6) = y(1:2,kk);

    end   
end