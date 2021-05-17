function [y,u0, deltaUs] = constrained_mpc_simulation(system, aug_system, xm, Xf, F, phi, phiTRs, Nc, Np, ref_signal,rw, N_sim,u,y, u2, u0, rate_constraint,amplitude_constraint, omega, K, parameters)
      
    Ad = system.Ad ;
    Bd1 = system.Bd1 ; %use first column of Bd matrix
    Bd2 = system.Bd2 ; %use first column of Bd matrix
    Cd = system.Cd ;

    u2(end + 1 : (numel(u2) + Nc)) = 0 ;  
     
    deltaW = zeros(Nc,1);
              
    
    for kk = 2 : N_sim 
        
        %%% calculation feedforward term
%         Kv = parameters.Lr * parameters.m / (2 * parameters.Cf * (parameters.Lr + parameters.Lf)) - parameters.Lf * parameters.m / (2 * parameters.Cr * (parameters.Lr + parameters.Lf)) ;
%         ay = parameters.Vx^2 * parameters.curv ;
%         u_ff =  (parameters.Lr + parameters.Lf) * parameters.curv + Kv * ay - K(1,3) * (parameters.Lr * parameters.curv - ((parameters.Lf / (2 * parameters.Cr)) * parameters.m * parameters.Vx^2 * parameters.curv / (parameters.Lf + parameters.Lr)));
%         
        
        
        for i = 1 : Nc
            deltaW(i,1) = u2(kk + i - 1) - u2(kk - 1 + i - 1) ;       
        end
       
       %%% constraints amplitude and rate of change
        for i = 1 :  Nc
            for j = 1 : Nc
                if i == j 
                    M(i , j) = 1 ;
                    gamma(i,1) =  rate_constraint;
                    M(i + Nc, j) = -1 ;
                    gamma(i + Nc, 1) = rate_constraint ;
                end
                
                if i >= j
                    M(i + 2 * Nc , j) = 1 ;
                    M(i + 3 * Nc , j) = -1 ;
                    gamma(i + 2 * Nc ,1) =  amplitude_constraint - u0(kk-1);
                    gamma(i + 3 * Nc ,1) = amplitude_constraint + u0(kk-1) ;
                end    
            end
        end
       
        E =  2 * (phi' * phi + rw * eye(Nc)) ;
        Fconstraint = - 2 *  (phiTRs * ref_signal(:,kk)  - phi' * F * Xf - phi' * omega * deltaW)  ; 
        deltaU = Qphild(E,Fconstraint,M,gamma) ;

            
         deltaUs(kk) = deltaU(1);
         u = u + deltaU(1); %+ u_ff ;
         u0(kk) = u;
         xm_old = xm ;
         xm = Ad * xm + Bd1 * u  +  Bd2 * u2(:,kk);
         y(:,kk) = Cd * xm ;

         Xf(1:4) = xm - xm_old ;
         Xf(5:6) = y(1:2,kk);
    end

end