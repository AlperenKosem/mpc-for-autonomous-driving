function [F, phi, omega, phiT_Rs, phiT_phi, phiT_F, phiT_omega] = calculate_mpc_gains(aug_system, Nc, Np )
    
    A_aug = aug_system.A_e ;
    B_aug = aug_system.B_e ;
    C_aug = aug_system.C_e ;
    Bw_aug = aug_system.Bw_e ;
    [output_size, X ] = size(C_aug);
    [X , number_of_inputs] = size(B_aug);

    aug_state_size = size(A_aug); % augmented state size

    [Fi, Fj] = size(C_aug * A_aug);

    F = zeros(Np * Fi ,Fj);
    phi = zeros(Np * Fi , Nc * number_of_inputs); %Phi = Np * Fi X Nc *number of inputs

    if output_size == 2 && number_of_inputs == 2
        %%%% calculating F
        m = 1 ;
        j = 1;
        pow = 1;

        for i = 1 : 2 : Np * Fi 

           dummy = C_aug * A_aug^pow ;

           while j <= aug_state_size(2)        
               for k = 0 : output_size - 1
                    F(i+k, j) = dummy(k + 1,j);
               end
                j = j + 1;
           end      
           j = 1 ;
           pow = pow +1;
        end
        
        %%%% calculating phi
        m = 0 ;
        j = 1 ;
        n = 0 ;
        while j <=  Nc * number_of_inputs
            for i = 1 : 2 : Np * Fi

                dummy = C_aug * A_aug^(m) * B_aug ;
                phi(n + i:n + i+1,j:j+1) = dummy(1:2,:);


                m = m + 1 ;
            end
            j = j + 2 ;
            m = 0 ;
            n = n + 2 ;
        end
        
         phi = phi(1:end - ((Nc * Fi) -2) ,:) ; 
         phiT_phi = phi' * phi ;
         phiT_F = phi' * F  ;

         R_bar = phiT_F(:,end + 1 - output_size:end) ;
         phiT_Rs = R_bar;  % may add Q weight matrix

    end
    
    if output_size == 2 && number_of_inputs == 1
        
                %%%% calculating F
        m = 1 ;
        j = 1;
        pow = 1;

        for i = 1 : 2 : Np * Fi 

           dummy = C_aug * A_aug^pow ;

           while j <= aug_state_size(2)        
               for k = 0 : output_size - 1

                    F(i+k, j) = dummy(k + 1,j);

               end
        %         F(i, j) = dummy(1,j);
        %         F(i + 1 , j) = dummy(2, j)

                j = j + 1;
           end      
           j = 1 ;
           pow = pow +1;
        end

        %%%% calculating phi
        i = 1 ; 
        k = 0 ;
        m = 1 ;
        for j = 1 : Nc * number_of_inputs
            while (i <= Np)

        %         phi(i,j) = Cnew*Anew^(i-1 - k) * Bnew;
                  dummy = C_aug * A_aug^(i-1 - k) * B_aug;
                  dummy2 = C_aug * A_aug^(i-1 - k) * Bw_aug;
                  
                  phi(m,j) = dummy(1);
                  phi(m + 1, j) = dummy(2);
                  omega(m,j) = dummy2(1);
                  omega(m + 1, j) = dummy2(2);
                  
                  i = i + 1;    
                  m = m + 2;
            end 
           i = j + 1 ;
           m = j + i ;
           k = k + 1 ;
        end

        phiT_phi = phi' * phi ;
        phiT_F = phi' * F  ;
        phiT_omega = phi' * omega ;
        phiT_Rs = phiT_F(:,end + 1 - output_size:end) ;  % last # of output of phi'*F columns ;
        
        
    end
         
         
         
end