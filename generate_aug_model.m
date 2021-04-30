function [ A_e, B_e, Bw_e, C_e ] = generate_aug_model(system)
        
    Ad = system.Ad ;
    Bd = system.Bd1 ;
    Cd = system.Cd ;
    Bdw = system.Bd2 ;
    Bw_e = 0; 
    %%%%%%%%%%%%%%%%
    %Augment state equations
    %%%%%%%%%%%%%%%%

    [ m1,n1 ]=size(Cd);
    [n1,n_in]=size(Bd);

    A_e = eye(n1+m1,n1+m1);
    A_e(1:n1,1:n1) = Ad;
    A_e(n1+1:n1+m1,1:n1) = Cd*Ad;
    B_e = zeros(n1+m1,n_in);

    B_e(1:n1,:) = Bd;
    B_e(n1+1:n1+m1,:) = Cd * Bd;
    
    Bw_e(1:n1,:) = Bdw;
    Bw_e(n1+1:n1+m1,:) = Cd * Bdw;
    
    C_e = zeros(m1,n1+m1) ;
    C_e(:,n1+1:n1+m1) = eye(m1,m1);
    
end
