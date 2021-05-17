function eta = Qphild(H,f,M,gamma)

% E = H ;
% F = f ;
% M = A_cons ;
% gamma = b;
% 
% eta = x ;

[n1,m1] = size(M);
eta = - H \ f ;
% eta = - f^-1 * H;
kk = 0;
for i = 1 : n1
    
    if(M(i,:)* eta > gamma(i))
        kk = kk + 1 ;
    
    end
end

if(kk == 0)
    return
end

P = M * (H \ M' );
d = (M * (H \ f) + gamma);
[n,m] = size(d) ;

x_ini = zeros(n,m) ;
lambda = x_ini ;
% al = 10;

for km = 1 : 38
    lambda_p = lambda ;
    for i = 1 : n
        w = P(i,:)* lambda - P(i,i)* lambda(i,1);
        w = w + d(i,1);
        la = -w / P(i,i);
        if max(la) < 0 
            lambda(i,1) = 0;
        else
            lambda(i,1) = max(la);
        end
        

    end
    al = (lambda - lambda_p)' * (lambda - lambda_p);
    if al < 10e-8
        break;
    end
end

eta = - H \ f - H \ M' * lambda;

end