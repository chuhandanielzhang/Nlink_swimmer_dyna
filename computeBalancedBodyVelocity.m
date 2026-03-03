function [ xi ] = computeBalancedBodyVelocity( alpha, d_alpha, beta, d_beta, leg_posi, L, Lleg, activation, K, F,force_model )

if sum(activation)==0
    xi=[0;0;0;0];
else
    
    xi0 = computeBodyVelocity( alpha, d_alpha, beta, d_beta, leg_posi, L, Lleg, activation, K, F );
    
    xi00=xi0;  % x-y coordinate
    xi=norm(xi0(1),xi0(2));
    xi0(2)=atan2(xi0(1),xi0(2));
    xi0(1)=xi; % cylindrical coordinate
    option=optimset('Display','off');
    [xi1,a1,exitflag] = fsolve(@(xi) (forceResidual( xi, alpha, d_alpha, beta, d_beta, leg_posi, L, Lleg, activation,force_model)),...
        xi0,option);
    
    xii=xi1(1);
    xi1(1)=xii*sin(xi1(2));
    xi1(2)=xii*cos(xi1(2));
    
    if exitflag ~= 1
        [xi2,a2,exitflag] = fsolve(@(xi) (forceResidual_1( xi, alpha, d_alpha, beta, d_beta, leg_posi, L, Lleg, activation,force_model)),...
            xi00,option);
        if exitflag ~=1
            display(['Non-convergence Exit Flag:',num2str(exitflag)]);
        else
            xi=xi2;
        end
    else
        xi=xi1;
        a2=inf;
        xi2=zeros(size(xi));
    end
    
    if exitflag ~= 1
        xi = computeBodyVelocity( alpha, d_alpha, beta, d_beta, leg_posi, L, Lleg, activation, K, F );
        xi(4)=1;
    else
        xi(4)=0;
    end
end
end