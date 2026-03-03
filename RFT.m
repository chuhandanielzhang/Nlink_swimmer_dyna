function [ wrench ] = RFT( xi,force_model )

%input a body velocity and return the resistive force of a small segment
if force_model == 0 % poppy seeds
    A = 0.27;
    B = -0.32;
    C = 0.66;
    Fto = 0.09;
    len = norm(xi(1:2));
    phi=abs(atan(xi(2)/xi(1)));
    sig1=2/(1+exp(5*xi(1)))-1;
    sig2=2/(1+exp(5*xi(2)))-1;
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = [sig1*((A*cos(phi)+B*(1-sin(phi)))+Fto);sig2*C*sin(phi);0];%no torque term
    end
elseif force_model == -1 % wheels
    a = 0.0011;com_ave
    b = 0.47;
    c = 1.1;
    d = -0.19;
    p1 = -7.8e-9;
    p2 = 1.2e-6;
    p3 = -8.7e-5;
    p4 = 0.003;
    p5 = 0.097;
    len = norm(xi(1:2));
    phi = rad2deg(abs(atan(xi(2)/xi(1))));
    sig_exp = 1e4;
    sig1 = 2/(1+exp(sig_exp*xi(1)))-1;
    sig2 = 2/(1+exp(sig_exp*xi(2)))-1;
    F_perp = @(phi) a*phi + b./(1+exp(-phi/c)) + d; % From Jennifer's PRE
    F_para = @(phi) p1*phi.^4 + p2*phi.^3 + p3*phi.^2 + p4*phi + p5;
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = [sig1*F_para(phi);sig2*F_perp(phi);0];%no torque term
    end
elseif force_model == -2 % viscous fluid
    len = norm(xi(1:2));
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = -diag([1;1;0])*xi;%no torque term
    end
elseif force_model == -3 % BB
    Cs = 3.21;
    Cf = 1.34;
    Cl = -0.82;
    gamma = 2.79;
    len = norm(xi(1:2));
    phi=abs(atan(xi(2)/xi(1)));
    sig1 = 2/(1+exp(5*xi(1)))-1;
    sig2 = 2/(1+exp(5*xi(2)))-1;
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = [sig1*(Cf*cos(phi)+Cl*(1-sin(phi)));sig2*Cs*sin(atan(gamma*sin(phi)));0];%no torque term
    end
elseif force_model == -3.5 % BB, head-cap
    Cs = 0.73;
    Cf = 0.30;
    Cl = -0.19;
    gamma = 0.52;
    len = norm(xi(1:2));
    phi=abs(atan(xi(2)/xi(1)));
    sig1 = (2/(1+exp(10*xi(1)))-1);
    sig2 = (2/(1+exp(10*xi(2)))-1);
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = [sig1*(Cf*cos(phi)+Cl*(1-sin(phi)));sig2*Cs*sin(atan(gamma*sin(phi)));0];%no torque term
    end
elseif force_model == 2
    A = 2;
    len = norm(xi(1:2))*2*pi; % 2*pi for freq
    phi=abs(atan(xi(2)/xi(1)));
    sig1 = sign(xi(1));
    sig2 = sign(xi(2));
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = [sig1*cos(phi);sig2*A*sin(phi);0]*len;%no torque term
    end
elseif force_model == 3  % aqua
    A = 2;
    len = norm(xi(1:2))*2*pi; % 2*pi for freq
    phi=abs(atan(xi(2)/xi(1)));
    sig1 = sign(xi(1));
    sig2 = sign(xi(2));
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = [sig1*cos(phi);sig2*A*sin(phi);0]*len^2;%no torque term
    end
elseif force_model == -4  % aqua
    parameters.r = 0.12; % animal radius in m
    parameters.z = 0.12; % depth in material in m
    parameters.paraMult = 1;%1*parameters.mu/0.2; % coefficient for Al plate = 0.2;
    parameters.perpMult = 1;%1;
    parameters.fbelly = 0;
    A = 685.1;
    B = 0.00013;
    C = 0.25;

    % para fit coefficients (for angles in radians)
    a0 = 24.85;
    a1 = 103.6;
    b1 = 58.9;
    a2 = -22.29;
    b2 = -11.79;
    c = 1.502;

    d = parameters.z;
    w = parameters.r*2;

    f_perp = @(x) A*w*d*(1+C./sqrt(B+sin(x).^2)).*sin(x);
    f_para = @(x) w*d*(a0+a1*cos(c*x)+b1*sin(c*x)+a2*cos(2*c*x)+b2*sin(2*c*x));

    len = norm(xi(1:2));


    vHat = xi(1:2)/len; %[cos(psi); sin(psi)];
    sig1 = (2/(1+exp(1e2*xi(1)))-1);
    sig2 = (2/(1+exp(1e2*xi(2)))-1);
    cosPsi = vHat(1); % v_hat dot t_hat = cos(psi)
    sinPsi = vHat(2); % v_hat dot n_hat = cos(pi/2-psi) = sin(psi)
    phi = abs(atan(vHat(2)/vHat(1))); % map phi back into first quadrant

    % Perrin's code does this (though it seems that -pi/2 <= psi <= pi/2):
    %ft = f_para(abs(acos(cosPsi)));%.*Tang(1);
    %fn = f_perp(Psi);%.*Norm(2);%A*w*d.*(1 + C./ sqrt(B + sinPsi.^2)).*sinPsi;

    % Baxi's code is consistent with this:
    %ft = sign(cosPsi)*(f_para(phi));
    %fn = sign(sinPsi)*(f_perp(phi));

    ft = sig1*abs(f_para(phi));
    fn = sig2*abs(f_perp(phi));

    wrench = -[ft*parameters.paraMult;fn*parameters.perpMult;0];%no torque term
elseif force_model == -5
    len = norm(xi(1:2));
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = -1 * diag([1;1;0]) * xi / len ;
    end

else
    A = force_model;
    len = norm(xi(1:2))*2*pi; % 2*pi for freq
    phi=abs(atan(xi(2)/xi(1)));
    sig1 = (2/(1+exp(1e1*xi(1)))-1);
    sig2 = (2/(1+exp(1e1*xi(2)))-1);



    % f = @(x) (2 ./ (1 + exp(-1e1 * x))) - 1;
    % 
    % fplot(f, [-1, 1]);
    % title('S');
    % xlabel('x');
    % ylabel('y');
    % legend('y = (2/(1+exp(-10x)))-1', 'Location');
    if len == 0
        wrench = zeros(3,1);
    else
        wrench = [sig1*cos(phi);sig2*A*sin(phi);0];%no torque term
    end
end


end

