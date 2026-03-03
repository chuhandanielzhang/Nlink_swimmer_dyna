%%
N = 6; % Number of link
% 1 for Coulomb friction; 0 for viscous
isRFT = 1;
% Number of gait cycles
period = 20; 
% lambda = 0.5 % number of waves
% Mapping from velocity to force on each link
K = diag([1;1;1]);
% amp = pi/8; %Amplitude
L = 2.3;
basis_paras.frac = 1; %fraction of body in contact
freq = 2*pi/5; %
% % % % % % ================================================================
basis_paras.lfs = 1 - 1/N;
basis_paras.spatial_phase = deg2rad(45);
basis_paras.Amp1 = deg2rad(20);
basis_paras.Amp2 = deg2rad(20);
paras.leg_posi = 1:N;
% paras.leg_posi = [1, 2, 3, 4];
% paras.leg_posi = [];
% paras.leg_posi = [];
paras.num_leg = length(paras.leg_posi);
paras.Aleg = deg2rad(10);
paras.Lleg = 4.5;
paras.dutyf = 0.5;
paras.symmmode = 1;
basis1 = sin((0:basis_paras.lfs:basis_paras.lfs*(N-2))*2*pi + basis_paras.spatial_phase)';
basis2 = cos((0:basis_paras.lfs:basis_paras.lfs*(N-2))*2*pi + basis_paras.spatial_phase)';
phase_max = 0;
phase = @(t) atan2(basis_paras.Amp1*sin(freq*t), basis_paras.Amp2*cos(freq*t)) - phase_max;

%%
beta0 = [pi/2*ones(1,paras.num_leg); -pi/2*ones(1,paras.num_leg)];
beta0 = beta0(:);
beta   = @(c) get_beta(c, paras.num_leg, paras.Aleg, paras.dutyf, basis_paras.lfs, paras.symmmode) + beta0;
d_beta = @(c) (beta(c + 0.01) - beta(c)) / 0.01;
activation = @(c) [ones(1, N)*0, get_act(c, basis_paras.lfs, paras.num_leg, paras.dutyf, paras.symmmode)];
%%1
%++++++++++

alpha   = @(t)  basis_paras.Amp1*sin(freq*t)*basis1 - basis_paras.Amp2*cos(freq*t)*basis2;
d_alpha = @(t)  basis_paras.Amp1*freq*cos(freq*t)*basis1 + basis_paras.Amp2*freq*sin(freq*t)*basis2;
%%
% current time
t = 0; 
T = 100;
dt = 2*pi/T;


ave_body = -pi/2;
g_h=eye(3)*[cos(ave_body) -sin(ave_body), 0;
    sin(ave_body)  cos(ave_body), 0;
    0              0, 1];

ave_body=sum(cumsum(alpha(t)))/N; % get average body rotation
[k, ~] = framesInHead(alpha(t), beta(phase(t)), paras.leg_posi, L, paras.Lleg); % find the initial body frames

com_ave = computeCOM(k); % find the CoM from k
gi_ave_body=g_h*[cos(ave_body) -sin(ave_body), com_ave(1);
    sin(ave_body)  cos(ave_body), com_ave(2);
    0              0, 1]; % set the initla config to origin
%%
% Generate code to make the plot
figure()
hold on;
axis equal;
axis([-15 15 -15 15]*2*L);
xlabel('X (cm)','fontsize',15);
ylabel('Y (cm)','fontsize',15);


[h, com_ini] = drawActiveFrameSnake(g_h, alpha(t), beta(phase(t)), paras.leg_posi, L, paras.Lleg, 0, activation(phase(t))); %% model switch
while t < period*2*pi
    % update time
    t = t + dt;
    alphas = alpha(t);
    dalphas = d_alpha(t);
    % betas = beta(t);
    % dbetas = d_beta(t);
    c = phase(t);
    betas = beta(c);
    dbetas = d_beta(c);
    activations = activation(c);
    [g, gLeg] = framesInHead(alphas, betas, paras.leg_posi, L, paras.Lleg);
    [J, JLeg] = spatialJacobian(alphas, paras.leg_posi, L);
    
    options = optimset('Display','off');
    % options = optimset('Display','off','MaxIter',40000000,'MaxFunEvals',80000000);
    if isRFT
        xi_vec = computeBalancedBodyVelocity(alphas, dalphas, betas, dbetas, paras.leg_posi, L, paras.Lleg, activations, K, [], 1); %%
        % disp(xi_vec)
        xi = xi_vec(1:3);
    else
        xi = computeBalancedBodyVelocity(alphas, dalphas, betas, dbetas, paras.leg_posi, L, paras.Lleg, activations, K, [0;0;0], 0);
    end


    g_h = g_h*expm(dt*xiHat(xi));
    


  for i = 1:size(h,2)
        delete(h{i});
    end

    [h, com] = drawActiveFrameSnake(g_h, alpha(t), beta(c), paras.leg_posi, L, paras.Lleg, 1, activations);
    
    drawnow;
   
end
ave_body=sum(cumsum(alpha(t)))/N;
gf_ave_body=g_h*[cos(ave_body) -sin(ave_body), com_ave(1);
    sin(ave_body)  cos(ave_body), com_ave(2);
    0              0,          1;];
rot_m_ave_body=inv(gi_ave_body)*gf_ave_body;

% [h,com] = drawActiveFrameSnake(g_h,alpha(w1,w2),L,activation,blackred);
x_displacement=rot_m_ave_body(1,3)/(2*L*(N-1))/basis_paras.frac;
y_displacement=rot_m_ave_body(2,3)/(2*L*(N-1))/basis_paras.frac;
theta_displacement=atan2(rot_m_ave_body(2,1),rot_m_ave_body(1,1))/pi*180;

displacement.x_displacement=x_displacement;
displacement.y_displacement=y_displacement;
displacement.xy = norm([com_ini - com])/(2*L*(N-1));
displacement.theta_displacement=theta_displacement;
close(gcf)
