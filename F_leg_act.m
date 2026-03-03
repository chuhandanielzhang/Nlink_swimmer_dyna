function leg_act_col=F_leg_act(time,dutyf)

leg_act_col=zeros(size(time));
time = mod(time, 2*pi);
leg_act_col(time < 2*pi-(1-dutyf)*pi) = 2;
leg_act_col(time < (1-dutyf)*pi) = 0;
