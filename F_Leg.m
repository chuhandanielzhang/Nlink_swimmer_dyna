function leg_act_col=F_Leg(Aleg,time,dutyf)
% leg kinematics
time = mod(time, 2*pi);
leg_act_col = -Aleg*cos((time-2*pi+(1-dutyf)*pi)/2/(1-dutyf));%swing cont
leg_act_col(time < 2*pi-(1-dutyf)*pi) = Aleg*cos((time-(1-dutyf)*pi)*(pi/(2*pi-2*(1-dutyf)*pi)));%stance
leg_act_col(time < (1-dutyf)*pi) = Aleg*sin(time/2/(1-dutyf));