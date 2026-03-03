function beta=get_beta(c,num_leg,Aleg,dutyf,lfs,symmode)
%% swing angle
if symmode==1 %c+pi 
    beta=zeros(2,num_leg);
    for ind=1:num_leg
        beta(:,ind)=[F_Leg(Aleg,c+(ind-1)*lfs*2*pi,dutyf), -F_Leg(Aleg,c+pi+(ind-1)*lfs*2*pi,dutyf)];
    end
    beta=beta(:);
elseif symmode==0
    beta=zeros(2,num_leg);
    for ind=1:num_leg
        beta(:,ind)=[F_Leg(Aleg,c+(ind-1)*lfs*2*pi,dutyf), -F_Leg(Aleg,c+(ind-1)*lfs*2*pi,dutyf)];
    end
    beta=beta(:);   
elseif symmode==2
    beta=zeros(2,num_leg);
    for ind=1:num_leg
        beta(:,ind)=[F_Leg(Aleg,c+(ind-1)*lfs*2*pi,dutyf), F_Leg(Aleg,c+pi+(ind-1)*lfs*2*pi,dutyf)];
    end
    beta=beta(:);  
end