function act=get_act(c,lfs,num_leg,dutyf,symmode)
% Get leg contact based on spatial frequency and duty factor
if symmode == 1 || symmode == 2
    act=zeros(2,num_leg);
    for ind=1:num_leg
        act(:,ind)=[F_leg_act(c+(ind-1)*lfs*2*pi,dutyf), F_leg_act(c+pi+(ind-1)*lfs*2*pi,dutyf)];% one more pi in the right
    end
    act=act(:)';
elseif symmode==0
    act=zeros(2,num_leg);
    for ind=1:num_leg
        act(:,ind)=[F_leg_act(c+(ind-1)*lfs*2*pi,dutyf), F_leg_act(c+(ind-1)*lfs*2*pi,dutyf)];%Synchronous
    end
    act=act(:)';
    
end
