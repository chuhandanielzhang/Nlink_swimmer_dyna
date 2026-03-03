function [ wrench ] = forceResidual_1( xi, alpha, d_alpha, beta, d_beta, leg_posi, L, Lleg, activation,force_model)

% xi=[v_x v_y v_theta] with x being the translation direction
% wrench=[F_x F_y F_theta] with x being the translation direction
% xi_t body velocity of each segment
N = size(alpha,1)+1;
wrench = zeros(3,1);
[g, gLeg]=framesInHead(alpha, beta, leg_posi, L, Lleg);
%map body velocity: myAdjoint(inv(g{4}))*J{4}*d_alpha(2)
[J, JLeg] = spatialJacobian(alpha, leg_posi, L);
for l = 1:N
    %for each of the links
    xi_t = myAdjoint(inv(g{l}))*xi+myAdjoint(inv(g{l}))*[J(:,1:l-1),zeros(3,N-l)]*d_alpha;
    wrench = wrench + activation(l)*myAdjoint(inv(g{l}))'*RFT_element(L,xi_t,force_model);    
end

wrench1=zeros(3,1);
for ind=1:length(leg_posi)
    if leg_posi(ind)==1
        xi_t=myAdjoint(inv(gLeg{1,ind}))*xi+myAdjoint(inv(gLeg{1,ind}))*JLeg{ind}*d_beta(ind*2-1);
        wrench1 = wrench1 + activation(N+ind*2-1)*myAdjoint(inv(gLeg{1,ind}))'*RFT_element(L/5,xi_t,force_model); 

        xi_t=myAdjoint(inv(gLeg{2,ind}))*xi+myAdjoint(inv(gLeg{2,ind}))*JLeg{ind}*d_beta(ind*2  );
        wrench1 = wrench1 + activation(N+ind*2)  *myAdjoint(inv(gLeg{2,ind}))'*RFT_element(L/5,xi_t,force_model); 
    else
        xi_t=myAdjoint(inv(gLeg{1,ind}))*xi+myAdjoint(inv(gLeg{1,ind}))*[J(:,1:leg_posi(ind)-1),zeros(3,N-leg_posi(ind)),JLeg{ind}]*[d_alpha;d_beta(ind*2-1)];
        wrench1 = wrench1 + activation(N+ind*2-1)*myAdjoint(inv(gLeg{1,ind}))'*RFT_element(L/5,xi_t,force_model); 

        xi_t=myAdjoint(inv(gLeg{2,ind}))*xi+myAdjoint(inv(gLeg{2,ind}))*[J(:,1:leg_posi(ind)-1),zeros(3,N-leg_posi(ind)),JLeg{ind}]*[d_alpha;d_beta(ind*2)];
        wrench1 = wrench1 + activation(N+ind*2)*myAdjoint(inv(gLeg{2,ind}))'*RFT_element(L/5,xi_t,force_model); 
    end
end
wrench=wrench+wrench1;

end

