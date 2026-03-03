function [ wrench ] = RFT_element( L, xi ,force_model)
% L = 3*L;
ds = L;
wrench = zeros(3,1);
for l = -L/2+ds/2:ds:L/2-ds/2
    g_t = [eye(2),[0;l];0 0 1];%location of the infestimal segment
    xi_t = myAdjoint(inv(g_t))*xi;
    wrench = wrench+myAdjoint(inv(g_t))'*RFT(xi_t,force_model)*ds;
end

end

