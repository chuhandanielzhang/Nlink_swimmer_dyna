function [ xi ] = computeBodyVelocity( alpha, d_alpha, beta, d_beta,leg_posi, L, Lleg, activation, K, F)
%compute the body velocity of the head frame
%compute all the frames relative to the head frame
if sum(activation)==0
    xi=[0;0;0];
else
    n=size(alpha,1)+1;
    [g, gLeg]=framesInHead(alpha, beta, leg_posi, L, Lleg);%get all the modules frame relative to head
    % k=10;
    % K=diag([k k 200*k]);
    %compute omega1
    %activation: body, front left, front right, back left, back right
    omega1 = activation(1) * (myAdjoint(inv(g{1}))).'*K*myAdjoint(inv(g{1}));
    for i=2:n
        omega1 = omega1+activation(i)*(myAdjoint(inv(g{i}))).'*K*myAdjoint(inv(g{i}));
    end
    for ind=1:length(leg_posi)
        omega1 = omega1 + activation(n+ind*2-1) * (myAdjoint(inv(gLeg{1,ind}))).'*K*myAdjoint(inv(gLeg{1,ind}));
        omega1 = omega1 + activation(n+ind*2) * (myAdjoint(inv(gLeg{2,ind}))).'*K*myAdjoint(inv(gLeg{2,ind}));
    end
    %compute omega2
    %compute spatial Jacobian, front legs, hind legs

    [J1, JLeg] = spatialJacobian(alpha, leg_posi, L);
    omega2=zeros(3,1);
    for i=1:n
        omega2=omega2+activation(i)*(myAdjoint(inv(g{i}))).'*K*myAdjoint(inv(g{i}))*[J1(:,1:i-1),zeros(3,n-i)]*d_alpha;
    end
    for ind=1:length(leg_posi)
        if leg_posi(ind)==1
            omega2=omega2+activation(n+ind*2-1)*(myAdjoint(inv(gLeg{1,ind}))).'*K*myAdjoint(inv(gLeg{1,ind}))*JLeg{ind}*d_beta(ind*2-1);
            omega2=omega2+activation(n+ind*2)*(myAdjoint(inv(gLeg{2,ind}))).'*K*myAdjoint(inv(gLeg{2,ind}))*JLeg{ind}*d_beta(ind*2);
        else
            omega2=omega2+activation(n+ind*2-1)*(myAdjoint(inv(gLeg{1,ind}))).'*K*myAdjoint(inv(gLeg{1,ind}))*[J1(:,1:leg_posi(ind)-1),zeros(3,n-leg_posi(ind)),JLeg{ind}]*[d_alpha;d_beta(ind*2-1)];
            omega2=omega2+activation(n+ind*2)*(myAdjoint(inv(gLeg{2,ind}))).'*K*myAdjoint(inv(gLeg{2,ind}))*  [J1(:,1:leg_posi(ind)-1),zeros(3,n-leg_posi(ind)),JLeg{ind}]*[d_alpha;d_beta(ind*2)];
        end
    end

    omega2 = omega2 ;
   
    xi = -omega1\omega2;
end
end


