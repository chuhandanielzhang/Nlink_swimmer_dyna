function [ J1, Jleg ] = spatialJacobian( alpha, leg_posi, L )
%alpha are joint angles
%spatial Jacobian computes the spatial manipulator Jacobian with the head
%module defined as the spatial reference frame
%alpha: n X 1, J: 3 X n
F = [eye(2),[L;0];0 0 1];
%compute joint frame position: note the first frame is the head not the
%joint
n = size(alpha,1)+1;
q = zeros(2,n-1);%joint positions
g_joint = jointsInHead(alpha,L);
for i = 1:n-1
    q(:,i) = g_joint{i+1}(1:2,3);
end
%get all the joint positions, and we now the rotational axis is [0 0 1]
%now, we construct the spatial jacobian
q = [q(2,:);-q(1,:)];%position part in the spatial Jacobian
J1 = [q;ones(1,n-1)];

%front legs
Jleg=cell(1,length(leg_posi));
for ind=1:length(leg_posi)
    q0 = g_joint{leg_posi(ind)}*F;
    q0 = q0(1:2, 3);
    q0 = [q0(2);-q0(1)];
    Jleg{ind} = [q0;1];
end

end

