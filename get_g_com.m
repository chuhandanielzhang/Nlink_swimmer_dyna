function g_com = get_g_com(alpha,L)
n=size(alpha,1)+1;

g = cell(1,n);
g{1} = eye(3);% module 0 configuration: identity in itself
x = zeros(1,n);
y = zeros(1,n);
t = zeros(1,n);
F = [eye(2),[L;0];0 0 1];%move forward
R = @(a)[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];%rotation frame
for i=2:n
    g{i} = g{i-1}*F*R(alpha(i-1))*F;
    x(i) = g{i}(1,3);
    y(i) = g{i}(2,3);
    t(i) = atan2(g{i}(2,1),g{i}(1,1));
end

x_bar = mean(x);
y_bar = mean(y);
t_bar = mean(t);

g_com = R(t_bar) + [0 0 x_bar; 0 0 y_bar; 0 0 0];
g_links = framesInHead(alpha, L);

com_pos = computeCOM(g_links);

g_comm = eye(3);
g_comm(1:2, 3) = com_pos;
% disp('...')
% disp(g_com)
% disp('+++')
% disp(g_comm)
% disp('...')



end
% 
% 
% function [ g_com ] = get_g_com( alpha, L )
% 
% g_links = framesInHead(alpha, L);
% 
% com_pos = computeCOM(g_links);
% 
% g_com = eye(3);
% g_com(1:2, 3) = com_pos;
% 
% end