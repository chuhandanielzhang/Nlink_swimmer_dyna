function [ h ] = drawActiveFrameCircle( g, radius, color )
% 绘制圆形的腿部
% g: 变换矩阵
% radius: 圆的半径
% color: 颜色

% 创建圆形的点
theta = linspace(0, 2*pi, 101);
circle_x = radius * cos(theta);
circle_y = radius * sin(theta);
circle_p = [circle_x; circle_y];

% 应用旋转变换
R = g(1:2, 1:2);
circle_p = R * circle_p;

% 绘制圆形
h = patch(g(1,3) + circle_p(1,:), g(2,3) + circle_p(2,:), color, 'linewidth', 1);

end
