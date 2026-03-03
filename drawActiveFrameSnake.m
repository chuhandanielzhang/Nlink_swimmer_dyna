function [ h, com ] = drawActiveFrameSnake( g, alpha, beta, leg_posi, L, Lleg, check_stat, activations)
    n = size(alpha,1)+1;
    h = cell(1,n);
    
    body_color = [0.8 0.8 0.8];
    leg_color_contact = [1 0 0];
    leg_color_lift = [0.2 0.2 0.2];

    er=0.2;
    
    if 1==check_stat
        h{1} = drawActiveFrameEllipse(g, L, body_color);
    end
    
    [k, gLeg] = framesInHead(alpha, beta, leg_posi, L, Lleg);
    k{1} = g;
    
    for i = 2:n
        k{i} = g*k{i};
        if 1==check_stat
            h{i} = drawActiveFrameEllipse(k{i}, L, body_color);
        end
    end
    
    com = computeCOM(k);
    
    if 1==check_stat
        er=1;
        for ind=1:length(leg_posi)
            % 获取当前腿部的接触状态
            if nargin >= 8 && ~isempty(activations)
                % 左腿接触状态 (activations中从第n+1个元素开始是腿部状态)
                left_leg_contact = activations(n + ind*2 - 1) > 0;
                right_leg_contact = activations(n + ind*2) > 0;
                
                % 根据接触状态选择颜色
                if left_leg_contact
                    left_leg_color = leg_color_contact;
                else
                    left_leg_color = leg_color_lift;
                end
                
                if right_leg_contact
                    right_leg_color = leg_color_contact;
                else
                    right_leg_color = leg_color_lift;
                end
            else
                % 如果没有提供激活状态，使用默认黑色
                left_leg_color = leg_color_lift;
                right_leg_color = leg_color_lift;
            end
            
            % 使用更大的圆形来绘制腿部
            leg_size = L/2.5;  % 增加腿部大小 (从L/5改为L/2.5)
            h{end+1} = drawActiveFrameCircle(g*gLeg{1,ind}, leg_size, left_leg_color);
            h{end+1} = drawActiveFrameCircle(g*gLeg{2,ind}, leg_size, right_leg_color);
        end
    end
end




