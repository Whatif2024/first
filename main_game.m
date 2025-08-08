function main_game()
% 主博弈迭代函数
% 实现EVCA和IESO之间的博弈过程

% 初始化全局变量
global P_PV;
global P_WT;

% 设置可再生能源出力数据（24小时）
P_PV = 200 + 100 * sin((1:24) * 2 * pi / 24); % 光伏出力曲线
P_WT = 150 + 50 * cos((1:24) * 2 * pi / 24); % 风电出力曲线

% 初始化EVCA购电比例
W = ones(1,24) * 0.5; % 初始购电比例50%

% 博弈参数
max_iterations = 10; % 最大博弈迭代次数
tolerance = 1e-3; % 收敛容差
W_history = zeros(max_iterations, 24); % 存储历史W值

fprintf('开始博弈迭代过程...\n');

for iter = 1:max_iterations
    fprintf('\n=== 博弈迭代 %d ===\n', iter);
    
    % 存储当前W值
    W_history(iter, :) = W;
    
    % IESO优化：给定W，求解MOPSO
    try
        mm = mopso(W);  % 传入当前W值
        
        % 检查MOPSO结果
        if isempty(mm.swarm)
            fprintf('警告: MOPSO未找到可行解\n');
            break;
        end
        
        % 选择一个代表性解（这里选择第一个非支配解）
        best_solution = mm.swarm(1);
        fprintf('IESO优化完成，找到 %d 个Pareto最优解\n', length(mm.swarm));
        fprintf('选择解的目标值: [%.4f, %.4f]\n', best_solution.cost(1), best_solution.cost(2));
        
        % EVCA响应：根据IESO的策略更新W
        % 这里实现一个简单的响应策略
        W_new = updateEVCAStrategy(W, best_solution, iter);
        
        % 检查收敛性
        convergence = norm(W_new - W);
        fprintf('收敛指标: %.6f\n', convergence);
        
        if convergence < tolerance
            fprintf('博弈达到收敛!\n');
            break;
        end
        
        % 更新W值
        W = W_new;
        
    catch ME
        fprintf('错误: %s\n', ME.message);
        break;
    end
end

% 输出最终结果
fprintf('\n=== 博弈结果 ===\n');
fprintf('最终EVCA购电比例:\n');
for h = 1:24
    fprintf('第%2d小时: %.4f\n', h, W(h));
end

% 绘制收敛过程
if iter > 1
    figure(2);
    plot(1:iter, W_history(1:iter, 1:6), 'LineWidth', 2);
    xlabel('博弈迭代次数');
    ylabel('EVCA购电比例');
    title('博弈收敛过程（前6小时）');
    legend('1时', '2时', '3时', '4时', '5时', '6时', 'Location', 'best');
    grid on;
end

end

function W_new = updateEVCAStrategy(W_current, ieso_solution, iteration)
% EVCA策略更新函数
% W_current: 当前购电比例
% ieso_solution: IESO的最优解
% iteration: 当前迭代次数

% 提取IESO解决方案中的电网购电功率
x = ieso_solution.position;
P_PV_actual = x(1:24);
P_WT_actual = x(25:48);
P_CHP = x(49:72);

% 计算电网购电需求
global P_PV;
global P_WT;
n_CHP_e = 0.92;
P_CHP_e = P_CHP * n_CHP_e;
P_load_e = 800 + 200 * sin((1:24) * 2 * pi / 24 + pi/4); % 电负荷
P_grid = P_load_e - P_PV_actual - P_WT_actual - P_CHP_e;

% EVCA响应策略：基于电网购电量和电价调整购电比例
electricity_price = 0.1 + 0.05 * sin((1:24) * 2 * pi / 24);

% 学习率（随迭代次数递减）
learning_rate = 0.3 / iteration;

% 策略调整
W_adjustment = zeros(1, 24);
for h = 1:24
    % 如果电价高且需要从电网购电，减少购电比例
    if P_grid(h) > 0 && electricity_price(h) > 0.12
        W_adjustment(h) = -learning_rate * 0.1;
    % 如果电价低且向电网售电，增加购电比例
    elseif P_grid(h) < 0 && electricity_price(h) < 0.08
        W_adjustment(h) = learning_rate * 0.1;
    end
end

% 更新购电比例
W_new = W_current + W_adjustment;

% 约束在[0,1]范围内
W_new = max(0, min(1, W_new));

fprintf('EVCA策略更新完成，平均调整幅度: %.6f\n', mean(abs(W_adjustment)));

end