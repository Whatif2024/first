function [cost, infeasablity] = prob(x, W)
% Multi-objective optimization problem for energy system
% x: decision variables [P_PV(1:24), P_WT(1:24), P_CHP(1:24), P_GB(1:24), P_HP(1:24), P_AC(1:24), P_EC(1:24)]
% W: EVCA购电比例 (24小时)
% cost: [负IESO收益, 碳排放成本] (双目标最小化)
% infeasablity: 约束违反程度

global P_PV;
global P_WT;

% 如果全局变量未定义，使用默认值
if isempty(P_PV)
    P_PV = 200 + 100 * sin((1:24) * 2 * pi / 24); % 示例光伏出力曲线
end
if isempty(P_WT)
    P_WT = 150 + 50 * cos((1:24) * 2 * pi / 24); % 示例风电出力曲线
end

% 提取决策变量
P_PV_actual = x(1:24);
P_WT_actual = x(25:48);
P_CHP = x(49:72);
P_GB = x(73:96);
P_HP = x(97:120);
P_AC = x(121:144);
P_EC = x(145:168);

% 系统参数
n_CHP_e = 0.92; % CHP发电效率
n_CHP_h = 0.69; % CHP热效率
n_GB = 0.92; % GB效率
n_HP = 4.5; % HP热泵系数
n_AC = 1.2; % 吸收式制冷机系数
n_EC = 6; % 电制冷系数

% 计算各设备输出
P_CHP_e = P_CHP * n_CHP_e; % CHP电输出
P_CHP_h = P_CHP * n_CHP_h; % CHP热输出
P_GB_h = P_GB * n_GB; % GB热输出
P_HP_h = P_HP * n_HP; % HP热输出
P_AC_c = P_AC * n_AC; % AC冷输出
P_EC_c = P_EC * n_EC; % EC冷输出

% 假设负荷需求（示例数据）
P_load_e = 800 + 200 * sin((1:24) * 2 * pi / 24 + pi/4); % 电负荷
P_load_h = 600 + 150 * cos((1:24) * 2 * pi / 24 + pi/6); % 热负荷
P_load_c = 400 + 100 * sin((1:24) * 2 * pi / 24 + pi/3); % 冷负荷

% 电力平衡
P_grid = P_load_e - P_PV_actual - P_WT_actual - P_CHP_e + P_HP + P_EC; % 电网购电

% 约束违反检查
infeasablity = 0;

% 电网购电约束
grid_violation = max(0, max(P_grid - 5000)) + max(0, max(-2000 - P_grid));
infeasablity = infeasablity + grid_violation;

% 热平衡约束
heat_balance_violation = sum(max(0, P_load_h - P_CHP_h - P_GB_h - P_HP_h));
infeasablity = infeasablity + heat_balance_violation;

% 冷平衡约束
cool_balance_violation = sum(max(0, P_load_c - P_AC_c - P_EC_c));
infeasablity = infeasablity + cool_balance_violation;

% 可再生能源约束
renewable_violation = sum(max(0, P_PV_actual - P_PV)) + sum(max(0, P_WT_actual - P_WT));
infeasablity = infeasablity + renewable_violation;

% 目标函数计算
% 目标1: 负IESO收益（最大化收益 = 最小化负收益）
electricity_price = 0.1 + 0.05 * sin((1:24) * 2 * pi / 24); % 电价（元/kWh）
gas_price = 0.08; % 天然气价格（元/kWh）
maintenance_cost = 0.02; % 维护成本（元/kWh）

% 收益计算
revenue = sum(P_grid .* electricity_price .* W); % 售电收益
cost_fuel = sum(P_CHP * gas_price + P_GB * gas_price); % 燃料成本
cost_maintenance = sum((P_CHP + P_GB + P_HP + P_AC + P_EC) * maintenance_cost); % 维护成本

IESO_profit = revenue - cost_fuel - cost_maintenance;
objective1 = -IESO_profit; % 负收益（最小化）

% 目标2: 碳排放成本
carbon_factor_grid = 0.8; % 电网碳排放因子（kg CO2/kWh）
carbon_factor_gas = 0.2; % 天然气碳排放因子（kg CO2/kWh）
carbon_price = 0.05; % 碳价（元/kg CO2）

carbon_emissions = sum(P_grid * carbon_factor_grid + P_CHP * carbon_factor_gas + P_GB * carbon_factor_gas);
objective2 = carbon_emissions * carbon_price; % 碳排放成本

% 返回目标函数值
cost = [objective1, objective2];

end