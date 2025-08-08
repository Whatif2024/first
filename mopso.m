function [REP]= mopso(W)
if nargin < 1
    W = ones(1,24)*0.5; % 默认EVCA购电比例
end

%% initialize parameters
global P_PV;
global P_WT;

%CHP设备
MinP_CHP = 0; %CHP最小输出功率
MaxP_CHP = 800; %CHP最大输出功率

%GB设备
MinP_GB_h = 0; %GB最小输出功率
MaxP_GB_h = 800; %GB最大输出功率

%HP设备
MinP_HP_h = 0;
MaxP_HP_h = 500;

%AC设备
MinP_AC_c = 0;
MaxP_AC_c = 600;

%EC设备
MinP_EC_c = 0;
MaxP_EC_c = 800;

%电网最大购电功率
MaxP_Grid = 5000;
%电网最小购电功率
MinP_Grid = -2000;

% 默认参数设置 - 移到函数开始部分，确保变量在整个函数中可用
c = [0.1,0.2]; % [cognitive acceleration, social acceleration] coefficients
iw = [0.5 0.001]; % [starting, ending] inertia weight
max_iter = 100; % maximum iterations
swarm_size = 100; %swarm size
rep_size = 100; %Repository Size
grid_size = 7; %Number of Grids per Dimension
alpha = 0.1; %Inflation Rate
beta = 2; %Leader Selection Pressure
gamma = 2; %Deletion Selection Pressure
mu = 0.1; %Mutation Rate

% 初始化边界数组
lower_bound = zeros(1, 168);
upper_bound = zeros(1, 168);

%%%%%%可以在此处添加需求侧响应%%%%%%
for j=1:168 %粒子长度为168（光伏，风电，CHP, GB, HP, AC, EC的7*24个小时出力）
     if j<=24  % 光伏 1-24
        lower_bound(j) = 0;
        upper_bound(j) = P_PV(j);
     elseif j<=48  % 风电 25-48
        lower_bound(j) = 0;
        upper_bound(j) = P_WT(j-24);
     elseif j<=72  % CHP 49-72
        lower_bound(j) = MinP_CHP;
        upper_bound(j) = MaxP_CHP;
     elseif j<=96  % GB 73-96
        lower_bound(j) = MinP_GB_h;
        upper_bound(j) = MaxP_GB_h;
     elseif j<=120  % HP 97-120
        lower_bound(j) = MinP_HP_h;
        upper_bound(j) = MaxP_HP_h;
     elseif j<=144  % AC 121-144
        lower_bound(j) = MinP_AC_c;
        upper_bound(j) = MaxP_AC_c;
     else  % EC 145-168
        lower_bound(j) = MinP_EC_c;
        upper_bound(j) = MaxP_EC_c;
     end
end

% 目标函数定义
problem = @(x) prob(x, W); % 绑定W参数

%% initialize particles
fprintf('Initializing swarm ...\n')
w = @(it) ((max_iter - it) * (iw(1) - iw(2)))/max_iter + iw(2);
pm = @(it) (1-(it-1)/(max_iter-1))^(1/mu);

% 初始化粒子群
swarm = []; % 初始化为空数组
for i = 1:swarm_size
    particle = Particle(lower_bound, upper_bound, problem);
    retry = 0;
    while particle.infeasablity > 0 && retry < 100
        particle = Particle(lower_bound, upper_bound, problem);
        retry = retry + 1;
    end
    swarm = [swarm, particle];
end

REP = Repository(swarm, rep_size, grid_size, alpha, beta, gamma);

%% Loop
fprintf('Starting the optimization loop ...\n')
for it=1:max_iter
    leader = REP.SelectLeader();
    wc = w(it); %current inertia weight
    pc = pm(it); %current mutation rate
    for i = 1:swarm_size %update particles
        swarm(i) = swarm(i).update(wc, c, pc, leader, problem);
    end
    REP = REP.update(swarm);
    Title = sprintf('Iteration %d, Number of Rep Members = %d', it, length(REP.swarm));
    PlotCosts(swarm, REP.swarm, Title)
    disp(Title);
end
end

function PlotCosts(swarm, rep, Title)
figure(1)
feasable_swarm = swarm([swarm.infeasablity]==0);
infeasable_swarm = swarm([swarm.infeasablity]>0);
LEG = {};

if ~isempty(feasable_swarm)
    swarm_costs = vertcat(feasable_swarm.cost);
    plot(-swarm_costs(:,1), swarm_costs(:,2), 'go')
    hold on
    LEG = {'Current feasable SWARM'};
    Title = sprintf([Title '\nfeasable swarm=%d'], length(feasable_swarm));
end

if ~isempty(infeasable_swarm)
    swarm_costs = vertcat(infeasable_swarm.cost);
    plot(-swarm_costs(:,1), swarm_costs(:,2), 'ro')
    hold on
    LEG = [LEG, 'Current infeasable SWARM'];
    if contains(Title, newline)
        Title = sprintf([Title ', infeasable swarm=%d'], length(infeasable_swarm));
    else
        Title = sprintf([Title '\ninfeasable swarm=%d'], length(infeasable_swarm));
    end
end

rep_costs = vertcat(rep.cost);
plot(-rep_costs(:,1), rep_costs(:,2), 'b*')
xlabel('IESO收益')
ylabel('碳排放成本')
grid on
hold off
title(Title)
legend([LEG ,'REPASITORY'], 'location', 'best')
drawnow
end