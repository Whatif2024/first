function [REP]= mopso(W)
if nargin < 1
    W = ones(1,24)*0.5; % 默认EVCA购电比例
end
%function [REP]= mopso(W,c,iw,max_iter,lower_bound,upper_bound,swarm_size,rep_size,grid_size,alpha,beta,gamma,mu,problem)
%mopso is an implementation of multi objective particle swarm optimization
%technique for a minimization problem
%% initialize parameters
global P_PV;
global P_WT;
%CHP设备
MinP_CHP = 0; %CHP最小输出功率
MaxP_CHP = 800; %CHP最大输出功率
%n_CHP_e = 0.92; %CHP发电效率
%n_CHP_h = n_CHP_e * 0.75;%热电比
%P_CHP_e + P_CHP_h == P_CHP;
%P_CHP_e = P_CHP * n_CHP_e;
%P_CHP_h = P_CHP * n_CHP_h;
%GB设备
MinP_GB_h = 0; %GB最小输出功率
MaxP_GB_h = 800; %GB最大输出功率
%n_GB = 0.92; %GB发电效率
%P_GB_h = P_g_GB * n_GB;
%HP设备
MinP_HP_h = 0;
MaxP_HP_h = 500;
%n_HP = 4.5;
%P_HP_h = P_e_GB * n_HP;
%AC设备
MinP_AC_c = 0;
MaxP_AC_c = 600;
%n_AC = 1.2;
%P_AC_c = P_h_AC * n_AC;
%EC设备
MinP_EC_c = 0;
MaxP_EC_c = 800;
%n_EC = 6;
%P_EC_c = P_e_EC * n_EC;
%电网最大购电功率
MaxP_Grid = 5000;
%电网最小购电功率
MinP_Grid = -2000;

% 初始化所有参数（移到前面避免变量作用域问题）
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

%%%%%%可以在此处添加需求侧响应%%%%%%
for j=1:168 %粒子长度为264/192/168（光伏，风电，CHP, GB, HP, AC, EC, 电网，冷热电负荷响应的8*24个小时出力）
     if j<=24
        lower_bound(j) = 0;
        upper_bound(j) = P_PV(j);
     elseif j<=48
        lower_bound(j) = 0;
        upper_bound(j) = P_WT(j-24);
     elseif j<=72
        lower_bound(j) = MinP_CHP;
        upper_bound(j) = MaxP_CHP;
     elseif j<=96
        lower_bound(j) = MinP_GB_h;
        upper_bound(j) = MaxP_GB_h;
     elseif j<=120
        lower_bound(j) = MinP_HP_h;
        upper_bound(j) = MaxP_HP_h;
     elseif j<=144
        lower_bound(j) = MinP_AC_c;
        upper_bound(j) = MaxP_AC_c;
     else % j<=168
        lower_bound(j) = MinP_EC_c;
        upper_bound(j) = MaxP_EC_c;
     end
     
%{
     if j>168
        lower_bound(j) = MinP_Grid;
        upper_bound(j) = MaxP_Grid;
     end

     if j>192&&j<217
        lower_bound(j) = -100;
        upper_bound(j) = 100;
     end
     if j>216&&j<241
        lower_bound (j)=-200;
        upper_bound(j)=200;
     end
     if j>240
        lower_bound (j)=-300;
        upper_bound(j)=300;
     end
%}
end

%problem=@prob; %objective function
problem=@(x) prob(x, W); % 绑定W参数

%% initialize particles
fprintf('Initializing swarm ...\n')
w = @(it) ((max_iter - it) - (iw(1) - iw(2)))/max_iter + iw(2);
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

REP = Repository(swarm,rep_size,grid_size,alpha,beta,gamma);
%% Loop
fprintf('Starting the optimization loop ...\n')
for it=1:max_iter
    leader = REP.SelectLeader();
    wc = w(it); %current inertia weight
    pc = pm(it); %current mutation rate
    for i =1:swarm_size %update particles
        swarm(i)=swarm(i).update(wc,c,pc,leader,problem);
    end
    REP = REP.update(swarm);
    Title = sprintf('Iteration %d, Number of Rep Members = %d',it,length(REP.swarm));
    PlotCosts(swarm,REP.swarm,Title)
    disp(Title);
end
end

function PlotCosts(swarm,rep,Title)
figure(1)
feasable_swarm = swarm([swarm.infeasablity]==0);
infeasable_swarm = swarm([swarm.infeasablity]>0);
LEG = {};
if ~isempty(feasable_swarm)
    swarm_costs=vertcat(feasable_swarm.cost);
    %plot3(1000000-swarm_costs(:,1),swarm_costs(:,2), swarm_costs(:,3),'go')
    plot(-swarm_costs(:,1),swarm_costs(:,2),'go')
    hold on
    LEG = {'Current feasable SWARM'};
    Title = sprintf([Title '\nfeasable swarm=%d'],length(feasable_swarm));
end
if ~isempty(infeasable_swarm)
    swarm_costs=vertcat(infeasable_swarm.cost);
    %plot3(1000000-swarm_costs(:,1),swarm_costs(:,2),swarm_costs(:,3),'ro')
    plot(-swarm_costs(:,1),swarm_costs(:,2),'ro')
    hold on
    LEG = [LEG, 'Current infeasable SWARM'];
    if contains(Title,newline)
        Title = sprintf([Title ', infeasable swarm=%d'],length(infeasable_swarm));
    else
        Title = sprintf([Title '\ninfeasable swarm=%d'],length(infeasable_swarm));
    end
end
rep_costs=vertcat(rep.cost);
%plot3(1000000-rep_costs(:,1),rep_costs(:,2), rep_costs(:,3),'b*')
plot(-rep_costs(:,1),rep_costs(:,2),'b*')
xlabel('IESO收益')
ylabel('碳排放成本')
%zlabel('3^{nd} Objective')
grid on
hold off
title(Title)
legend([LEG ,'REPASITORY'],'location','best')
drawnow
end