% ==============================================
% Gurobi求解线性规划（LP）示例（修正sense格式）
% ==============================================
clear; clc; close all;

% 1. 初始化Gurobi模型
model = struct();

% 2. 设置目标函数：min f = 2x1 + 3x2 + x3
model.obj = [2, 3, 1];  % 目标函数系数（对应x1,x2,x3）
model.modelsense = 'min';  % 优化方向：'min'（最小化）/ 'max'（最大化）

% 3. 设置约束条件
% 约束矩阵A：每一行对应一个约束，列对应变量x1,x2,x3
model.A = sparse([1, 1, 1; 2, 1, 0]);  % 稀疏矩阵（2行3列，2个约束）
model.rhs = [10, 8];  % 约束右侧值（2个元素，对应2个约束）
% 关键修正：用单元格数组定义约束类型（长度=约束数=2）
model.sense = {'>=', '>='};  % 替代原错误的['>=', '>=']

% 4. 设置变量上下界（默认下界0，上界无穷，可省略）
model.lb = [0, 0, 0];  % 变量下界 x1,x2,x3 ≥ 0
model.ub = [inf, inf, inf];  % 变量上界（无穷大）

% 5. （可选）设置求解参数
params = struct();
params.OutputFlag = 1;  % 1=显示求解日志，0=不显示
params.TimeLimit = 300;  % 最大求解时间（秒）

% 6. 调用Gurobi求解
result = gurobi(model, params);

% 7. 输出结果
fprintf('==================== 求解结果 ====================\n');
if strcmp(result.status, 'OPTIMAL')  % 最优解存在
    fprintf('最优目标函数值：%.4f\n', result.objval);
    fprintf('最优变量值：\n');
    fprintf('x1 = %.4f\n', result.x(1));
    fprintf('x2 = %.4f\n', result.x(2));
    fprintf('x3 = %.4f\n', result.x(3));
    
    % （可选）输出约束对偶值
    fprintf('\n约束对偶值（影子价格）：\n');
    fprintf('约束1（x1+x2+x3≥10）：%.4f\n', result.pi(1));
    fprintf('约束2（2x1+x2≥8）：%.4f\n', result.pi(2));
else
    fprintf('求解状态：%s\n', result.status);
    fprintf('无最优解，请检查模型！\n');
end