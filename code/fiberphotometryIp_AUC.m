clear all; clc;
cd('E:\1.Scn2a_mDAN\fiberphotometry\Scn2a-HE_SCH23390\AUCanalysis\AUC');
cwd = pwd;

fidx = dir('20241224-292R_SCH_Cal580.mat');
load(fidx.name);

%% 设置参数
fs = 40; % 采样频率
inject_point = 1260; % 腹腔注射时间点，单位秒
%offset = 7.744337806
%% 处理 NaN 值
nan_num = sum(isnan(data(:, 3)));
if isnan(data(1, 3))
    data(:, 3) = circshift(data(:, 3), -nan_num);
end

%% 设置 AUC 计算时间范围
timepoint = fs * inject_point;
fit_start = fs * (inject_point + 800) + 1;
fit_end = fs * (inject_point + 1000);
AUC_start = fs * (inject_point - 300) + 1;
AUC_end = timepoint;

a = detrend(data(fit_start:fit_end, 3));
curve = data(fit_start:fit_end, 3) - a; % 拟合出来的一段衰减曲线
k = (curve(end) - curve(1)) / length(curve);
b = curve(1) - k * fit_start;
%b = offset;
t = data(:, 1);
curve_all = k * t * fs + b;
%F0 = mean(curve_all(fit_start:fit_end));
%data_dff = (data(:, 3) - curve_all) ./ F0;
%data_dff = (data(:, 3) - curve_all) ./ mean(curve_all(fit_start:fit_end));
data_dff = (data(:, 3) - curve_all) ./ curve_all;
AUC = sum(data_dff(AUC_start:AUC_end) * (1 / fs));

%% 绘制信号和拟合曲线
plot_range = AUC_start-12000:fit_end;
%plot_range = AUC_start:fit_end;
% 绘制原始信号和拟合曲线
figure(1), clf;
hold on;
plot(data(plot_range, 1), data(plot_range, 3), 'b', 'DisplayName', '原始信号');  % 蓝色线
plot(data(plot_range, 1), curve_all(plot_range), 'r--', 'DisplayName', '拟合曲线', 'LineWidth', 3);  % 红色虚线
xlabel('时间 (s)');
ylabel('信号强度');
title('拟合前后信号及拟合曲线');
legend('show');
grid on;

% 绘制去趋势后的信号
figure(2), clf;
plot(data(plot_range, 1), data_dff(plot_range));
xlabel('时间 (s)');
ylabel('dF/F');
title('去趋势后的信号');
grid on;

%% 写入 Excel
% 提取文件名（去除路径）
[~, filename, ~] = fileparts(fidx.name);

% 创建一个包含表头和数据的cell数组
header = {'文件名', 'AUC'};  % 表头
results = {filename, AUC};   % 数据

% 写入表头到Excel
writecell(header, '0Results.xlsx', 'Sheet', 'AUC', 'Range', 'A1');

% 将文件名和AUC写入Excel文件，从第二行开始
writecell(results, '0Results.xlsx', 'Sheet', 'AUC', 'Range', 'A2');
