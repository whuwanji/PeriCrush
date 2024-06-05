clc; clear
% 本程序可以展示任意时刻的变量散点图
fileDir    = '..\res\exam01\';                                   % 输出结果的文件夹
partName   = 'S1';                                               % 模型中哪个Part的名字
model      = readModel(fileDir, partName);                       % 读入Part的信息     
stepNumber = 160000;                                             % 第几步
dmg        = readStepVariable(model, stepNumber, 'damage');      % 读入第stepNumber步的损伤
dis        = readStepVariable(model, stepNumber, 'displacement');% 读入第stepNumber步的位移

figure(2); clf  % 画损伤图
scatters(model.Coordinate*1e3, dmg(:,1), 15, 20); caxis([0,1])
xlabel('X(mm)'), ylabel('Y(mm)'),zlabel('Z(mm)');
title(['damage:  time=', num2str(stepNumber*model.dt*1e6), '\mus'])
set(gca, 'fontsize', 16, 'fontname', 'times new roman')

figure(3); clf  % 画位移图
scatters(model.Coordinate*1e3, dis(:,3), 15, 20);
xlabel('X(mm)'), ylabel('Y(mm)'),zlabel('Z(mm)');
title(['u_3: time=', num2str(stepNumber*model.dt*1e6), '\mus'])
set(gca, 'fontsize', 16, 'fontname', 'times new roman')