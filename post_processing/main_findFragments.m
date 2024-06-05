clc; clear
% 本程序可以展示任意时刻的变量散点图
fileDir    = '..\res\exam01\';                                   % 输出结果的文件夹
partName   = 'S1';                                               % 模型中哪个Part的名字
model      = readModel(fileDir, partName);                       % 读入Part的信息     
stepNumber = 168000;                                              % 第几步
fail       = readStepVariable(model, stepNumber, 'fail');        % 读入第stepNumber步的损伤
sr = zeros(model.pn+1,1);
sr(1) = 1;
for i = 1:1:model.pn
    sr(i+1) = sr(i) + model.HorizonParticleNumber(i);
end
frag1 = findFragments( fail, sr, model.Horizon );                 % 计算Fragments
[frag2, fragVol] = volumeDistribution(model, stepNumber, 1.5);    % 计算Fragments以及fragMents对应的体积