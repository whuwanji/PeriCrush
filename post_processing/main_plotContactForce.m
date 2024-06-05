clc; clear
% 本程序画接触力图
fileDir    = '..\res\exam01\';                                   % 输出结果的文件夹
partName   = 'S1';                                               % 模型中哪个Part的名字
model      = readModel(fileDir, partName);                       % 读入Part的信息     
cf = load([fileDir, 'contactForce.txt']);
cf3 = cf(:,3);
u3 = (1:size(cf,1))'*100*model.dt*0.1*1e6;
plot(u3, -cf3)
xlabel('U_3(\mum)')
ylabel('Contact Force(N)')
set(gca, 'fontsize', 16, 'fontname', 'times new roman')