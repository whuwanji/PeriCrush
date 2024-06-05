function model = readModel(fileDir, partName)
% 根据part的名字以及结果所在文件夹读入模型信息
model.ndim = 3;
tb = readtable([fileDir, partName, 'GeometryMaterialConstant.dat']);
model.dx          = tb{strcmp(tb{:,1},'delta_x'),2};
model.delta       = tb{strcmp(tb{:,1},'delta'),2};
model.length      = tb{strcmp(tb{:,1},'length'),2};
model.width       = tb{strcmp(tb{:,1},'width'),2};
model.height      = tb{strcmp(tb{:,1},'height'),2};
model.xn          = tb{strcmp(tb{:,1},'xn'),2};
model.yn          = tb{strcmp(tb{:,1},'yn'),2};
model.zn          = tb{strcmp(tb{:,1},'zn'),2};
model.dy          = tb{strcmp(tb{:,1},'dy'),2};
model.dz          = tb{strcmp(tb{:,1},'dz'),2};
model.mvalue      = tb{strcmp(tb{:,1},'mvalue'),2};
model.pn          = tb{strcmp(tb{:,1},'particleNumber'),2};
model.hpn         = tb{strcmp(tb{:,1},'horizonParticleNumber'),2};
model.emod        = tb{strcmp(tb{:,1},'YoungsModulus'),2};
model.LameConstant = tb{strcmp(tb{:,1},'LameConstant'),2};
model.smod         = tb{strcmp(tb{:,1},'ShearModulus'),2};
model.bmod         = tb{strcmp(tb{:,1},'BulkModulus'),2};
model.pr           =     tb{strcmp(tb{:,1},'PoissonsRatio'),2};
model.spc          = tb{strcmp(tb{:,1},'SpringConstant'),2};
model.dens         =  tb{strcmp(tb{:,1},'MassDensity'),2};
model.scr          =  tb{strcmp(tb{:,1},'CriticalStretch'),2};
model.energyReleaseRate =  tb{strcmp(tb{:,1},'energyReleaseRate'),2};
model.dt                =  tb{strcmp(tb{:,1},'dt'),2};
% 读入Coordinate
fid = fopen([fileDir,partName,'Coordinate.bin'],'rb');
[Coordinate,~] = fread(fid,'float64');
model.ndim = numel(Coordinate)/model.pn;
model.Coordinate = reshape(Coordinate, model.pn, model.ndim);
fclose(fid);
% 读入 RepresentativeRadius
fid = fopen([fileDir,partName,'RepresentativeRadius.bin'],'rb');
[RepresentativeRadius,~] = fread(fid,'float64');
model.RepresentativeRadius = reshape(RepresentativeRadius, model.pn, 1);
fclose(fid);
% 读入 RepresentativeSideLength
fid = fopen([fileDir,partName,'RepresentativeSideLength.bin'],'rb');
[RepresentativeSideLength,~] = fread(fid,'float64');
model.RepresentativeSideLength = reshape(RepresentativeSideLength, model.pn, 1);
fclose(fid);
% 读入 ParticleVolume
fid = fopen([fileDir,partName,'ParticleVolume.bin'],'rb');
[ParticleVolume,~] = fread(fid,'float64');
model.ParticleVolume = reshape(ParticleVolume, model.pn, 1);
fclose(fid);
% 读入 Horizon
fid = fopen([fileDir,partName,'Horizon.bin'],'rb');
[Horizon,~] = fread(fid,'int32');
model.Horizon = Horizon(:);
fclose(fid);
% 读入 HorizonParticleStartNumber
fid = fopen([fileDir,partName,'HorizonParticleStartNumber.bin'],'rb');
[HorizonParticleStartNumber,~] = fread(fid,'int32');
model.HorizonParticleStartNumber = HorizonParticleStartNumber(:);
fclose(fid);
% 读入 HorizonParticleEndNumber
fid = fopen([fileDir,partName,'HorizonParticleEndNumber.bin'],'rb');
[HorizonParticleEndNumber,~] = fread(fid,'int32');
model.HorizonParticleEndNumber = HorizonParticleEndNumber(:);
fclose(fid);
%读入 HorizonParticleNumber
fid = fopen([fileDir,partName,'HorizonParticleNumber.bin'],'rb');
[HorizonParticleNumber,~] = fread(fid, 'int32');
model.HorizonParticleNumber = HorizonParticleNumber(:);
fclose(fid);
model.partName = partName;
model.fileDir  = fileDir;
end