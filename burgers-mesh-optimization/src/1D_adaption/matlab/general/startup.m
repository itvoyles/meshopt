function startup
% Script for starting examination of the mesh opt problem

% Set path

working_dir = '/home/grad4/itvoyles/Meshopt/Repository/burgers-mesh-optimization/';
% working_dir = 'Z:/Meshopt/Repository/burgers-mesh-optimization/';
exe_path = [working_dir,'bin'];
path(path,exe_path);
exe_path = [working_dir,'src/1D_adaption/matlab/burgers'];
path(path,exe_path);
exe_path = [working_dir,'src/1D_adaption/matlab/q1d'];
path(path,exe_path);
exe_path = [working_dir,'src/1D_adaption/matlab/spring-f'];
path(path,exe_path);
exe_path = [working_dir,'src/1D_adaption/matlab/spring-k'];
path(path,exe_path);

% Open functions
edit runner.m
edit general_adaption_scipt.m
edit q1d_TE.m
edit exactsol.m
edit super.m
edit gridmake

end
