function burgers2d_solver(nodes,reynoldsnumber)

%setup executable path
exe_path = '~/Software/burgers-mesh-adaption/library/bin';
path(path,exe_path);

exe = ['burgers2d_n',num2str(nodes),'_re',num2str(reynoldsnumber),'.exe'];

system(exe);

end