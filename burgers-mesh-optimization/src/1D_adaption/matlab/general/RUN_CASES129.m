% RUN_CASES129.m
global sol lim_type functional flux_type imax

imax = 129;
III = 55:56; % Cases to be run

%      sol      lim_type     functional  flux_type    
var = {'sub', 'vanalbada','cont','roe'   % 1
       'sub', 'vanalbada','xmtm', 'roe'
       'sub', 'vanalbada','energy', 'roe'
       'sub', 'vanalbada','rss','roe'
       'super', 'vanalbada','cont','roe'
       'super', 'vanalbada','xmtm','roe'
       'super', 'vanalbada','energy','roe'
       'super', 'vanalbada','rss','roe'
       'shock', 'vanalbada','cont','roe'
       'shock', 'vanalbada','xmtm','roe'
       'shock', 'vanalbada','energy','roe'
       'shock', 'vanalbada','rss','roe'   % 12
       'sub', 'vanleer','cont','roe'
       'sub', 'vanleer','xmtm','roe'
       'sub', 'vanleer','energy','roe'
       'sub', 'vanleer','rss','roe'
       'super', 'vanleer','cont','roe'
       'super', 'vanleer','xmtm','roe'
       'super', 'vanleer','energy','roe'
       'super', 'vanleer','rss','roe'
       'shock', 'vanleer','cont','roe'
       'shock', 'vanleer','xmtm','roe'
       'shock', 'vanleer','energy','roe'
       'shock', 'vanleer','rss','roe'     % 24
       'sub', 'off','cont','roe'
       'sub', 'off','xmtm','roe'
       'sub', 'off','energy','roe'
       'sub', 'off','rss','roe'
       'super', 'off','cont','roe'
       'super', 'off','xmtm','roe'
       'super', 'off','energy','roe'
       'super', 'off','rss','roe'
       'shock', 'off','cont','roe'
       'shock', 'off','xmtm','roe'
       'shock', 'off','energy','roe'
       'shock', 'off','rss','roe'    %36
       'sub', 'vanalbada','cont','vanleer'
       'sub', 'vanalbada','xmtm', 'vanleer'
       'sub', 'vanalbada','energy', 'vanleer'
       'sub', 'vanalbada','rss','vanleer'
       'super', 'vanalbada','cont','vanleer'
       'super', 'vanalbada','xmtm','vanleer'
       'super', 'vanalbada','energy','vanleer'
       'super', 'vanalbada','rss','vanleer'
       'shock', 'vanalbada','cont','vanleer'
       'shock', 'vanalbada','xmtm','vanleer'
       'shock', 'vanalbada','energy','vanleer'
       'shock', 'vanalbada','rss','vanleer' %48
       'sub', 'vanleer','cont','vanleer'
       'sub', 'vanleer','xmtm','vanleer'
       'sub', 'vanleer','energy','vanleer'
       'sub', 'vanleer','rss','vanleer'
       'super', 'vanleer','cont','vanleer'
       'super', 'vanleer','xmtm','vanleer'
       'super', 'vanleer','energy','vanleer'
       'super', 'vanleer','rss','vanleer'
       'shock', 'vanleer','cont','vanleer'
       'shock', 'vanleer','xmtm','vanleer'
       'shock', 'vanleer','energy','vanleer'
       'shock', 'vanleer','rss','vanleer' % 60
       'sub', 'off','cont','vanleer'
       'sub', 'off','xmtm','vanleer'
       'sub', 'off','energy','vanleer'
       'sub', 'off','rss','vanleer'
       'super', 'off','cont','vanleer'
       'super', 'off','xmtm','vanleer'
       'super', 'off','energy','vanleer'
       'super', 'off','rss','vanleer'
       'shock', 'off','cont','vanleer'
       'shock', 'off','xmtm','vanleer'
       'shock', 'off','energy','vanleer'
       'shock', 'off','rss','vanleer'};   %72
   
for ii = III(1:end),

    sol = var{ii,1};
    lim_type = var{ii,2};
    functional = var{ii,3};
    flux_type = var{ii,4};
    
    general_adaption_scipt
    fclose all;
    mkdir([strcat('Case',num2str(imax)) int2str(ii) ]);
    
   
     movefile('TE.dat',[strcat('Case',num2str(imax)) int2str(ii) ]);
      movefile('TEuniform.dat',[strcat('Case',num2str(imax)) int2str(ii) ]);
       movefile('config.mat',[strcat('Case',num2str(imax)) int2str(ii) ]);
        movefile('exact.dat',[strcat('Case',num2str(imax)) int2str(ii) ]);
         movefile('faces.dat',[strcat('Case',num2str(imax)) int2str(ii) ]);
          movefile('grid.grd',[strcat('Case',num2str(imax)) int2str(ii) ]);
           movefile('history.dat',[strcat('Case',num2str(imax)) int2str(ii) ]);
            movefile('opt-hist-1d.out',[strcat('Case',num2str(imax)) int2str(ii) ]);
    
end

exit
   