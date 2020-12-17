% script that prints the results of an SPM contrats (obtained via the GUI on SPM)
% to a csv file
function [result_dir] = spm2csv(hReg,xSwE)
TabDat = swe_list('List',xSwE,hReg);

%% Print a csv files of the results
% if there are any results
if size(TabDat.dat,1) > 0
    myfile = 'results.csv';
    fid = fopen (myfile, 'w');
    %%
    % print header
    for j=1:size(TabDat.hdr,2)
        fprintf(fid, '%s', TabDat.hdr{1,j});
        if strcmp(TabDat.hdr{1,11},' ') fprintf(fid, '_'); end
        fprintf(fid, '%s,', TabDat.hdr{2,j});
    end
    fprintf(fid,'\n');
    %%
    % print data
    for i=1:size(TabDat.dat,1)
        for j=1:size(TabDat.dat,2)
            fprintf(fid, '%6.3f,', TabDat.dat{i,j});
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end
%%
myfile2 = 'info.csv';
fid2 = fopen (myfile2, 'w');
% print info text
currentfolder = pwd;
fprintf(fid2, '%s,', currentfolder(69:end));
fprintf(fid2, '%s,', xSwE.title);
fprintf(fid2, '%s', TabDat.tit);
fprintf(fid2, '\n');
for i = 1:size(TabDat.ftr,1)
    fprintf(fid2, TabDat.ftr{i,1},TabDat.ftr{i,2});
    fprintf(fid2, '\n\n');
end

print('...done.')
fclose(fid2);

result_dir = currentfolder;
end
