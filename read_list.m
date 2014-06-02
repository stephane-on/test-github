function list_sacfiles=read_list(fich)

fid=fopen(fich,'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~exist('list_sacfiles')
        list_sacfiles=tline;
    else
        list_sacfiles=char(list_sacfiles,tline);
    end
end
fclose(fid);