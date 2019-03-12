function ncdf2mat(ncname)
x=ncinfo(ncname);
x_cell=squeeze(struct2cell(x.Variables));
var_name=x_cell(1,:);

for i=1:length(var_name);
    var_here=var_name{i};
    data_here=ncread(ncname,var_here);
    save(var_here,'data_here');
end