function [R, ctx] = import_ctx_data(id)

    location = id(1:3);
    ptnum = str2num(id(4:7));
    
    
    %BWH M100 and SHS
    if (strcmp(location,'BWH') && ptnum <= 2266) || strcmp(location,'SHS')
        prefix1 = '049';
        prefix2 = '052';
    %BWH M1010
    elseif strcmp(location,'BWH') || strcmp(location,'SSH')
        prefix1 = '';
        prefix2 = '';

%        disp('cannot load data from BWH M1010 yet');
%        return;

    elseif strcmp(location,'CCF')
        prefix1 = '001';
        prefix2 = '051';
    end

    % load the raw data
    [A,B] = read_ctx_data(id, location, prefix1, prefix2);
    
    % correct the toco/IUPC values
    if (strcmp(location,'BWH') && ptnum <= 2266) || strcmp(location,'CCF')
        %ctx = ((2.048*(B(29,:) - 32767)/32767)+0.862)*(100/1.2);
        ctx = 0.0052*B(29,:) - 98.551;
    elseif strcmp(location,'BWH')
        disp('do not know how to scale these ctx values');
        ctx = (B(29,:)-min(B(29,:)))/(max(B(29,:))-min(B(29,:)))*120+10;
    else
        disp('do not know how to scale these ctx values');
        ctx = 0.0052*B(29,:) - 98.551;
    end
    
    % correct the abdominal values
    R = 2.048*(A-32767)/32767;

end

function [A,B] = read_ctx_data(id, location, prefix1, prefix2)

%    local_name1 = ['CTXDevelData\' prefix1 id '000.dat'];
    local_name1 = [id '/' prefix1 id '000.dat'];

    if exist(local_name1,'file')
        A = load_dat(local_name1);
    else
        disp('No local copy (049)');
        A = load_dat(['Y:\' location '\' id '\' prefix1 id '000.dat']);
    end
    
%    local_name2 = ['CTXDevelData\' prefix2 id '000.dat'];
    local_name2 = [id '/' prefix2 id '000.dat'];

    if exist(local_name2,'file')
        
        B = load_dat(local_name2);
    else
        disp('No local copy (052)');
        B = load_dat(['Y:\' location '\' id '\' prefix2 id '000.dat']);
    end
    
end