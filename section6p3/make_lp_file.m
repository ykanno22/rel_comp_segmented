function [dummy] = make_lp_file(lpfile,vec_object,con_ax,con_bx,lb,ub,vartype,ctype)
%
dummy = 1;
%
    fid = fopen(lpfile, 'w');
    fprintf(fid, 'Minimize\n');
    fprintf(fid, ' obj: ');
    for j=1:size(con_ax,2)
        if abs(vec_object(j)) > 10^(-10)
            if vec_object(j) > 0
                fprintf(fid, '  + %.14f x_%g \n', full(vec_object(j)), j);
            else
                fprintf(fid, '  - %.14f x_%g \n', full(abs(vec_object(j))), j);
            end
        end
    end
    fprintf(fid, '\nSubject to\n');
    for i=1:size(con_ax,1)
        fprintf(fid, 'r_%g: ', i);
        for j=1:size(con_ax,2)
            if abs(con_ax(i,j)) > 10^(-10)
                if con_ax(i,j) > 0
                    fprintf(fid, '  + %.14f x_%g \n', full(con_ax(i,j)), j);
                else
                    fprintf(fid, '  - %.14f x_%g \n', abs(full(con_ax(i,j))), j);
                end
            end
        end
        if ctype(i) == 'U'
            fprintf(fid, '    <= %.14f\n', full(con_bx(i)));
        elseif ctype(i) == 'S'
            fprintf(fid, '    = %.14f\n', full(con_bx(i)));
        end
    end
    fprintf(fid, '\nBounds\n');
    for j=1:size(con_ax,2)
        fprintf(fid, ' %.14f <= x_%g <= %.14f\n', lb(j), j, ub(j));
    end
    fprintf(fid, '\nBinaries\n');
    for j=1:size(con_ax,2)
        if vartype(j) == 'B'
            fprintf(fid, ' x_%g\n', j);
        end
    end
    fprintf(fid, '\nEnd\n');
    fclose(fid);
