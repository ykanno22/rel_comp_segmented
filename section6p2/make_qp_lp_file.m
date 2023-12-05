function [dummy] = make_qp_lp_file(lpfile,vec_obj,con_ax,con_bx,lb,ub,...
    con_ax_s,con_bx_s,vartype,ctype,con_Qx_q,con_rx_q,con_bx_q)
%
dummy = 1;
%
fid = fopen(lpfile, 'w');
fprintf(fid, 'Minimize\n');
fprintf(fid, ' obj: ');
for j=1:size(con_ax,2)
    if abs(vec_obj(j)) > 10^(-10)
        if vec_obj(j) > 0
            fprintf(fid, '  + %.14f x_%g \n', full(vec_obj(j)), j);
        else
            fprintf(fid, '  - %.14f x_%g \n', full(abs(vec_obj(j))), j);
        end
    end
end

fprintf(fid, '\n');
fprintf(fid, 'Subject to\n');
for i=1:size(con_ax,1)
    fprintf(fid, 'r_%g: ', i);
    for j=1:size(con_ax,2)
        if abs(con_ax(i,j)) > 10^(-10)
            if con_ax(i,j) > 0
                fprintf(fid, '  + %.14f x_%g \n', full(con_ax(i,j)), j);
            else
                fprintf(fid, '  - %.14f x_%g \n', full(abs(con_ax(i,j))), j);
            end
        end
    end
    if ctype(i) == 'U'
        fprintf(fid, '    <= %.14f\n', full(con_bx(i)));
    elseif ctype(i) == 'S'
        fprintf(fid, '    = %.14f\n', full(con_bx(i)));
    end
end

for i=1:length(con_ax_s)
    fprintf(fid, 'r_%g: ', length(con_bx)+i);
    fprintf(fid, '  [');
    [jjj,kkk,sss] = find(con_ax_s{i});
    for l=1:length(jjj)
        if sss(l) > 0
            fprintf(fid, '\n  + %.14f x_%g * x_%g', sss(l), jjj(l), kkk(l));
        else
            fprintf(fid, '\n  - %.14f x_%g * x_%g', abs(sss(l)), jjj(l), kkk(l));
        end
    end
    fprintf(fid, '  ]');
    fprintf(fid, '    <= %.14f\n', full(con_bx_s(i)));
end

for i=1:length(con_Qx_q)
    fprintf(fid, 'r_%g: ', length(con_bx)+length(con_ax_s)+i);
    fprintf(fid, '  [');
    [jjj,kkk,sss] = find(con_Qx_q{i});
    for l=1:length(jjj)
        if sss(l) > 0
            fprintf(fid, '\n  + %.14f x_%g * x_%g', sss(l), jjj(l), kkk(l));
        else
            fprintf(fid, '\n  - %.14f x_%g * x_%g', abs(sss(l)), jjj(l), kkk(l));
        end
    end
    fprintf(fid, '  ]');
    [jjj,~,sss] = find(con_rx_q{i});
    if ~isempty(jjj)
        for l=1:length(jjj)
            if sss(l) > 0
                fprintf(fid, '\n  + %.14f x_%g', sss(l), jjj(l));
            else
                fprintf(fid, '\n  - %.14f x_%g', abs(sss(l)), jjj(l));
            end
        end
    end
    fprintf(fid, '    <= %.14f\n', full(con_bx_q(i)));
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
