function [dll,matH,coord_x,ir,irr,ird] = member(dummy)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coord_x (nk, 2)
% --->
coord_x = [
    1, 0;
    3, 0;
    5, 1;
    0, 1;
    2, 2;
    0, 3;
    0, 5;
    2, 5];
coord_x =  1.0 * coord_x;
nk = size(coord_x,1);

%%%% irr(nm, 2)
irr = [
    1, 2;
    1, 4;
    1, 5;
    2, 3;
    2, 5;
    3, 5;
    4, 5;
    4, 6;
    5, 6;
    5, 8;
    6, 7;
    6, 8];
nm = size(irr,1);

%%%% ird(nk, 3)
ird = ones(nk,3);
ird(end-1,:) = [0, 0, 0];
ird(end,:)   = [0, 0, 0];
nd = sum(sum(ird));
%
ii = 0;
for j=1:nk
    for k=1:3
        if ird(j,k) == 1
            ii = ii + 1;
            ird(j,k) = ii;
        else
            ird(j,k) = nd +1;
        end
    end
end

%%%% ir(nm, 6)
ir = zeros(nm,6);
for i=1:nm
    for j=1:3
        ir(i,j)   = ird(irr(i,1), j);
        ir(i,j+3) = ird(irr(i,2), j);
    end
end

%%%% matH(nd, nm)
dll = zeros(nm,1);
matH{1} = zeros(nd+1,nm);
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    dx = coord_x(j2,1) - coord_x(j1,1);
    dy = coord_x(j2,2) - coord_x(j1,2);
    dll(i) = norm([dx; dy], 2);
    dir_cos(1) = dx/dll(i);
    dir_cos(2) = dy/dll(i);
    for j=1:2
        if abs(dir_cos(j)) < 10^(-16)
            dir_cos(j) = 0;
        end
        matH{1}(ir(i,j),i)    = -dir_cos(j);
        matH{1}(ir(i,j+3), i) =  dir_cos(j);
    end
    clear dir_cos
end
matH{1} = sparse(matH{1}(1:nd,:));

matH{2} = zeros(nd+1,nm);
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    dx = coord_x(j2,1) - coord_x(j1,1);
    dy = coord_x(j2,2) - coord_x(j1,2);
    dir_cos2(1) =  2*dy/(dll(i)^2);
    dir_cos2(2) = -2*dx/(dll(i)^2);
    for j=1:2
        if abs(dir_cos2(j)) < 10^(-16)
            dir_cos2(j) = 0;
        end
        matH{2}(ir(i,j),i)    = -dir_cos2(j);
        matH{2}(ir(i,j+3), i) =  dir_cos2(j);
    end
    matH{2}(ir(i,3), i) = 1;
    matH{2}(ir(i,6), i) = 1;
end
matH{2} = sparse(matH{2}(1:nd,:));

matH{3} = zeros(nd+1,nm);
for i=1:nm
    matH{3}(ir(i,3), i) = -1;
    matH{3}(ir(i,6), i) =  1;
end
matH{3} = sparse(matH{3}(1:nd,:));

