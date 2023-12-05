function [dll,matH,coord_x,ir,irr,ird] = member(dummy)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coord_x (nk, 2)
% --->
x_num{1} = 3;
x_num{2} = 3;

for jj=1:size(x_num,2)
    xx{jj} = [];
end

for j=1:x_num{1}
    xx{1} = [xx{1};...
        1.0 * cos(pi*(((2*(j-1))/x_num{1}) + (0/12))),...
        1.0 * sin(pi*(((2*(j-1))/x_num{1}) + (0/12)))];
end
for j=1:x_num{2}
    xx{2} = [xx{2};...
        1.0 * cos(pi*(((2*(j-1))/x_num{2}) + (-1/4))),...
        1.0 * sin(pi*(((2*(j-1))/x_num{2}) + (-1/4)))];
end

coord_x = [];
zz = 0;
for jj=1:size(x_num,2)
    coord_x = [coord_x;...
        xx{jj}, zz*ones(x_num{jj},1)];
    zz = zz + 1.5;
end

nk = size(coord_x,1);

%%%% irr(nm, 2)
irr = [
    1, 2;
    2, 3;
    3, 1;
    4, 5;
    5, 6;
    6, 4;
    1, 4;
    2, 5;
    3, 6;
    1, 5;
    2, 6;
    3, 4;
    1, 6;
    2, 4;
    3, 5];

nm = size(irr,1);

%%%% ird(nk, 3)
ird = ones(nk,3);
ird(1,:) = [0, 0, 0];
ird(2,:) = [1, 0, 0];
ird(3,:) = [1, 1, 0];
%
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
%
%%%% ir(nm, 6)
ir = zeros(nm,6);
for i=1:nm
    for j=1:3
        ir(i,j)   = ird(irr(i,1), j);
        ir(i,j+3) = ird(irr(i,2), j);
    end
end
%
%%%% matH(nd, nm)
dll = zeros(nm,1);
matH = sparse(zeros(nd+1,nm));
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    dx = coord_x(j2,1) - coord_x(j1,1);
    dy = coord_x(j2,2) - coord_x(j1,2);
    dz = coord_x(j2,3) - coord_x(j1,3);
    dll(i) = norm([dx; dy; dz], 2);
    dir_cos(1) =-dx/dll(i);
    dir_cos(2) =-dy/dll(i);
    dir_cos(3) =-dz/dll(i);
    dir_cos(4) = dx/dll(i);
    dir_cos(5) = dy/dll(i);
    dir_cos(6) = dz/dll(i);
    for j=1:6
        if abs(dir_cos(j)) < 10^(-15)
            dir_cos(j) = 0;
        end
        matH(ir(i,j),i) = dir_cos(j);
    end
end
%
matH = matH(1:nd,:);

