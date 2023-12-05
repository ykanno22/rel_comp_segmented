function F = ref_eq_residual(vec_u)

global vec_p matL matN ref_Young ref_back intersec_pt

nm = size(matL, 1);

vec_eps = matL * vec_u;

vec_sig = zeros(nm,1);
for i=1:nm
    strain = vec_eps(i);
    if strain <= intersec_pt{1}(1)
        vec_sig(i) = (ref_Young(1) * strain) + ref_back(1);
    elseif strain >= intersec_pt{2}(1)
        vec_sig(i) = (ref_Young(3) * strain) + ref_back(3);
    else
        vec_sig(i) = (ref_Young(2) * strain) + ref_back(2);
    end
end

F = (matN * vec_sig) - vec_p;

