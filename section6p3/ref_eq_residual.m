function F = ref_eq_residual(vec_u)

global vec_p matL matN ref_Young ref_back intersec_pt
global ref_Young_no ref_back_no vec_init_strain intersec_pt_no

nm = size(matL, 1);

vec_eps = (matL * vec_u) + vec_init_strain;

vec_sig = zeros(nm,1);
for i=1:(nm-3)
    strain = vec_eps(i);
    if strain <= intersec_pt_no{1}(1)
        vec_sig(i) = (ref_Young_no(1) * strain) + ref_back_no(1);
    elseif strain >= intersec_pt_no{2}(1)
        vec_sig(i) = (ref_Young_no(3) * strain) + ref_back_no(3);
    else
        vec_sig(i) = (ref_Young_no(2) * strain) + ref_back_no(2);
    end
end
for i=(nm-2):nm
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

