function r = sample_POS_ONLY_update( t_elapsed,D, n_sigma )
% t_elapsed and D may be vectors.
% Assumes system has no boundaries, r ~ [0,inf] available.

% NO REJECTION SAMPLING NEEDED.
% SAMPLE r ~ [original probabiilty density]
r = abs( normrnd(0, (6.*D.*t_elapsed).^(1/2)) );  

% n_sigma = 3;
R_diff_sphere = n_sigma.*sqrt(6.*D.*t_elapsed);

% Testing: throw error if needed. After testing, ccomment this out
% if r > R_diff_sphere
%    display(r)
%    display(R_diff_sphere)
%    assert(0==1)
%     
% end


% Restrict r to be less than R_diff_sphere boundary
list = r > R_diff_sphere;
if sum(list)>0
    r(list) = 0.99 * R_diff_sphere(list);
end
% if r > R_diff_sphere
%     r = 0.97 * R_diff_sphere;
% end


end

