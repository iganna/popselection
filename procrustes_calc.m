function Z = procrustes_calc(Y, params)

n = size(Y, 2);

n_angles = n*(n-1)/2;
angles = params(1:n_angles);
rotation_zero = diag(ones(n, 1));
R = rotation_zero;
id = 1;
for i = 1:n
    for j = (i+1):n
        rot_tmp = rotation_zero;
        rot_tmp([i j], [i, j]) = [cos(angles(id)) -sin(angles(id));...
            sin(angles(id)), cos(angles(id))];
        id = id + 1;
        R = rot_tmp * R;
    end
end

% scale_params = diag(params(7:9));

t = params(n_angles + (1:n));

% uni_scale = params(n_angles + n + 1)';
uni_scale = 1;

for i = 1:size(Y,1)
    Z(i,:) = uni_scale * (R*Y(i,:)')' + t;
end