function obj_func = GaussMix2Distr_regul(params,X,Y,sigma_x, sigma_y, fx, fy, reg)


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




n1 = size(X, 1);
n2 = size(Y, 1);

Z = convert_optimal(Y, params);

d_xz = [];
for i = 1:size(X, 1)
    d_xz(i,:) = pdist_eucledian(Z, X(i,:)) .* fy * fx(i);
end

d_xz = sum(min(d_xz) .* fy') + sum(min(d_xz') .* fx');



obj_func = 0;
for i = 1:n1
for j = 1:n2
    point = (X(i,:) - Z(j,:))';
    sigma_joint = sigma_x(:,:,i) + R * sigma_y(:,:,j) * R';
%     sigma_joint = sigma_x(:,:,i) + R * sigma_y(:,:,j) * R';
    obj_func = obj_func - 2 * fx(i)*fy(j) *...
        mvnpdf(point,zeros(n, 1),sigma_joint);
end
end

obj_func_x = 0;

for i = 1:n1
for j = 1:n1
    point = (X(i,:) - X(j,:))';
    sigma_joint = sigma_x(:,:,i) + sigma_x(:,:,j);
%     sigma_joint = sigma_x(:,:,i) + R * sigma_y(:,:,j) * R';
    obj_func_x = obj_func_x + fx(i)*fx(j) *...
        mvnpdf(point,zeros(n, 1),sigma_joint);
end
end
obj_func_x;

obj_func_y = 0;
for i = 1:n2
for j = 1:n2
    point = (Z(i,:) - Z(j,:))';
    sigma_joint = sigma_y(:,:,i) + sigma_y(:,:,j);
%     sigma_joint = sigma_x(:,:,i) + R * sigma_y(:,:,j) * R';
    obj_func_y = obj_func_y + fy(i)*fy(j) *...
        mvnpdf(point,zeros(n, 1),sigma_joint);
end
end
obj_func_y;

% obj_func =  obj_func + obj_func_x + obj_func_y; 


