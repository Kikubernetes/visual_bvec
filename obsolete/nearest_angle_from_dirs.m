function nn = nearest_angle_from_dirs(axis_dirs)
% axis_dirs : 3×M の単位ベクトル（voronoi_shell_bipolar の stats.axis_dirs）
%
% nn : struct
%   .theta      : 最近傍角度 [deg] (M×1)
%   .min, .max
%   .mean, .median
%
    % 念のため正規化
    axis_dirs = axis_dirs ./ vecnorm(axis_dirs, 2, 1);

    % 内積行列 G(i,j) = v_i · v_j
    G = axis_dirs.' * axis_dirs;  % M×M
    M = size(G,1);

    % 自身 (対角) を除外
    G(1:M+1:end) = -Inf;

    % 最近傍の cosθ
    [cos_nearest, ~] = max(G, [], 2);   % M×1
    theta = acosd(cos_nearest);        % deg

    nn.theta  = theta;
    nn.min    = min(theta);
    nn.max    = max(theta);
    nn.mean   = mean(theta);
    nn.median = median(theta);
end

