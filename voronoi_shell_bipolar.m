function stats = voronoi_shell_bipolar(bvec_file, bval_file, target_b_list, tol)
% voronoi_shell_bipolar  軸対称 (v/-v) を考慮した Voronoi 面積評価
%
%  stats = voronoi_shell_bipolar(bvec_file, bval_file, target_b_list, tol)
%
%  bvec_file     : bvec ファイル名 (例 'DWI_AP.bvec')
%  bval_file     : bval ファイル名 (例 'DWI_AP.bval')
%  target_b_list : 評価したい b 値 (スカラー or ベクトル) 例: 2000, [1500 3000]
%  tol           : b 値許容幅 (省略時 50)
%
%  戻り値 stats:
%    .N_axis      : 軸の本数 (M)
%    .vertex_area : 各軸方向の Voronoi 面積 (M×1)
%    .total_area  : sum(vertex_area) (≈ 4π)
%    .mean, .median, .std, .CV, .min, .max, .maxmin_ratio
%

    if nargin < 4 || isempty(tol)
        tol = 50;
    end

    % ----- load -----
    bvec = load(bvec_file);   % 3×N or N×3
    bval = load(bval_file);   % 1×N or N×1

    % bvec -> 3×N
    if size(bvec,1) ~= 3 && size(bvec,2) == 3
        bvec = bvec.';        % N×3 -> 3×N
    elseif size(bvec,1) ~= 3
        error('bvec must be 3×N or N×3');
    end

    % bval -> 1×N
    if size(bval,1) > 1 && size(bval,2) == 1
        bval = bval.';        % 列 -> 行
    end

    % ----- b シェルのマスク (複数対応) -----
    mask_shell = false(size(bval));
    for tb = target_b_list(:).'
        mask_shell = mask_shell | (abs(bval - tb) < tol);
    end

    dirs = bvec(:, mask_shell);   % 3×N_shell

    if isempty(dirs)
        error('No directions found for b = [%s] with tol=%.0f.', ...
              num2str(target_b_list), tol);
    end

    % ゼロ列除外 & 正規化
    dirs = dirs(:, any(dirs ~= 0,1));
    dirs = dirs ./ vecnorm(dirs, 2, 1);

    M = size(dirs,2);  % 軸本数
    fprintf('Voronoi (b=[%s], tol=%.0f): M = %d axis directions\n', ...
            num2str(target_b_list), tol, M);

    % ----- Voronoi 用に ±方向を用いた点集合 -----
    dirs_full = [dirs, -dirs];   % 3×(2M)
    nV = size(dirs_full,2);

    % 球面上の凸包
    faces = convhull(dirs_full(1,:).', dirs_full(2,:).', dirs_full(3,:).');  % F×3
    F = size(faces,1);

    vertex_area_full = zeros(nV,1);

    for f = 1:F
        idx = faces(f,:);
        v1  = dirs_full(:, idx(1));
        v2  = dirs_full(:, idx(2));
        v3  = dirs_full(:, idx(3));

        % 大円距離
        a = acos(max(min(dot(v2,v3),1),-1));
        b = acos(max(min(dot(v1,v3),1),-1));
        c = acos(max(min(dot(v1,v2),1),-1));

        % L'Huilier の公式で球面三角形面積
        s = (a + b + c) / 2;
        t = tan(s/2) .* tan((s-a)/2) .* tan((s-b)/2) .* tan((s-c)/2);
        t = max(t, 0);  % 数値誤差対策
        E = 4 * atan(sqrt(t));   % 面積 (単位球)

        vertex_area_full(idx) = vertex_area_full(idx) + E/3;
    end

    % ----- 軸として +v と -v をペアにして統合 -----
    % dirs_full = [dirs(1..M), -dirs(1..M)] という並びなので
    vertex_area_axis = vertex_area_full(1:M) + vertex_area_full(M+1:end);

    total_area = sum(vertex_area_axis);
    mean_area  = mean(vertex_area_axis);
    std_area   = std(vertex_area_axis);
    CV         = std_area / mean_area;
    min_area   = min(vertex_area_axis);
    max_area   = max(vertex_area_axis);
    maxmin_ratio = max_area / min_area;

    fprintf('sum(area)     = %.4f (4*pi = %.4f)\n', total_area, 4*pi);
    fprintf('area min/max  = %.4e / %.4e\n', min_area, max_area);
    fprintf('area mean     = %.4e (ideal=%.4e)\n', mean_area, 4*pi/M);
    fprintf('area median   = %.4e\n', median(vertex_area_axis));
    fprintf('area std      = %.4e\n', std_area);
    fprintf('CV (std/mean) = %.4f\n', CV);
    fprintf('max/min ratio = %.2f\n\n', maxmin_ratio);

    % 出力 struct
    stats.N_axis       = M;
    stats.vertex_area  = vertex_area_axis;
    stats.total_area   = total_area;
    stats.mean         = mean_area;
    stats.median       = median(vertex_area_axis);
    stats.std          = std_area;
    stats.CV           = CV;
    stats.min          = min_area;
    stats.max          = max_area;
    stats.maxmin_ratio = maxmin_ratio;
    % ★ 追加：軸方向（3×M）も返す
    stats.axis_dirs    = dirs;
end

