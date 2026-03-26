function plot_multi_shell_heatmap_bipolar_colored( ...
    bvec_file, bval_file, target_b_list, tol, sigma_deg, use_antipodal)
% plot_multi_shell_heatmap_bipolar_colored
%
%  複数シェルの方向をまとめて密度ヒートマップ化し、
%  点はシェルごとに色分けして表示する（密度計算は軸対称 v/-v に対応可）。
%
%  bvec_file     : 'DWI_AP.bvec' など
%  bval_file     : 'DWI_AP.bval' など
%  target_b_list : [b1 b2 ...] 例: [1500 3000], [700 2000 3000]
%  tol           : b値許容幅 (省略時 50)
%  sigma_deg     : ガウス幅 (省略時 20 度)
%  use_antipodal : true なら密度計算に -bvec も入れる (省略時 true)
%
%  例:
%   plot_multi_shell_heatmap_bipolar_colored('DWI_AP.bvec','DWI_AP.bval',[1500 3000]);

    % ---- デフォルト ----
    if nargin < 4 || isempty(tol)
        tol = 50;
    end
    if nargin < 5 || isempty(sigma_deg)
        sigma_deg = 20;
    end
    if nargin < 6 || isempty(use_antipodal)
        use_antipodal = true;
    end

    % ---- 読み込み ----
    bvec = load(bvec_file);
    bval = load(bval_file);

    % bvec -> 3×N に整形
    if size(bvec,1) ~= 3 && size(bvec,2) == 3
        bvec = bvec.';                     % N×3 -> 3×N
    elseif size(bvec,1) ~= 3
        error('bvec must be 3×N or N×3');
    end

    % bval -> 1×N に整形
    if size(bval,1) > 1 && size(bval,2) == 1
        bval = bval.';                     % 列 -> 行
    end

    % ---- 複数シェルの方向を集める ----
    N = size(bvec,2);
    mask_all = false(1,N);
    for tb = target_b_list(:).'          % 行ベクトルとしてループ
        mask_all = mask_all | (abs(bval - tb) < tol);
    end

    idx_all = find(mask_all);            % 元インデックス
    dirs_all = bvec(:, idx_all);         % 3×N_all

    if isempty(dirs_all)
        warning('No directions found for b=[%s] (tol=%.0f).', ...
                num2str(target_b_list), tol);
        return;
    end

    % ゼロ列除去
    nonzero_mask = any(dirs_all ~= 0, 1);
    dirs_all = dirs_all(:, nonzero_mask);
    idx_all  = idx_all(nonzero_mask);

    % 正規化
    dirs_all = dirs_all ./ vecnorm(dirs_all, 2, 1);
    M = size(dirs_all,2);
    fprintf('Heatmap (b=[%s], tol=%.0f): %d directions\n', ...
            num2str(target_b_list), tol, M);

    % ---- 密度計算用の方向集合 ----
    if use_antipodal
        dirs_density = [dirs_all, -dirs_all];
        tag = ' (antipodal extended)';
    else
        dirs_density = dirs_all;
        tag = '';
    end

    % ---- 球面グリッド ----
    res = 60;
    [XS,YS,ZS] = sphere(res);
    verts = [XS(:), YS(:), ZS(:)];
    K = size(verts,1);

    density = zeros(K,1);
    sigma = deg2rad(sigma_deg);

    for j = 1:K
        p  = verts(j,:).';               % 3×1
        dp = dirs_density.' * p;         % M×1
        dp = max(min(dp,1),-1);
        ang = acos(dp);
        w   = exp(-0.5*(ang/sigma).^2);
        density(j) = sum(w);
    end

    % ★ ここで平均で割って「理想=1」にする
    density_rel = density / mean(density);

    density_grid = reshape(density_rel, size(XS));

    % ---- 描画開始 ----
    figure('Position',[100 100 650 650]);

    % ヒートマップ
    surf(XS,YS,ZS, density_grid, ...
         'FaceAlpha',0.95, ...
         'EdgeColor','none');
    hold on;

    % ---- シェルごとの色分け（点プロット＋凡例用ダミー） ----
    nShell      = numel(target_b_list);
    cmap        = lines(max(nShell,3));
    leg_handles = gobjects(0);
    leg_labels  = {};

    for i = 1:nShell
        tb = target_b_list(i);
        mask_b_local = abs(bval(idx_all) - tb) < tol;
        dirs_b = dirs_all(:, mask_b_local);

        if ~isempty(dirs_b)
            % ① 実データ点
            plot3(dirs_b(1,:), dirs_b(2,:), dirs_b(3,:), ...
                  '.', 'MarkerSize',18, 'Color',cmap(i,:));

            % ② 反対方向（bipolar）の点
            plot3(-dirs_b(1,:), -dirs_b(2,:), -dirs_b(3,:), ...
                  '.', 'MarkerSize',18, 'Color',cmap(i,:));

            % 凡例用ダミー（1 回だけ表示）
            h_leg = plot3(nan, nan, nan, '.', ...
                          'MarkerSize',18, 'Color',cmap(i,:));
            leg_handles(end+1) = h_leg;
            leg_labels{end+1}  = sprintf('b = %.0f (bipolar)', tb);

        end
    end

    colormap(turbo);
    cb = colorbar;

    % 例：±10% のズレを表現するスケールに固定
    caxis([0.9 1.1]);
    cb.Label.String = 'relative density ( / mean )';

    % カラーバーの位置を右端に移動（タイトルと重ならない）
    cb.Position = [0.88 0.20 0.03 0.60];

    axis([-1 1 -1 1 -1 1]);
    axis vis3d;
    daspect([1 1 1]);
    pbaspect([1 1 1]);
    camproj('orthographic');
    grid on;

    xlabel x; ylabel y; zlabel z;
    title(sprintf('Directional density b=[%s]%s (\\sigma=%d^\\circ)', ...
           num2str(target_b_list), tag, sigma_deg));

    if ~isempty(leg_handles)
        legend(leg_handles, leg_labels, 'Location','bestoutside');
    end

    view(30,20);
    rotate3d on;

    % ---- 軽く1回転（元のアニメーション）----
    for az = 0:10:360
        if ~ishandle(gcf), break; end
        view(az,20);
        drawnow;
        pause(0.05);
    end

    % ==== ここから「最後に付け加える」スクリーンショット部分 ====

    % 見栄えのする斜め上アングルに固定（例：45°/30°）
    view(45,30);

    % b値リストからファイル名を自動生成
    out_png = sprintf('heatmap_b[%s].png', ...
        strrep(num2str(target_b_list),' ','_'));

    % 高解像度で自動保存（exportgraphics があれば使用）
    try
        exportgraphics(gcf, out_png, 'Resolution',300);
    catch
        % 古い MATLAB 用のフォールバック
        print(gcf, out_png, '-dpng', '-r300');
    end

    fprintf('Saved screenshot: %s\n', out_png);

end
