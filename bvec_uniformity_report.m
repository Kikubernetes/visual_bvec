function report = bvec_uniformity_report( ...
    dataset_name, bvec_file, bval_file, target_b_list, tol, Lmax, html_name)
% UNIFORMITY_REPORT_ONE_DATASET (HTML version)
%  球面ボロノイ・SH・最近傍角度・electrostatic energy の指標を
%  table + HTML レポートで出力する。
%
% 入力:
%   dataset_name  : データセット名
%   bvec_file     : bvec
%   bval_file     : bval
%   target_b_list : 評価したい b 値
%   tol           : b 値許容幅 (default 50)
%   Lmax          : SH 最大次数 (default 8)
%   html_name     : 出力ベース名
%                   例) 'dir98_uniformity' →
%                       ./dir98_uniformity/dir98_uniformity.html
%
% 出力:
%   report : base.html
%            graphs.png
%            table.xlsx
%  

    if nargin < 5 || isempty(tol)
        tol = 50;
    end
    if nargin < 6 || isempty(Lmax)
        Lmax = 8;
    end
    if nargin < 7
        html_name = '';
    end

    target_b_list = target_b_list(:).';

    rows          = [];
    shell_results = [];
    all_dirs_axis = [];   % 合算シェル用に全方向を貯める

    % ==========================
    % 出力フォルダの決定
    % ==========================
    html_path = '';
    out_dir   = '';
    out_base  = '';

    if ~isempty(html_name)
        [parent_dir, base, ~] = fileparts(html_name);
        if isempty(parent_dir)
            parent_dir = '.';
        end
        if isempty(base)
            error('html_name が不正です（ベース名が空）。');
        end

        out_dir = fullfile(parent_dir, base);
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end

        out_base  = base;
        html_path = fullfile(out_dir, [base '.html']);
    end

    % ==========================
    % 各 shell の計算
    % ==========================
    for ib = 1:numel(target_b_list)
        tb = target_b_list(ib);

        vs = voronoi_shell_bipolar(bvec_file, bval_file, tb, tol);
        dirs_axis = vs.axis_dirs;          % 3×M
        all_dirs_axis = [all_dirs_axis, dirs_axis]; %#ok<AGROW>

        [Lvals, Pl, UI] = sh_uniformity_from_dirs(dirs_axis, Lmax);
        nn = nearest_angle_from_dirs(dirs_axis);
        es = electrostatic_energy_from_dirs(dirs_axis, 1);

        % --- table 行 ---
        row = struct();
        row.Dataset        = string(dataset_name);
        row.Shell_b        = tb;
        row.N_axis         = vs.N_axis;
        row.Vor_CV         = vs.CV;
        row.Vor_MaxMin     = vs.maxmin_ratio;
        row.SH_UI          = UI;
        row.NN_min_deg     = nn.min;
        row.NN_mean_deg    = nn.mean;
        row.NN_median_deg  = nn.median;
        row.NN_max_deg     = nn.max;
        row.ES_Emean       = es.E_meanpair;

        rows = [rows; row]; %#ok<AGROW>

        % --- HTML / PNG 用情報 ---
        info = struct();
        info.Shell_b        = tb;
        info.N_axis         = vs.N_axis;
        info.Vor_CV         = vs.CV;
        info.Vor_MaxMin     = vs.maxmin_ratio;
        info.SH_UI          = UI;
        info.NN_min_deg     = nn.min;
        info.NN_mean_deg    = nn.mean;
        info.NN_median_deg  = nn.median;
        info.NN_max_deg     = nn.max;
        info.ES_Emean       = es.E_meanpair;
        info.voronoi_png    = '';
        info.sh_png         = '';
        info.nn_png         = '';
        info.is_combined    = false;

        if ~isempty(html_path)
            prefix = fullfile(out_dir, sprintf('%s_shell%02d_b%d', ...
                out_base, ib, round(tb)));

            % Voronoi
            fig = figure('Visible','off');
            histogram(vs.vertex_area, 20)
            xlabel('Voronoi cell area')
            ylabel('Count')
            title(sprintf('%s  b≈%g (CV=%.3f)', dataset_name, tb, vs.CV))
            grid on
            vor_png = [prefix '_voronoi.png'];
            exportgraphics(fig, vor_png, 'Resolution', 200);
            close(fig);

            % SH
            fig = figure('Visible','off');
            stem(Lvals, Pl, 'filled')
            xlabel('Degree l')
            ylabel('Power P_l')
            title(sprintf('%s  b≈%g  SH (UI=%.3f)', dataset_name, tb, UI))
            grid on
            sh_png = [prefix '_sh.png'];
            exportgraphics(fig, sh_png, 'Resolution', 200);
            close(fig);

            % NN
            fig = figure('Visible','off');
            histogram(nn.theta, 0:2:60)
            xlabel('Nearest-neighbor angle (deg)')
            ylabel('Count')
            title(sprintf('%s b≈%g  NN(mean=%.2f°)', dataset_name, tb, nn.mean))
            grid on
            nn_png = [prefix '_nn.png'];
            exportgraphics(fig, nn_png, 'Resolution', 200);
            close(fig);

            info.voronoi_png = vor_png;
            info.sh_png      = sh_png;
            info.nn_png      = nn_png;
        end

        shell_results = [shell_results; info]; %#ok<AGROW>
    end

    % ==========================
    % 合算シェル（全方向）を追加
    % ==========================
    if ~isempty(all_dirs_axis) && numel(target_b_list) > 1
        % Voronoi は target_b_list をまとめて指定して再計算
        vs_all = voronoi_shell_bipolar(bvec_file, bval_file, target_b_list, tol);
        dirs_all = vs_all.axis_dirs;

        [Lvals_all, Pl_all, UI_all] = sh_uniformity_from_dirs(dirs_all, Lmax);
        nn_all = nearest_angle_from_dirs(dirs_all);
        es_all = electrostatic_energy_from_dirs(dirs_all, 1);

        row = struct();
        row.Dataset        = string(dataset_name) + " (combined)";
        row.Shell_b        = "combined";
        row.N_axis         = vs_all.N_axis;
        row.Vor_CV         = vs_all.CV;
        row.Vor_MaxMin     = vs_all.maxmin_ratio;
        row.SH_UI          = UI_all;
        row.NN_min_deg     = nn_all.min;
        row.NN_mean_deg    = nn_all.mean;
        row.NN_median_deg  = nn_all.median;
        row.NN_max_deg     = nn_all.max;
        row.ES_Emean       = es_all.E_meanpair;

        rows = [rows; row]; %#ok<AGROW>

        if ~isempty(html_path)
            prefix = fullfile(out_dir, sprintf('%s_combined', out_base));

            % Voronoi（合算）
            fig = figure('Visible','off');
            histogram(vs_all.vertex_area, 20)
            xlabel('Voronoi cell area')
            ylabel('Count')
            title(sprintf('%s combined Voronoi (CV=%.3f)', dataset_name, vs_all.CV))
            grid on
            vor_png_all = [prefix '_voronoi.png'];
            exportgraphics(fig, vor_png_all, 'Resolution', 200);
            close(fig);

            % SH（合算）
            fig = figure('Visible','off');
            stem(Lvals_all, Pl_all, 'filled')
            xlabel('Degree l')
            ylabel('Power P_l')
            title(sprintf('%s  combined SH (UI=%.3f)', dataset_name, UI_all))
            grid on
            sh_png_all = [prefix '_sh.png'];
            exportgraphics(fig, sh_png_all, 'Resolution', 200);
            close(fig);

            % NN（合算）
            fig = figure('Visible','off');
            histogram(nn_all.theta, 0:2:60)
            xlabel('Nearest-neighbor angle (deg)')
            ylabel('Count')
            title(sprintf('%s combined NN(mean=%.2f°)', dataset_name, nn_all.mean))
            grid on
            nn_png_all = [prefix '_nn.png'];
            exportgraphics(fig, nn_png_all, 'Resolution', 200);
            close(fig);

            info = struct();
            info.Shell_b        = -1;
            info.N_axis         = vs_all.N_axis;
            info.Vor_CV         = vs_all.CV;
            info.Vor_MaxMin     = vs_all.maxmin_ratio;
            info.SH_UI          = UI_all;
            info.NN_min_deg     = nn_all.min;
            info.NN_mean_deg    = nn_all.mean;
            info.NN_median_deg  = nn_all.median;
            info.NN_max_deg     = nn_all.max;
            info.ES_Emean       = es_all.E_meanpair;
            info.voronoi_png    = vor_png_all;
            info.sh_png         = sh_png_all;
            info.nn_png         = nn_png_all;
            info.is_combined    = true;

            shell_results = [shell_results; info]; %#ok<AGROW>
        end
    end

    % ==========================
    % table 化
    % ==========================
    report = struct2table(rows);

    % Excel 出力
    if ~isempty(html_name)
        excel_path = fullfile(out_dir, [out_base '.xlsx']);
    else
        excel_path = [char(string(dataset_name)) '_uniformity_report.xlsx'];
    end
    writetable(report, excel_path);

    % ==========================
    % HTML レポート作成
    % ==========================
    if ~isempty(html_path)
        write_uniformity_html_report(html_path, dataset_name, report, shell_results);
    end

    disp(report);
end

% ============================================================
% サブ関数: HTML レポート作成
% ============================================================
function write_uniformity_html_report(html_name, dataset_name, report, shell_results)
    if isempty(html_name)
        return;
    end

    fid = fopen(html_name, 'w');
    if fid == -1
        warning('HTML ファイル %s を開けませんでした。', html_name);
        return;
    end

    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '<!DOCTYPE html>\n<html lang="ja">\n<head>\n');
    fprintf(fid, '<meta charset="UTF-8">\n');
    fprintf(fid, '<title>%s – Directional uniformity report</title>\n', dataset_name);

    fprintf(fid, '<style>\n');
    fprintf(fid, 'body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 20px; line-height: 1.5; }\n');
    fprintf(fid, 'h1 { font-size: 1.8em; margin-bottom: 0.2em; }\n');
    fprintf(fid, 'h2 { margin-top: 1.4em; font-size: 1.3em; }\n');
    fprintf(fid, 'table { border-collapse: collapse; margin-top: 0.5em; }\n');
    fprintf(fid, 'th, td { border: 1px solid #ccc; padding: 4px 8px; font-size: 0.85em; text-align: right; }\n');
    fprintf(fid, 'th { background-color: #f5f5f5; text-align: center; }\n');
    fprintf(fid, '.note { font-size: 0.85em; color: #555; max-width: 900px; }\n');
    fprintf(fid, '.panel { display: flex; flex-wrap: wrap; gap: 16px; margin-top: 0.5em; }\n');
    fprintf(fid, '.panel figure { flex: 1 1 260px; margin: 0; }\n');
    fprintf(fid, '.panel img { max-width: 100%%; height: auto; border: 1px solid #ddd; }\n');
    fprintf(fid, 'figcaption { font-size: 0.8em; color: #444; margin-top: 4px; }\n');
    fprintf(fid, '</style>\n</head>\n<body>\n');

    fprintf(fid, '<h1>%s – Directional uniformity report</h1>\n', dataset_name);

    % Summary table
    fprintf(fid, '<h2>Summary table</h2>\n');
    fprintf(fid, '<table>\n');

    varNames = report.Properties.VariableNames;
    fprintf(fid, '<tr>');
    for i = 1:numel(varNames)
        fprintf(fid, '<th>%s</th>', varNames{i});
    end
    fprintf(fid, '</tr>\n');

    for r = 1:height(report)
        fprintf(fid, '<tr>');
        for c = 1:numel(varNames)
            val = report.(varNames{c})(r);
            if isstring(val) || ischar(val)
                fprintf(fid, '<td style="text-align:left;">%s</td>', string(val));
            else
                if isnan(val)
                    fprintf(fid, '<td>N/A</td>');
                else
                    fprintf(fid, '<td>%g</td>', val);
                end
            end
        end
        fprintf(fid, '</tr>\n');
    end
    fprintf(fid, '</table>\n');

    fprintf(fid, '<p class="note"><strong>How to read the metrics (rule-of-thumb):</strong><br>\n');
    fprintf(fid, 'Vor_CV: Voronoi cell area の変動係数 (std/mean)。小さいほど良い。<br>\n');
    fprintf(fid, 'Vor_MaxMin: Voronoi cell area の最大/最小比。1 に近いほど良い。<br>\n');
    fprintf(fid, 'SH_UI: spherical harmonics に基づく均等性指標。小さいほど良い。<br>\n');
    fprintf(fid, 'NN angles: 最近傍角度分布。極端に小さい角度が少なく、平均・中央値が適度に保たれるほど良い。<br>\n');
    fprintf(fid, 'ES_Emean: electrostatic energy の平均ペア値。小さいほど全体として方向間距離がよく保たれている。</p>\n');

    % Per-shell plots
    fprintf(fid, '<h2>Per-shell plots</h2>\n');

    for i = 1:numel(shell_results)
        s = shell_results(i);

        is_combined = isfield(s, 'is_combined') && s.is_combined;

        if is_combined
            fprintf(fid, '<h3>Combined (all shells; N axis = %g)</h3>\n', s.N_axis);
        else
            fprintf(fid, '<h3>b ≈ %g (N axis = %g)</h3>\n', s.Shell_b, s.N_axis);
        end

        fprintf(fid, '<p class="note">Vor_CV = %.3f, Vor_MaxMin = %.2f, SH_UI = %.3f, NN mean = %.2f°, ES_Emean = %.4f</p>\n', ...
            s.Vor_CV, s.Vor_MaxMin, s.SH_UI, s.NN_mean_deg, s.ES_Emean);

        fprintf(fid, '<div class="panel">\n');

        if ~isempty(s.voronoi_png)
            fprintf(fid, '<figure>\n');
            fprintf(fid, '<img src="%s" alt="Voronoi areas">\n', local_relpath(html_name, s.voronoi_png));
            fprintf(fid, '<figcaption>Voronoi cell area histogram.</figcaption>\n');
            fprintf(fid, '</figure>\n');
        end

        if ~isempty(s.sh_png)
            fprintf(fid, '<figure>\n');
            fprintf(fid, '<img src="%s" alt="SH power spectrum">\n', local_relpath(html_name, s.sh_png));
            fprintf(fid, '<figcaption>SH power spectrum.</figcaption>\n');
            fprintf(fid, '</figure>\n');
        end

        if ~isempty(s.nn_png)
            fprintf(fid, '<figure>\n');
            fprintf(fid, '<img src="%s" alt="Nearest-neighbor angles">\n', local_relpath(html_name, s.nn_png));
            fprintf(fid, '<figcaption>Nearest-neighbor angle distribution.</figcaption>\n');
            fprintf(fid, '</figure>\n');
        end

        fprintf(fid, '</div>\n');
    end

    fprintf(fid, '</body></html>\n');
end

% ============================================================
% サブ関数: HTML から PNG への相対パス
% ============================================================
function rel = local_relpath(~, target_path)
    [~, png_name, png_ext] = fileparts(target_path);
    rel = [png_name png_ext];
end