function report = uniformity_report_one_dataset( ...
    dataset_name, bvec_file, bval_file, target_b_list, tol, Lmax, pdf_name)
% UNIFORMITY_REPORT_ONE_DATASET
%   1つのデータセットについて、球面ボロノイ・SH・最近傍角度による
%   均等性指標をまとめて table にし、必要なら PDF に図を出力する。
%
% 入力:
%   dataset_name  : データセット名 (例 'dir98_AP')
%   bvec_file     : bvec ファイル名 (例 'dMRI_dir98_AP.bvec')
%   bval_file     : bval ファイル名 (例 'dMRI_dir98_AP.bval')
%   target_b_list : 評価したい b 値 (例 [1500 3000])
%   tol           : b 値許容幅 (省略時 50)
%   Lmax          : SH の最大次数 (省略時 8)
%   pdf_name      : 図を出力する PDF 名 (例 'dir98_uniformity.pdf')
%                   空文字や省略時は PDF 出力なし
%
% 出力:
%   report : table
%     列: Dataset, Shell_b, N_axis, Vor_CV, Vor_MaxMin,
%         SH_UI, NN_min_deg, NN_mean_deg, NN_median_deg, NN_max_deg

    % --- オプション引数のデフォルト ---
    if nargin < 5 || isempty(tol)
        tol = 50;
    end
    if nargin < 6 || isempty(Lmax)
        Lmax = 8;
    end
    if nargin < 7
        pdf_name = '';
    end

    target_b_list = target_b_list(:).';  % 行ベクトルに揃える

    rows = [];              % struct の配列として蓄積
    created_pdf = false;    % PDF への初回書き込みフラグ

    % ==========================
    % 各 shell ごとの評価ループ
    % ==========================
    for ib = 1:numel(target_b_list)
        tb = target_b_list(ib);

          % 1) 球面ボロノイ（軸対称 v/-v を考慮）
        vs = voronoi_shell_bipolar(bvec_file, bval_file, tb, tol);
        dirs_axis = vs.axis_dirs;  % 3×M の単位ベクトル

        % 2) SH による均等性 (UI)
        [Lvals, Pl, UI] = sh_uniformity_from_dirs(dirs_axis, Lmax);

          % 3) 最近傍角度分布
        nn = nearest_angle_from_dirs(dirs_axis);

        % 4) table 用の 1 行を struct で作成
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

        rows = [rows; row]; %#ok<AGROW>

        % 5) PDF 出力（指定がある場合）: 各 shell の図を追加
        if ~isempty(pdf_name)
            % (a) Voronoi 面積のヒストグラム
            fig1 = figure('Visible','off');
            histogram(vs.vertex_area, 20);
            xlabel('Voronoi cell area');
            ylabel('Count');
            title(sprintf('%s  b≈%g  Voronoi areas (CV=%.3f)', ...
                  dataset_name, tb, vs.CV));
            drawnow;

            if ~created_pdf
                exportgraphics(fig1, pdf_name, 'Resolution', 300);
                created_pdf = true;
            else
                exportgraphics(fig1, pdf_name, 'Resolution', 300, 'Append', true);
            end
            close(fig1);

            % (b) SH パワースペクトル
            fig2 = figure('Visible','off');
            stem(Lvals, Pl, 'filled');
            xlabel('Degree l');
            ylabel('Power P_l');
            title(sprintf('%s  b≈%g  SH power (UI=%.3f)', ...
                  dataset_name, tb, UI));
            grid on;
            drawnow;
            exportgraphics(fig2, pdf_name, 'Resolution', 300, 'Append', true);
            close(fig2);

            % (c) 最近傍角度分布
            fig3 = figure('Visible','off');
            histogram(nn.theta, 0:2:60);
            xlabel('Nearest-neighbor angle (deg)');
            ylabel('Count');
            title(sprintf('%s  b≈%g  NN angles (mean=%.2f°)', ...
                  dataset_name, tb, nn.mean));
            grid on;
            drawnow;
            exportgraphics(fig3, pdf_name, 'Resolution', 300, 'Append', true);
            close(fig3);
        end
    end

    % ==========================
    % 6) struct 配列 → table へ変換
    % ==========================
    report = struct2table(rows);

    % 7) PDF に表ページを追加（あれば）
    if ~isempty(pdf_name)
        txt = evalc('disp(report)');   % table のテキスト表現

        fig0 = figure('Visible','off');
        axis off
        text(0, 1, txt, ...
            'FontName','Courier', ...   % 等幅フォント
            'FontSize',10, ...
            'Interpreter','none', ...
            'VerticalAlignment','top');

        if ~created_pdf
            exportgraphics(fig0, pdf_name, 'Resolution',300);
        else
            exportgraphics(fig0, pdf_name, 'Resolution',300,'Append',true);
        end
        close(fig0);
    end

    % 画面にも軽く表示
    disp(report);
end

