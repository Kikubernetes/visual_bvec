function eval_sh_uniformity_multishell(bvec_file, bval_file, shell_vals)
% EVAL_SH_UNIFORMITY_MULTISHELL
%   bvec, bval のファイル名と、評価したい 2 つの b-shell をユーザーが指定する。
%
% 使い方:
%   eval_sh_uniformity_multishell('DWI_AP.bvec','DWI_AP.bval',[1000 2000]);
%
% 入力:
%   bvec_file : FSL形式の3×N bvecファイル名（例 'DWI_AP.bvec'）
%   bval_file : 1×N bvalファイル名（例 'DWI_AP.bval'）
%   shell_vals: [b1 b2] の 1×2 ベクトル（例 [1000 2000]）
%               → bval が多少ずれていても「近似一致」で抽出する
%
% 出力:
%   各 shell の SH-based Uniformity Index (UI)
%   SH パワースペクトルの図を表示

    % ==== 入力チェック ====
    if numel(shell_vals) ~= 2
        error('shell_vals は [b1 b2] の 1×2 ベクトルで指定してください。');
    end
    b1 = shell_vals(1);
    b2 = shell_vals(2);

    fprintf('Target shells: b ≈ %g and %g\n', b1, b2);

    % ==== bvec/bval 読み込み ====
    bvec = load(bvec_file);   % 3×N を想定
    bval = load(bval_file);   % 1×N を想定

    % ---- bvec: N×3 に統一 ----
    if size(bvec,1)==3 && size(bvec,2)~=3
        bvec = bvec.';   % 3×N → N×3
    end
    if size(bvec,2) ~= 3
        error('bvec の形式が不正です。N×3 または 3×N のはずです。');
    end

    % ---- bval: 1×N に統一 ----
    if size(bval,1) > 1 && size(bval,2)==1
        bval = bval.';   % Nx1 → 1×N
    end
    if size(bval,1) ~= 1
        error('bval の形式が不正です。1×N である必要があります。');
    end

    N = size(bvec,1);
    if length(bval) ~= N
        error('bvec の本数と bval の本数が一致しません。');
    end

    % ==== シェル抽出（近似一致）====
    % bval が ±5％ 程度ずれることを許容
    tol1 = 0.05 * b1;
    tol2 = 0.05 * b2;

    idx1 = abs(bval - b1) <= tol1;
    idx2 = abs(bval - b2) <= tol2;

    fprintf('Shell1 (b≈%g): detected %d dirs\n', b1, nnz(idx1));
    fprintf('Shell2 (b≈%g): detected %d dirs\n', b2, nnz(idx2));

    if nnz(idx1)==0 || nnz(idx2)==0
        error('指定した b-shell の方向が見つかりませんでした。tol を調整してください。');
    end

    dirs1 = bvec(idx1, :);
    dirs2 = bvec(idx2, :);

    % ==== SH による均等性評価 ====
    Lmax = 8;

    [L1, P1, UI1] = sh_uniformity_sphere(dirs1, Lmax);
    [L2, P2, UI2] = sh_uniformity_sphere(dirs2, Lmax);

    fprintf('\n===== Spherical Harmonics Uniformity =====\n');
    fprintf('Shell1  b≈%g : UI = %.4g\n', b1, UI1);
    fprintf('Shell2  b≈%g : UI = %.4g\n', b2, UI2);

    % ==== 図示 ====
    figure;
    subplot(1,2,1);
    stem(L1, P1, 'filled');
    grid on; xlabel('Degree l'); ylabel('Power P_l');
    title(sprintf('Shell1  b≈%g (UI=%.3g)', b1, UI1));

    subplot(1,2,2);
    stem(L2, P2, 'filled');
    grid on; xlabel('Degree l'); ylabel('Power P_l');
    title(sprintf('Shell2  b≈%g (UI=%.3g)', b2, UI2));
end


% ================================================================
% 単一 shell の SH 均等性評価（前回提供したものと同等）
% ================================================================
function [Lvals, Pl, UI] = sh_uniformity_sphere(dirs, Lmax)

    if size(dirs,2)==3
        % OK
    elseif size(dirs,1)==3
        dirs = dirs.';
    else
        error('dirs must be N×3 or 3×N.');
    end

    norms = sqrt(sum(dirs.^2,2));
    dirs  = dirs(norms>0,:);
    norms = norms(norms>0);

    dirs = dirs ./ norms;

    x = dirs(:,1);
    y = dirs(:,2);
    z = dirs(:,3);

    theta = acos(z);
    phi   = atan2(y,x); phi(phi<0)=phi(phi<0)+2*pi;

    Lvals = (0:Lmax).';
    Pl    = zeros(size(Lvals));

    ct = cos(theta(:).');

    for li = 1:numel(Lvals)
        l = Lvals(li);
        Plm_all = legendre(l, ct);

        c_lm = zeros(2*l+1,1);
        idx  = 1;

        for m = -l:l
            if m < 0
                mp = -m;
                Nlm = sqrt((2*l+1)/(4*pi) * factorial(l-mp)/factorial(l+mp));
                Plm = squeeze(Plm_all(mp+1,:)).';
                Ylmp = Nlm * Plm .* exp(1i*mp*phi);
                Ylm  = (-1)^mp * conj(Ylmp);
                c_lm(idx) = sum(Ylm);
            else
                Nlm = sqrt((2*l+1)/(4*pi) * factorial(l-m)/factorial(l+m));
                Plm = squeeze(Plm_all(m+1,:)).';
                Ylm = Nlm * Plm .* exp(1i*m*phi);
                c_lm(idx) = sum(Ylm);
            end
            idx = idx+1;
        end

        Pl(li) = sum(abs(c_lm).^2);
    end

    P0 = Pl(1);
    UI = sum(Pl(2:end)) / P0;
end

