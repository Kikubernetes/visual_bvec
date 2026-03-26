function est = electrostatic_energy_from_dirs(axis_dirs, p, bipolar)
% ELECTROSTATIC_ENERGY_FROM_DIRS
%   球面上の方向ベクトルから electrostatic energy を計算する
%
% 入力:
%   axis_dirs : 3xM または Mx3
%   p         : 指数（省略時 1）
%   bipolar   : true なら ±v を両方含めて計算（省略時 true）
%
% 出力 est:
%   .N             : 使用点数
%   .E_total       : 総 energy
%   .E_meanpair    : 平均 pair energy
%   .D             : 距離行列

    if nargin < 2 || isempty(p)
        p = 1;
    end
    if nargin < 3 || isempty(bipolar)
        bipolar = true;
    end

    if size(axis_dirs,1) == 3
        dirs = axis_dirs;
    elseif size(axis_dirs,2) == 3
        dirs = axis_dirs.';
    else
        error('axis_dirs must be 3xM or Mx3');
    end

    % ゼロ除外・正規化
    mask = any(dirs ~= 0, 1);
    dirs = dirs(:, mask);
    dirs = dirs ./ vecnorm(dirs, 2, 1);

    % bipolar 対応
    if bipolar
        dirs = [dirs, -dirs];
    end

    N = size(dirs, 2);

    % 距離行列
    G = dirs.' * dirs;                    % cos に相当
    G = max(min(G, 1), -1);
    D = sqrt(max(0, 2 - 2*G));           % ||vi-vj|| = sqrt(2-2 vi·vj)

    % 対角を除外
    D(1:N+1:end) = Inf;

    % 上三角だけ使って energy 計算
    iu = triu(true(N), 1);
    dvals = D(iu);

    E_total = sum(1 ./ (dvals .^ p));
    E_meanpair = mean(1 ./ (dvals .^ p));

    est.N = N;
    est.E_total = E_total;
    est.E_meanpair = E_meanpair;
    est.D = D;
end