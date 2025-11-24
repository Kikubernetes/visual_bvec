function [Lvals, Pl, UI] = sh_uniformity_from_dirs(axis_dirs, Lmax)
% axis_dirs : 3×M or M×3 の単位ベクトル
% Lmax      : 最大次数 (例 8)
%
% Lvals : 0:Lmax
% Pl    : 各次数のパワー
% UI    : Σ_{l>=1} P_l / P_0

    if nargin < 2
        Lmax = 8;
    end

    % 形を M×3 に統一
    if size(axis_dirs,2)==3
        dirs = axis_dirs;
    elseif size(axis_dirs,1)==3
        dirs = axis_dirs.';
    else
        error('axis_dirs must be 3×M or M×3');
    end

    % ゼロ除外 & 正規化
    norms = sqrt(sum(dirs.^2,2));
    dirs  = dirs(norms>0,:);
    norms = norms(norms>0);
    dirs  = dirs ./ norms;

    x = dirs(:,1);
    y = dirs(:,2);
    z = dirs(:,3);

    theta = acos(z);
    phi   = atan2(y,x); phi(phi<0) = phi(phi<0) + 2*pi;

    Lvals = (0:Lmax).';
    Pl    = zeros(size(Lvals));
    ct    = cos(theta(:).');

    for li = 1:numel(Lvals)
        l = Lvals(li);

        Plm_all = legendre(l, ct);   % (m+1)×M
        c_lm = zeros(2*l+1,1);
        idx  = 1;

        for m = -l:l
            if m < 0
                mp  = -m;
                Nlm = sqrt((2*l+1)/(4*pi) * factorial(l-mp)/factorial(l+mp));
                Plm = squeeze(Plm_all(mp+1,:)).';
                Ymp = Nlm * Plm .* exp(1i*mp*phi);
                Ylm = (-1)^mp * conj(Ymp);
                c_lm(idx) = sum(Ylm);
            else
                Nlm = sqrt((2*l+1)/(4*pi) * factorial(l-m)/factorial(l+m));
                Plm = squeeze(Plm_all(m+1,:)).';
                Ylm = Nlm * Plm .* exp(1i*m*phi);
                c_lm(idx) = sum(Ylm);
            end
            idx = idx + 1;
        end

        Pl(li) = sum(abs(c_lm).^2);
    end

    P0 = Pl(1);
    UI = sum(Pl(2:end)) / P0;
end

