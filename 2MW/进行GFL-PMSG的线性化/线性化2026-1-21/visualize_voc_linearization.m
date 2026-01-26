function visualize_voc_linearization(sys, eigenvalues, P, Q, V, xi)
%VISUALIZE_VOC_LINEARIZATION
%   对 VOC 线性化结果做可视化，包括：
%   1) 特征值在复平面上的分布（极点图）
%   2) 电流输出的小信号阶跃响应（示例）
%   3) IEEE 风格的阻抗伯德图 Z(s) = v_d(s) / i_{gD}(s)
%
%   输入：
%       sys         状态空间模型 ss(A,B,C,D)，来自线性化
%       eigenvalues A 的特征值（列向量）
%       P, Q, V, xi 当前工况（标幺制功率、电压与相角）

    % 为当前工况生成一个易读的标题
    titleStr = sprintf('P = %.2f pu, Q = %.2f pu, V = %.2f pu, \\xi = %.1f^\\circ', ...
                       P, Q, V, rad2deg(xi));

    %% 图 1：极点分布（特征值复平面）
    figure('Name', ['VOC eigenvalues: ', titleStr], ...
           'NumberTitle', 'off', ...
           'Color', 'w', ...
           'Position', [200 200 600 450]);

    plot(real(eigenvalues), imag(eigenvalues), 'x', ...
         'MarkerSize', 8, ...
         'LineWidth', 1.5);
    hold on; grid on; box on;

    xline(0, '--', 'LineWidth', 1.2);
    yline(0, ':',  'LineWidth', 1.0);

    xlabel('Real Part', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel('Imag Part [rad/s]', 'FontSize', 12, 'FontName', 'Times New Roman');
    title(['Eigenvalues of A  |  ', titleStr], ...
          'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    set(gca, 'FontName', 'Times New Roman', ...
             'FontSize', 11, ...
             'LineWidth', 1);

    re = real(eigenvalues);
    imv = imag(eigenvalues);
    if ~isempty(re)
        dx = max(0.5, 0.2 * range(re));
        dy = max(1.0, 0.2 * range(imv));
        xlim([min(re)-dx, max(re)+dx]);
        ylim([min(imv)-dy, max(imv)+dy]);
    end

    %% 图 3：IEEE 风格阻抗伯德图 Z(s) = v_d / i_{gD}
    try
        % 1️⃣ 先从 MIMO 系统中抽取"电压→电流"的那一条通道：
        %     输入：v_d（第 1 个输入）
        %     输出：i_{gD}（第 1 个输出）
        G_iD_vd = sys(1, 1);   % i_{gD}(s) / v_d(s)

        % 2️⃣ 阻抗 = 电压 / 电流 = 1 / 导纳
        %     Z_vd(s) = v_d(s) / i_{gD}(s) = 1 / G_iD_vd(s)
        Z_vd = inv(G_iD_vd);   % SISO 系统，inv 就是 1/G

        % 3️⃣ 选频率范围：根据特征值的虚部自动选范围
        w_imag = abs(imag(eigenvalues));
        if isempty(w_imag) || max(w_imag) == 0
            w_max = 1e4;      % rad/s
        else
            w_max = 10 * max(w_imag);
        end
        w_min = w_max / 1e4;
        w_min = max(w_min, 1e-1);
        w_max = max(w_max, 1);

        Nw = 400;
        w   = logspace(log10(w_min), log10(w_max), Nw);  % rad/s
        [mag, phase] = bode(Z_vd, w);   % mag: |Z(jw)|, phase: ∠Z(jw) (deg)

        mag   = squeeze(mag);
        phase = squeeze(phase);
        f_Hz  = w / (2*pi);             % rad/s -> Hz

        figure('Name', ['VOC impedance Bode: ', titleStr], ...
               'NumberTitle', 'off', ...
               'Color', 'w', ...
               'Position', [300 100 700 500]);

        % 上：幅频（阻抗幅值，单位：dB）
        subplot(2,1,1);
        semilogx(f_Hz, 20*log10(mag), 'LineWidth', 1.5);
        grid on; box on;
        ylabel('|Z(j\omega)| (dB)', 'FontSize', 12, 'FontName', 'Times New Roman');
        title(['Bode of Z_{dd}(s) = v_{d}(s) / i_{gD}(s)  |  ', titleStr], ...
              'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
        set(gca, 'FontName', 'Times New Roman', ...
                 'FontSize', 11, ...
                 'LineWidth', 1);

        % 下：相频
        subplot(2,1,2);
        semilogx(f_Hz, phase, 'LineWidth', 1.5);
        grid on; box on;
        xlabel('Frequency (Hz)', 'FontSize', 12, 'FontName', 'Times New Roman');
        ylabel('Phase (deg)',   'FontSize', 12, 'FontName', 'Times New Roman');
        set(gca, 'FontName', 'Times New Roman', ...
                 'FontSize', 11, ...
                 'LineWidth', 1);

    catch ME
        warning('绘制阻抗伯德图失败：%s', ME.message);
    end

end
