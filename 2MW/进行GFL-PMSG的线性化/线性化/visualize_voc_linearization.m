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

    % %% 图 1：极点分布（特征值复平面）
    % figure('Name', ['VOC eigenvalues: ', titleStr], ...
    %        'NumberTitle', 'off', ...
    %        'Color', 'w', ...
    %        'Position', [200 200 600 450]);
    % 
    % plot(real(eigenvalues), imag(eigenvalues), 'x', ...
    %      'MarkerSize', 8, ...
    %      'LineWidth', 1.5);
    % hold on; grid on; box on;
    % 
    % xline(0, '--', 'LineWidth', 1.2);
    % yline(0, ':',  'LineWidth', 1.0);
    % 
    % xlabel('Real Part', 'FontSize', 12, 'FontName', 'Times New Roman');
    % ylabel('Imag Part [rad/s]', 'FontSize', 12, 'FontName', 'Times New Roman');
    % title(['Eigenvalues of A  |  ', titleStr], ...
    %       'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    % 
    % set(gca, 'FontName', 'Times New Roman', ...
    %          'FontSize', 11, ...
    %          'LineWidth', 1);
    % 
    % re = real(eigenvalues);
    % imv = imag(eigenvalues);
    % if ~isempty(re)
    %     dx = max(0.5, 0.2 * range(re));
    %     dy = max(1.0, 0.2 * range(imv));
    %     xlim([min(re)-dx, max(re)+dx]);
    %     ylim([min(imv)-dy, max(imv)+dy]);
    % end

    %% 图 3：IEEE 风格阻抗伯德图 Z(s) = v_d / i_{gD}
    try
        % ================== 📌 可手动调节的坐标范围区域 ==================
        % 频率范围（Hz，幅值 & 相位共用）
        f_min_plot = 1;          % 10^0 Hz
        f_max_plot = 1e3;        % 10^3 Hz

        % 幅值纵轴范围（dB）：[] 表示自适应
        mag_ylim   = [];         % 例如 [20 80]；留空 [] = 自动

        % 相位纵轴范围（deg）：[] 表示自适应  ✅ 修改1：注释+默认值统一为[]
        phase_ylim = [];         % 例如 [-200 200]；留空 [] = 自动（原固定值：[-200 200]）

        % 频率采样点个数
        Nw = 400;
        % =========================================================

        % 1️⃣ 从 MIMO 系统中抽取"电压→电流"的那一条通道：
        %     输入：v_d（第 1 个输入）
        %     输出：i_{gD}（第 1 个输出）
        G_iD_vd = sys(2, 2);   % i_{gD}(s) / v_d(s)

        % 2️⃣ 阻抗 = 电压 / 电流 = 1 / 导纳
        %     Z_vd(s) = v_d(s) / i_{gD}(s) = 1 / G_iD_vd(s)
        Z_vd = inv(G_iD_vd);   % SISO 系统，inv 就是 1/G

        % 3️⃣ 频率向量（采用可调范围）
        f_Hz = logspace(log10(f_min_plot), log10(f_max_plot), Nw);  % Hz
        w    = 2*pi*f_Hz;                                           % rad/s

        [mag, phase] = bode(Z_vd, w);   % mag: |Z(jw)|, phase: ∠Z(jw) (deg)
        mag   = squeeze(mag);
        phase = squeeze(phase);

        % 将相位折叠到 [-180, 180] 区间内
        % phase = mod(phase + 180, 360) - 180;


        % 4️⃣ 绘图（稍微美化一下）
        figure('Name', ['VOC impedance Bode: ', titleStr], ...
            'NumberTitle', 'off', ...
            'Color', 'w', ...
            'Position', [300 100 800 550]);

        % ------- 上：幅频（阻抗幅值，单位：dB） -------
        subplot(2,1,1);
        semilogx(f_Hz, 20*log10(mag), 'LineWidth', 1.8);
        grid on; box on;
        xlim([f_min_plot, f_max_plot]);    % 横轴范围 10^0 ~ 10^3
        if ~isempty(mag_ylim)
            ylim(mag_ylim);                % 如设置了 mag_ylim，则应用
        end
        ylabel('|Z(j\omega)| (dB)', 'FontSize', 12, 'FontName', 'Times New Roman');
        title(['Bode of Z_{dd}(s) = v_{d}(s) / i_{gD}(s)  |  ', titleStr], ...
            'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
        set(gca, 'FontName', 'Times New Roman', ...
            'FontSize', 11, ...
            'LineWidth', 1, ...
            'XMinorTick', 'on', ...
            'YMinorTick', 'on');

        % ------- 下：相频 -------
        subplot(2,1,2);
        semilogx(f_Hz, phase, 'LineWidth', 1.8);
        grid on; box on;
        xlim([f_min_plot, f_max_plot]);    % 横轴范围 10^0 ~ 10^3
        % ✅ 修改2：相位纵轴添加非空判断（和幅值逻辑一致）
        if ~isempty(phase_ylim)
            ylim(phase_ylim);              % 仅当设置了值时，才固定相位范围
        end
        xlabel('Frequency (Hz)', 'FontSize', 12, 'FontName', 'Times New Roman');
        ylabel('Phase (deg)',   'FontSize', 12, 'FontName', 'Times New Roman');
        set(gca, 'FontName', 'Times New Roman', ...
            'FontSize', 11, ...
            'LineWidth', 1, ...
            'XMinorTick', 'on', ...
            'YMinorTick', 'on');

    catch ME
        warning('绘制阻抗伯德图失败：%s', ME.message);
    end


end




% function visualize_voc_linearization(sys, eigenvalues, P, Q, V, xi)
% %VISUALIZE_VOC_LINEARIZATION
% %   基于 IEEE 论文常见画图风格，绘制：
% %   1) 特征值在复平面上的分布（极点图）
% %   2) 阻抗伯德图 Z(s) = v_d(s) / i_{gD}(s)
% %
% %   输入：
% %       sys         状态空间模型 ss(A,B,C,D)，来自线性化
% %       eigenvalues A 的特征值（列向量）
% %       P, Q, V, xi 当前工况（标幺制功率、电压与相角）
% 
%     % 当前工况标题
%     titleStr = sprintf('P = %.2f pu, Q = %.2f pu, V = %.2f pu, \\xi = %.1f^\\circ', ...
%                        P, Q, V, rad2deg(xi));
% 
%     %% ---------------- 图 1：极点分布（特征值复平面） ----------------
%     figure;
%     re  = real(eigenvalues(:));
%     imv = imag(eigenvalues(:));
% 
%     plot(re, imv, 'rx', 'LineWidth', 1.2, 'MarkerSize', 7);  % 黑色叉号，常见 IEEE 风格
%     hold on;
%     xline(0, 'k--', 'LineWidth', 1.0);   % 实轴、虚轴细虚线
%     yline(0, 'k--', 'LineWidth', 1.0);
%     hold off;
% 
%     grid on;
%     box on;
% 
%     xlabel('Real Part', 'FontName', 'Times New Roman', 'FontSize', 11);
%     ylabel('Imag Part [rad/s]', 'FontName', 'Times New Roman', 'FontSize', 11);
%     title(['Eigenvalues of A  |  ', titleStr], ...
%           'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'normal');
% 
%     set(gca, 'FontName', 'Times New Roman', ...
%              'FontSize', 10, ...
%              'LineWidth', 0.8);
% 
%     % 稍微放宽坐标范围，便于打印 / 观测
%     if ~isempty(re)
%         dx = max(0.5, 0.2 * range(re));
%         dy = max(1.0, 0.2 * range(imv));
%         if dx == 0, dx = 1; end
%         if dy == 0, dy = 5; end
%         xlim([min(re)-dx, max(re)+dx]);
%         ylim([min(imv)-dy, max(imv)+dy]);
%     end
% 
%     %% ---------------- 图 2：阻抗伯德图 Z(s) = v_d / i_{gD} ----------------
%     % IEEE 风格：双子图，频率 Hz（对数坐标），幅值 dB，相位 deg
%     try
%         % 1) 取通道：i_{gD} <- v_d
%         G_iD_vd = sys(1, 1);      % i_{gD}(s) / v_d(s)
% 
%         % 2) Z(s) = v_d(s) / i_{gD}(s) = 1 / G_iD_vd(s)
%         Z_vd = inv(G_iD_vd);      % SISO 时等效为 1/G
% 
%         % 3) 根据特征值虚部选频率范围
%         w_imag = abs(imag(eigenvalues(:)));
%         if isempty(w_imag) || max(w_imag) == 0
%             w_max = 1e4;          % rad/s
%         else
%             w_max = 10 * max(w_imag);
%         end
%         w_min = w_max / 1e4;
%         w_min = max(w_min, 1e-1);
%         w_max = max(w_max, 1);
% 
%         Nw   = 600;
%         w    = logspace(log10(w_min), log10(w_max), Nw);  % rad/s
%         [mag, phase] = bode(Z_vd, w);    % mag: |Z(jw)|，phase: deg
% 
%         mag   = squeeze(mag);
%         phase = squeeze(phase);
%         f_Hz  = w / (2*pi);              % rad/s -> Hz
% 
%         figure;
% 
%         % ---- 上子图：幅频特性（dB） ----
%         subplot(2,1,1);
%         semilogx(f_Hz, 20*log10(mag), 'b-', 'LineWidth', 1.2);
%         grid on; box on;
%         ylabel('|Z(j\omega)| (dB)', 'FontName', 'Times New Roman', 'FontSize', 11);
%         set(gca, 'FontName', 'Times New Roman', ...
%                  'FontSize', 10, ...
%                  'LineWidth', 0.8, ...
%                  'XTickLabel', []);      % 上图不显示 x 轴刻度标签
% 
%         % ---- 下子图：相频特性（deg） ----
%         subplot(2,1,2);
%         semilogx(f_Hz, phase, 'b-', 'LineWidth', 1.2);
%         grid on; box on;
%         xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 11);
%         ylabel('Phase (deg)',    'FontName', 'Times New Roman', 'FontSize', 11);
%         set(gca, 'FontName', 'Times New Roman', ...
%                  'FontSize', 10, ...
%                  'LineWidth', 0.8);
% 
%         % 整体标题（放在 figure 顶部）
%         sgtitle(['Bode of Z_{dd}(s) = v_{d}(s) / i_{gD}(s)  |  ', titleStr], ...
%                 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'normal');
% 
%     catch ME
%         warning('绘制阻抗伯德图失败：%s', ME.message);
%     end
% end
