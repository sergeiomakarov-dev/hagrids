if flag_plot_iota
    figure;
    axes('FontName','Times','FontSize',24);
    Xlabel='$\Psi_t$';
    Ylabel='$\iota$';
    Gtitle='Iota vs toroidal magnetic flux';
    
    hold on
    plot(psi_t_test,iota_test,'LineWidth',3,'Color',Color_matlab{1});
    for i_line=1:n_line
        plot(psi_t_0(i_line),iota_out{i_line},SYMBOL{2},'LineWidth',3,'Color',Color_matlab{2});
    end
    set(gca,'FontName','Times','FontSize',36);
    
    xlabel(Xlabel,'interpreter','latex','FontName','Times','FontSize',36);
    ylabel(Ylabel,'interpreter','latex','FontName','Times','FontSize',36);
    title(Gtitle,'interpreter','latex','FontName','Times','FontSize',36);
    grid on;
    grid minor;
    set(gca,'LineWidth',3);
end

% figure;
% axes('FontName','Times','FontSize',24);
% Xlabel='$\Psi_t$';
% Ylabel='$\iota$';
% Gtitle='Iota vs toroidal magnetic flux';
% 
% hold on
% plot(phi_out,yout(:,1),SYMBOL{1},'LineWidth',3,'Color',Color_matlab{1});
% plot(phi_out,yout(:,2),SYMBOL{1},'LineWidth',3,'Color',Color_matlab{2});
% set(gca,'FontName','Times','FontSize',36);
% 
% xlabel(Xlabel,'interpreter','latex','FontName','Times','FontSize',36);
% ylabel(Ylabel,'interpreter','latex','FontName','Times','FontSize',36);
% title(Gtitle,'interpreter','latex','FontName','Times','FontSize',36);
% grid on;
% grid minor;
% set(gca,'LineWidth',3);
%%
if flag_plot_line_tracing && flag_store_lines
    figure;
    axes('FontName','Times','FontSize',24);
    Xlabel='X';
    Ylabel='Y';
    Zlabel='Z';
    Gtitle='Line tracing';
    
    i_line=1;
    
    plot3(X_out_line{i_line},Y_out_line{i_line},Z_out_line{i_line},'LineWidth',1);
    hold on
    for i_line=2:n_line
        plot3(X_out_line{i_line},Y_out_line{i_line},Z_out_line{i_line},'LineWidth',1);
    end
    set(gca,'FontName','Times','FontSize',36);
    axis(gca,'equal');
    
    xlabel(Xlabel,'interpreter','latex','FontName','Times','FontSize',36);
    ylabel(Ylabel,'interpreter','latex','FontName','Times','FontSize',36);
    zlabel(Zlabel,'interpreter','latex','FontName','Times','FontSize',36);
    title(Gtitle,'interpreter','latex','FontName','Times','FontSize',36);
    grid on;
    grid minor;
    set(gca,'LineWidth',3);
end
%%
if flag_plot_phase_tragectories && flag_store_lines
    figure;
    axes('FontName','Times','FontSize',24);
    Xlabel='$\gamma=m\alpha$';
    Ylabel='$p_\gamma=\sqrt{\iota_a}(\psi-\psi_0)$';
    Gtitle='Phase portarait';
    
    if flag_analyt
        i_line=1;
        plot(m_global*alpha_out_line{i_line},sqrt(iota_b_global)*(psi_tor_out_line{i_line}-psi_t_iota_res(1)),'LineWidth',1);
        hold on
        for i_line=2:n_line
            plot(m_global*alpha_out_line{i_line},sqrt(iota_b_global)*(psi_tor_out_line{i_line}-psi_t_iota_res(1)),'LineWidth',1);
        end
    end
    
    if flag_runge_kutta
        i_line=1;
        plot(m_global*alpha_out_line_num{i_line},sqrt(iota_b_global)*(psi_tor_out_line_num{i_line}-psi_t_iota_res(1)),'LineWidth',1);
        hold on
        for i_line=2:n_line
            plot(m_global*alpha_out_line_num{i_line},sqrt(iota_b_global)*(psi_tor_out_line_num{i_line}-psi_t_iota_res(1)),'LineWidth',1);
        end
    end
    set(gca,'FontName','Times','FontSize',36);
    
    xlabel(Xlabel,'interpreter','latex','FontName','Times','FontSize',36);
    ylabel(Ylabel,'interpreter','latex','FontName','Times','FontSize',36);
    title(Gtitle,'interpreter','latex','FontName','Times','FontSize',36);
    grid on;
    grid minor;
    set(gca,'LineWidth',3);
end
%%
if flag_plot_phase_tragectories_ani && flag_store_lines
    figure;
    axes('FontName','Times','FontSize',24);
    Xlabel='$\gamma=m\alpha$';
    Ylabel='$p_\gamma=\sqrt{\iota_a}(\psi-\psi_0)$';
    Gtitle='Phase portarait';
    
    i_line=1;
    plot(m_global*alpha_out_line{i_line},sqrt(iota_b_global)*(psi_tor_out_line{i_line}-psi_t_iota_res(1)),'LineWidth',1,'Color',Color_matlab{1});
    hold on

    
    %comet(m_global*alpha_out_line{i_line},sqrt(iota_b_global)*(psi_tor_out_line{i_line}-psi_t_iota_res(1)));
    set(gca,'FontName','Times','FontSize',36);
    
    xlabel(Xlabel,'interpreter','latex','FontName','Times','FontSize',36);
    ylabel(Ylabel,'interpreter','latex','FontName','Times','FontSize',36);
    title(Gtitle,'interpreter','latex','FontName','Times','FontSize',36);
    grid on;
    grid minor;
    set(gca,'LineWidth',3);
    for j=1:floor(length(alpha_out_line{i_line})./length_ani_line)
        h = animatedline('Color',Color_matlab{mod(j,size(Color_matlab,2))+1},'LineWidth',2);
        for i = (j-1)*length_ani_line+1:j*length_ani_line
            addpoints(h,m_global*alpha_out_line{i_line}(i),sqrt(iota_b_global)*(psi_tor_out_line{i_line}(i)-psi_t_iota_res(1)));
            if mod(i,speed)==0
                drawnow
            end
        end
    end
end

%%
if flag_plot_poincare
    figure;
    axes('FontName','Times','FontSize',24);
    Xlabel='R';
    Ylabel='Z';
    if flag_analyt
        Gtitle=['Analytical, poincare plot, $\phi$ resolution: ', num2str(points)];
    elseif flag_runge_kutta
        Gtitle=['Runge–Kutta, poincare plot, $\phi$ resolution: ', num2str(points)];
    else
        Gtitle=['Poincare plot, $\phi$ resolution: ', num2str(points)];
    end
    
    
    i_line=1;
    hold on
    plot(R_out_poicare{i_line}(1:n_I_phi),Z_out_poicare{i_line}(1:n_I_phi),'.','Color',Color_matlab{2});
    for i_line=2:n_line
        plot(R_out_poicare{i_line}(1:n_I_phi),Z_out_poicare{i_line}(1:n_I_phi),'.','Color',Color_matlab{2});
    end
    for i=1:length(iota_res)
        plot(R_iota_res(i,:),Z_iota_res(i,:),'-','LineWidth',2,'Color',Color_matlab{1});
    end
    
    set(gca,'FontName','Times','FontSize',36);
    axis(gca,'equal');
    
    xlabel(Xlabel,'interpreter','latex','FontName','Times','FontSize',36);
    ylabel(Ylabel,'interpreter','latex','FontName','Times','FontSize',36);
    title(Gtitle,'interpreter','latex','FontName','Times','FontSize',36);
    grid on;
    grid minor;
    set(gca,'LineWidth',3);
    ylim([-0.7 0])
    xlim([1.6 2.5])
end
