%% Preparation step
clear all

set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultFigureWindowStyle','docked')
format compact

Color_matlab={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% global Color_matlab
% 
% global SYMBOL;
% global COLOUR;
SYMBOL={'-o' '-s' '-*' '-h' '-^' '->' '-<' '-v' '-+' '-o' '-s' '-*' '-h' '-^' '->' '-<' '-v' '-+' '-o' '-s' '-*' '-h' '-^' '->' '-<' '-v' '-o' '-s' '-*' '-h' '-^' '->' '-<' '-v'};
COLOUR={[1,0,0] [0,0,1] [0.45,0.1,0.05] [0.75,0.2,0.15] [0.15,0.05,0.5] [0.25,0.5,0.25] [1,0.6,0] [1,0.05,0.75] ...
    [0.5,0.5,1] [0.0,0.5,1] [0.5,0.1,0] [0.6,0.05,0.4] [0.35,0.65,0.15] [0.8,0.2,0.8] [0.2,0.2,0.5] [0.4,0.4,0.55] ...
    [0.5,0.5,0.5] [0.2,0.8,0.8] [0.8,0.7,0.3] [0.3,0.9,0.5] [0,1,0.8] [0.25,0.25,0.85] [0.85,0.2,0.2] [0.3,0,0.9] [0.75,0.75,0.05] ...
    [0.9,0.2,0.2] [0.2,0.9,0.9] [0.15,0.9,0.15] [0.2,0.15,0.9] [0.2,0.2,0.65]};

%% input parameters
R_0=1.75;
% R=linspace(0.5,1.0,50);
% Z=linspace(-0.25,0.25,50);
% 
% r=linspace(0,2,50);
% theta=linspace(-pi,pi,50);
% 
% phi = linspace(0,2*pi,50);

points=16000;
n_turns=100;
B_0=1;

phi_span = linspace(0,n_turns*2*pi,n_turns*points);


n_line=100;
n_surf=4;

% n_line=1;
% n_surf=1;

psi_t_min=0.5;
psi_t_max=1.2;

%theta_0 = 2*pi*(rand(n_line,1)-0.5);
for i=1:n_surf
   theta_0((i-1)*round(n_line/n_surf)+1:i*round(n_line/n_surf))= 2*pi*(linspace(0.001,0.999,round(n_line/n_surf))-0.5);
end

flag_store_lines=0;


flag_plot_iota=0;
flag_plot_line_tracing=0;
flag_plot_poincare=1;


%% Culculations
%preparation run
i_line=1;
psi_t_0(i_line) = (psi_t_max-psi_t_min)*1/round(n_line/n_surf)*floor((i_line-1)/n_surf)+psi_t_min;
y0 = [psi_t_0(i_line); theta_0(i_line)];

% Solver
[ phi_out, yout ] = ode45( @fun, phi_span, y0);

for i_turns=0:n_turns
    temp=find(phi_out>=i_turns*2*pi);
    I_phi(i_turns+1)=temp(1);    
end
n_I_phi=length(I_phi);

%production run
for i_line=1:n_line
    psi_t_0(i_line) = (psi_t_max-psi_t_min)*1/round(n_line/n_surf)*floor((i_line-1)/n_surf)+psi_t_min;
    y0 = [psi_t_0(i_line); theta_0(i_line)];
    
    % Solver
    [ phi_out, yout ] = ode45( @fun, phi_span, y0);
    
    psi_tor_out(:,1)=yout(:,1);
    theta_out(:,1)=yout(:,2);
    r_out(:,1)=sqrt(psi_tor_out(:,1)./(B_0*pi));
    
    R_out(:,1)=R_0+r_out(:,1).*cos(theta_out(:,1));
    Z_out(:,1)=-r_out(:,1).*sin(theta_out(:,1));
    
    [X_out(:,1),Y_out(:,1),Z_out(:,1)]=pol2cart(phi_out(:,1),R_out(:,1),Z_out(:,1));
    
    
    if flag_store_lines
        theta_out_line{i_line}=theta_out;
        R_out_line{i_line}=R_out;
        X_out_line{i_line}=X_out;
        Y_out_line{i_line}=Y_out;
        Z_out_line{i_line}=Z_out;
    end

    R_out_poicare{i_line}(1:n_I_phi)=R_out(I_phi,1);
    Z_out_poicare{i_line}(1:n_I_phi)=Z_out(I_phi,1);
    
    iota_out{i_line}=(theta_out(I_phi(2),1)-theta_out(I_phi(1),1))./(2*pi);
end
%%

global iota_0_global;
global iota_a_global;

psi_t_test=linspace(min(psi_t_0),max(psi_t_0),20);
% iota_test=iota_a_global.*psi_t_test.^2+iota_0_global;
iota_test=iota_a_global.*psi_t_test+iota_0_global;

iota_res=[10/10, 10/11];
points_res=100;
R_iota_res=zeros(length(iota_res),points_res);
Z_iota_res=zeros(length(iota_res),points_res);
for i=1:length(iota_res)
    psi_t_iota_res(i)=(iota_res(i)-iota_0_global)./iota_a_global;
    r_iota_res(i)=sqrt(psi_t_iota_res(i)./(B_0*pi));
    [R_iota_res(i,:),Z_iota_res(i,:)]=pol2cart(linspace(0,2*pi,points_res),r_iota_res(i));
    R_iota_res(i,:)=R_iota_res(i,:)+R_0;
end
%% Plotting
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
if flag_plot_poincare
    figure;
    axes('FontName','Times','FontSize',24);
    Xlabel='R';
    Ylabel='Z';
    Gtitle=['Runge$-$Kutta, poincare plot, $\phi$ resolution: ', num2str(points)];
    
    
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
%% Function for ode45

function dydt = fun (t, y)
iota_0=0.7;
iota_a=0.3;

psi_power=5;
gauss_mean=1.000;
gauss_div=0.04;

n1=10;
m1=10;
A1=-0.004;

n2=10;
m2=11;
A2=-0.000;

global iota_0_global;
global iota_a_global;

iota_0_global=iota_0;
iota_a_global=iota_a;

rhs1=0;
%rhs1=A.*sin(y(2)-1./(iota_res).*t);
rhs1=A1.*m1.*sin(m1.*y(2)-n1.*t)+A2.*m2.*sin(m2.*y(2)-n2.*t);
%rhs1=A1.*y(1).^psi_power.*sin(m1.*y(2)-n1.*t)+A2.*y(1).^psi_power.*sin(m2.*y(2)-n2.*t);
%rhs1=A1./(gauss_div.*sqrt(2*pi)).*exp(-(y(1)-gauss_mean).^2./(2*gauss_div.^2)).*sin(m1.*y(2)-n1.*t);

% rhs2=((iota_a.*y(1).^2+iota_0)+2*(iota_a.*y(1).^2));
rhs2=iota_a.*y(1)+iota_0;
%rhs2=iota_a.*y(1)+iota_0-A1.*(y(1)-gauss_mean)./(gauss_div.^3.*sqrt(2*pi)).*exp(-(y(1)-gauss_mean).^2./(2*gauss_div.^2)).*cos(m1.*y(2)-n1.*t);
%rhs=((iota_a.*y(1).^2+iota_0)+2*(iota_a.*y(1).^2))-2*iota_a.*y(1).*t/(iota_a.*y(1).^2+iota_0).^2.*sin(y(2)-t./(iota_a.*y(1).^2+iota_0));

dydt = [rhs1; rhs2];
end