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

points=1;
n_turns=100;
B_0=1;

phi_span = linspace(0,n_turns*2*pi,n_turns*points+1);


n_line=2000;
n_surf=100;

% n_line=1;
% n_surf=1;

psi_t_min=0.7;
psi_t_max=1.2;
% psi_t_min=1.46000000000000;
% psi_t_min=0.06;
% psi_t_min=1;


speed_n=50;
speed=floor(speed_n/n_turns);
length_ani_line=40;

iota_res=1;
global iota_res_global
iota_res_global=iota_res;

%theta_0 = 2*pi*(rand(n_line,1)-0.5);


flag_store_lines=0;


flag_plot_iota=0;
flag_plot_line_tracing=0;
flag_plot_poincare=1;
flag_plot_phase_tragectories=0;
flag_plot_phase_tragectories_ani=0;

flag_runge_kutta=0;
flag_analyt=1;

m=10;

iota_b=0.7;
iota_a=0.3;

A=-0.004;

prec=1e0;

global prec_global;
global m_global;
global iota_b_global;
global iota_a_global;
global A_global;

prec_global=prec;
m_global=m;
iota_b_global=iota_b;
iota_a_global=iota_a;
A_global=A;

for i=1:n_surf
   alpha_0((i-1)*round(n_line/n_surf)+1:i*round(n_line/n_surf))= 2*pi*(linspace(0.001,0.999,round(n_line/n_surf))-0.5);
end
% for i_line=1:n_line
%    if mod(alpha_0(i_line),pi/m_global)==0
%        alpha_0(i_line)=alpha_0(i_line)+0.0001;
%    end
% end


% alpha_0=zeros(n_line,1);
% %  alpha_0=ones(n_line,1)*pi/5;
% % alpha_0=[2.92493109127326];
% % alpha_0=-pi;
% alpha_0=-3.4879e-16;
%% Culculations
run_time_indepndent_magnetic_hamiltonian

%% Plotting
plot_time_indepndent_magnetic_hamiltonian