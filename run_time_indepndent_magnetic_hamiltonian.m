

y0=[0; 0];
[dummy1, dummy2] = ode45( @fun, phi_span(1:2), y0);

phi_out=phi_span';

for i_turns=0:n_turns
    temp=find(phi_out>=i_turns*2*pi);
    I_phi(i_turns+1)=temp(1);    
end
n_I_phi=length(I_phi);

if flag_runge_kutta
    %production run
    for i_line=1:n_line
        psi_t_0(i_line) = (psi_t_max-psi_t_min)*1/round(n_line/n_surf)*floor((i_line-1)/n_surf)+psi_t_min;
        y0 = [psi_t_0(i_line); alpha_0(i_line)];
        
        % Solver
        [ phi_out, yout ] = ode45( @fun, phi_span, y0);
        
        psi_tor_out(:,1)=yout(:,1);
        alpha_out(:,1)=yout(:,2);
        theta_out(:,1)=alpha_out(:,1)+iota_res*phi_out(:,1);
        r_out(:,1)=sqrt(psi_tor_out(:,1)./(B_0*pi));
        
        R_out(:,1)=R_0+r_out(:,1).*cos(theta_out(:,1));
        Z_out(:,1)=-r_out(:,1).*sin(theta_out(:,1));
        
        [X_out(:,1),Y_out(:,1),Z_out(:,1)]=pol2cart(phi_out(:,1),R_out(:,1),Z_out(:,1));
        
        
        if flag_store_lines
            psi_tor_out_line_num{i_line}=psi_tor_out;
            alpha_out_line_num{i_line}=alpha_out;
            theta_out_line_num{i_line}=theta_out;
            R_out_line_num{i_line}=R_out;
            X_out_line_num{i_line}=X_out;
            Y_out_line_num{i_line}=Y_out;
            Z_out_line_num{i_line}=Z_out;
        end
        
        R_out_poicare{i_line}(1:n_I_phi)=R_out(I_phi,1);
        Z_out_poicare{i_line}(1:n_I_phi)=Z_out(I_phi,1);
        
        iota_out{i_line}=(theta_out(I_phi(2),1)-theta_out(I_phi(1),1))./(2*pi);
    end
end

%analytical
if flag_analyt
    for i_line=1:n_line
        psi_t_0(i_line) = (psi_t_max-psi_t_min)*1/round(n_line/n_surf)*floor((i_line-1)/n_surf)+psi_t_min;
        y0 = [psi_t_0(i_line); alpha_0(i_line)];

        % Solver
        [yout] = pendulum(phi_span, y0);

        psi_tor_out(:,1)=yout(:,1);
        alpha_out(:,1)=yout(:,2);
        theta_out(:,1)=alpha_out(:,1)+iota_res*phi_out(:,1);
        r_out(:,1)=sqrt(psi_tor_out(:,1)./(B_0*pi));

        R_out(:,1)=R_0+r_out(:,1).*cos(theta_out(:,1));
        Z_out(:,1)=-r_out(:,1).*sin(theta_out(:,1));

        [X_out(:,1),Y_out(:,1),Z_out(:,1)]=pol2cart(phi_out(:,1),R_out(:,1),Z_out(:,1));


        if flag_store_lines
            psi_tor_out_line{i_line}=psi_tor_out;
            alpha_out_line{i_line}=alpha_out;
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
end
%%

global m_global;
global iota_b_global;
global iota_a_global;

psi_t_test=linspace(min(psi_t_0),max(psi_t_0),20);
% iota_test=iota_a_global.*psi_t_test.^2+iota_b_global;
iota_test=iota_a_global.*psi_t_test+iota_b_global;


points_res=100;
R_iota_res=zeros(length(iota_res),points_res);
Z_iota_res=zeros(length(iota_res),points_res);
for i=1:length(iota_res)
    psi_t_iota_res(i)=(iota_res(i)-iota_b_global)./iota_a_global;
    r_iota_res(i)=sqrt(psi_t_iota_res(i)./(B_0*pi));
    [R_iota_res(i,:),Z_iota_res(i,:)]=pol2cart(linspace(0,2*pi,points_res),r_iota_res(i));
    R_iota_res(i,:)=R_iota_res(i,:)+R_0;
end

%% Function for ode45

function dydt = fun (t, y)

global m_global;
global iota_b_global;
global iota_a_global;
global A_global;

iota_b=iota_b_global;
iota_a=iota_a_global;
m=m_global;
A=A_global;

% rhs1=0;
%rhs1=A.*sin(y(2)-1./(iota_res).*t);
%rhs1=A1.*sin(m1.*y(2)-n1.*t)+A2./(y(1)).*sin(m2.*y(2)-n2.*t);
%rhs1=A1.*y(1).^psi_power.*sin(m1.*y(2)-n1.*t)+A2.*y(1).^psi_power.*sin(m2.*y(2)-n2.*t);
% rhs1=A1./(gauss_div.*sqrt(2*pi)).*exp(-(y(1)-gauss_mean1).^2./(2*gauss_div.^2)).*sin(m1.*y(2)-n1.*t)+A2./(gauss_div.*sqrt(2*pi)).*exp(-(y(1)-gauss_mean2).^2./(2*gauss_div.^2)).*sin(m2.*y(2)-n2.*t);

rhs1=m_global.*A.*sin(m.*y(2));

% rhs2=((iota_a.*y(1).^2+iota_b)+2*(iota_a.*y(1).^2));
rhs2=iota_a.*y(1)+iota_b-1;
%rhs2=iota_a.*y(1)+iota_b-A1.*(y(1)-gauss_mean)./(gauss_div.^3.*sqrt(2*pi)).*exp(-(y(1)-gauss_mean).^2./(2*gauss_div.^2)).*cos(m1.*y(2)-n1.*t);
%rhs=((iota_a.*y(1).^2+iota_b)+2*(iota_a.*y(1).^2))-2*iota_a.*y(1).*t/(iota_a.*y(1).^2+iota_b).^2.*sin(y(2)-t./(iota_a.*y(1).^2+iota_b));

dydt = [rhs1; rhs2];
end

%% Function pendulum

function y = pendulum(t, y0)

global prec_global;
global m_global;
global iota_b_global;
global iota_a_global;
global A_global;
global iota_res_global;

g=-A_global;

phi_0=(iota_res_global-iota_b_global)./iota_a_global;
% 
% E=(iota_a_global.*(y0(1)-phi_0).^2)/2+A_global.*cos(m_global.*y0(2))+phi_0.*(iota_a_global.*phi_0/2+iota_b_global);
% 
% E_tilde=E+g-phi_0.*(iota_a_global.*phi_0/2+iota_b_global);

gamma_1_tot=m_global*y0(2);
m_period=floor((gamma_1_tot+pi)/(2*pi));
gamma_1=gamma_1_tot-2*pi*m_period;

% if gamma_1==-pi
%    gamma_1=gamma_1+2*pi;
%    m_period=m_period-1;
% end

if gamma_1>=pi || gamma_1<=-pi
   gamma_1 
   error('gamma_1=>pi or gamma_1<-pi')
end

p_1=sqrt(iota_a_global).*(y0(1)-phi_0);

E_tilde=p_1^2/2+g*(1-cos(gamma_1));

k=sqrt(E_tilde./(2.*g));

if k<1
    uT=4*ellipke(k,eps*prec_global);
else
    uT=2/k*ellipke(1/k,eps*prec_global);
end

T=sqrt(1/g)*uT;

if k<1
    xi_1=asin(sin(gamma_1/2)/k);
    t_phase=ellipticF(xi_1,k)/sqrt(g);
else
    xi_1=gamma_1/2;
    t_phase=ellipticF(xi_1,1/k)/(sqrt(g)*k);    
end


if p_1>0
    flag_p_p=1;
    p_p_init=1;
else
    flag_p_p=0;
    p_p_init=-1;
end

u=(p_p_init*(t-t(1))+t_phase)*sqrt(g);
% flag_p_p=0;


counter_even=0;
counter_period=0;
flag_neibour=1;
flag_neibour_2=1;

% y0(1)
% y0(2)
if isreal(xi_1)
    if k<1
        for i=1:length(t)
            [SN,CN,DN] = ellipj(u(i),k,eps*prec_global);
            gamma(i)=2*asin(k*SN);
            if p_p_init==1
                t_period=floor(u(i)/uT);
                if mod(u(i),uT)<=uT/4 || mod(u(i),uT)>=3*uT/4
                    p(i)=sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
                else
                    p(i)=-sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
                end
            elseif p_p_init==-1
                t_period=floor(u(i)/uT);
                if mod(u(i),uT)<=uT/4 || mod(u(i),uT)>=3*uT/4
                    p(i)=-sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
                else
                    p(i)=sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
                end
            end
            
            if isreal(p(i))
                yout(i,1)=p(i)/sqrt(iota_a_global)+phi_0;
            else
                yout(i,1)=nan;
            end
            yout(i,2)=(gamma(i)+2*pi*m_period)/m_global;
        end
    else
        for i=1:length(t)
            [SN,CN,DN] = ellipj(k*u(i),1/k,eps*prec_global);
            gamma_period(i)=2*asin(SN);
            t_period=floor((u(i)+uT/2)/uT);
            if mod(t_period,2)==0
                gamma(i)=2*pi*(t_period)+gamma_period(i);
            else
                gamma(i)=2*pi*(t_period)-gamma_period(i);
            end
            
            if p_p_init==1
                p(i)=sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
            elseif p_p_init==-1
                p(i)=-sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
            end
            
            if isreal(p(i))
                yout(i,1)=p(i)/sqrt(iota_a_global)+phi_0;
            else
                yout(i,1)=nan;
            end
            yout(i,2)=(gamma(i)+2*pi*m_period)/m_global;
        end
    end
else 
    yout(1:length(t),1)=nan;
    yout(1:length(t),2)=nan;
end


% if k<1
% %     for i=1:length(t)
% %         [SN,CN,DN] = ellipj(u(i),k);
% %         gamma(i)=2*asin(k*SN);
% %         if flag_p_p
% %             p(i)=sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
% %             if i>=3
% %                 if  (abs(p(i-1)-p(i-2))>=abs(p(i-1))) && flag_neibour_2
% %                     p(i)=-sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
% %                     flag_p_p=0;
% %                     flag_neibour_2=0;
% %                 else
% %                     flag_neibour_2=1;
% %                 end
% %             end
% %         else
% %             p(i)=-sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
% %             if i>=3
% %                 if  (abs(p(i-1)-p(i-2))>=abs(p(i-1))) && flag_neibour_2
% %                     p(i)=sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
% %                     flag_p_p=1;
% %                     flag_neibour_2=0;
% %                 else
% %                     flag_neibour_2=1;
% %                 end
% %             end
% %         end
% %         
% %         yout(i,1)=p(i)/sqrt(iota_a_global)+phi_0;
% %         yout(i,2)=(gamma(i)+2*pi*m_period)/m_global;
% %     end
% else
% %     for i=1:length(t)
% %         [SN,CN,DN] = ellipj(k*u(i),1/k);
% %         temp(i)=SN;
% %         gamma_period(i)=2*asin(SN);
% %         if mod(counter_period,2)==0
% %             gamma(i)=p_p_init*2*pi*(counter_period)+gamma_period(i);
% %         else
% %             gamma(i)=p_p_init*2*pi*(counter_period)-gamma_period(i);
% %         end
% %         
% %         if i>=3
% %             if (abs(gamma_period(i-1)-gamma_period(i-2))>abs(gamma_period(i-1)-pi)) && flag_neibour
% %                 counter_period=counter_period+1;
% %                 flag_neibour=0;
% %             elseif (abs(gamma_period(i-1)-gamma_period(i-2))>abs(gamma_period(i-1)+pi)) && flag_neibour
% %                 counter_period=counter_period+1;
% %                 flag_neibour=0;
% %             else
% %                 flag_neibour=1;
% %             end
% %         end
% %         
% %         if flag_p_p
% %             p(i)=sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
% %             if i>=3
% %                 if sign(-(p(i)-p(i-1))*(p(i-1)-p(i-2)))>0 && (abs(p(i-1)-p(i-2))>abs(p(i-1)))
% %                     p(i)=-sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
% %                     flag_p_p=0;
% %                 end
% %             end
% %         else
% %             p(i)=-sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
% %             if i>=3
% %                 if (sign(-(p(i)-p(i-1))*(p(i-1)-p(i-2)))>0) && (abs(p(i-1)-p(i-2))>abs(p(i-1)))
% %                     p(i)=sqrt(2*(E_tilde-g*(1-cos(gamma(i)))));
% %                     flag_p_p=1;
% %                 end
% %             end
% %         end
% %         
% %         yout(i,1)=p(i)/sqrt(iota_a_global)+phi_0;
% %         yout(i,2)=(gamma(i)+2*pi*m_period)/m_global;
% %     end
% end
y=yout;
end