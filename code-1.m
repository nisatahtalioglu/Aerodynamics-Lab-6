clc; clear;

% Lab_6---------------------------------------------------------------%

c = 0.15;                  % chord length [m]
p_atm=150*9.80665;
p_static = 125*9.80665;
q_max = 200;                % [Pa]
V_inf = 18.3811;            % calculated in the report
rho = 1.1839;               % given in the report

 
N_l = readmatrix("Aerodynamics_Lab6.xlsx","Range","N7:W21");   % [mmH20] lower surface pressure 
N_u = readmatrix("Aerodynamics_Lab6.xlsx","Range","D7:M21");    % [mmH20] upper surface pressure 


x_l = readmatrix("Aerodynamics_Lab6.xlsx","Range","D35:D44")';     % [mm] lower surface x coord
x_u = readmatrix("Aerodynamics_Lab6.xlsx","Range","D25:D34")';      % [mm] upper surface x coord

x_c_u = x_u / c;
x_c_l = x_l / c;

y_l = readmatrix("Aerodynamics_Lab6.xlsx","Range","E35:E44");
y_u = readmatrix("Aerodynamics_Lab6.xlsx","Range","E25:E34");

bound = linspace(0,0.15,10);    % integral boundaries

for i = 1:15
    for j = 1:10
        C_p_l(i,j) = (N_l(i,j)*9.80665 + p_atm - p_static) / q_max;     % lower surface Cp 
        C_p_u(i,j) = (N_u(i,j)*9.80665 + p_atm - p_static) / q_max;     % upper surface Cp 
       
    end

    integral_2(i,1) = trapz(bound,C_p_l(i,:) - C_p_u(i,:));
   
    C_n(i) = 1/c * integral_2(i);       % normal force coeff 
    

    integral_4(i,1) = trapz(bound,C_p_u(i,:)*(transpose(y_u)/x_u) - C_p_l(i,:)*(transpose(y_l)/x_l));
    
    C_a(i) = 1/c * integral_4(i);       % axial force coeff 
    
    integral_6(i,1) = trapz(bound, C_p_u(i,:).*x_u - C_p_l(i,:).*x_l);
    integral_8(i,1) = trapz(bound, C_p_u(i,:).*(transpose(y_u)/x_u).*transpose(y_u));
    integral_10(i,1)= trapz(bound, -C_p_l(i,:).*(transpose(y_l)/x_l).*transpose(y_l));

    C_m_le(i) = 1/(c^2) * (integral_6(i,1) + integral_8(i,1) + integral_10(i,1));
end


 plot(x_c_l,C_p_l(15,:))
 hold on
 plot(x_c_u,C_p_u(15,:))
 xlabel("x/c")
 ylabel("C_p")
 title("C_p vs x/c for \alpha = 14 [deg]")
 legend("Lower surface","Upper surface")
 hold off


alpha = readmatrix("Aerodynamics_Lab6.xlsx","Range","C7:C21")';


C_d_theo = readmatrix("Aerodynamics_Lab6.xlsx","Range","H25:H39");    % theoretical cd from xfoil 
C_d = C_n.*sind(alpha) + C_a.*cosd(alpha);



%plot(alpha,C_d)
%hold on
%plot(alpha,C_d_theo')
%xlabel("\alpha (degrees)")
%ylabel("C_D")
%legend("Experimental","Theoretical")
%title("C_D vs \alpha")
%hold off


C_l = C_n.*cosd(alpha) - C_a.*sind(alpha);
C_l_theo = readmatrix("Aerodynamics_Lab6.xlsx", "Range", "I25:I39");

%plot(alpha,C_l)
%hold on 
%plot(alpha,C_l_theo)
%xlabel("\alpha (degrees)")
%ylabel("C_l")
%legend("Experimental","Theoretical")
%title("C_l vs \alpha")
%hold off


C_m_le_theo = readmatrix("Aerodynamics_Lab6.xlsx", "Range", "J25:J39");

%plot(alpha,C_m_le)
%hold on
%plot(alpha,C_m_le_theo)
%xlabel("\alpha (degrees)")
%ylabel("C_m_,LE")
%legend("Experimental", "Theoretical")
%title("C_m_,LE vs \alpha")
%hold off

% plot(C_d,C_l)
% hold on
% plot(C_d_theo,C_l_theo)
% xlabel("C_d")
% ylabel("C_l")
% legend("Experimental", "Theoretical")
% title("C_l vs C_d")
% hold off
