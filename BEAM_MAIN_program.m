 function amp_list = BEAM_MAIN_program(Young_mat_1,Den_mat_1,Young_mat_2,Den_mat_2)
%Young_mat_1 = 210e9;  Den_mat_1 = 7850;
%Young_mat_2 = 210e9;  Den_mat_2 = 7850;

E = [Young_mat_2,Young_mat_1]; q = 0; P = 100; M = 0; density = [Den_mat_2,Den_mat_1];b = [0.05,0.05]; h = [0.01,0.01]; I =b.*(h.^3)/12;
damp_ratio = 0.1;
Final_freq = 350;
init_freq = 1;
Point_Study = 21;

%le = 1;
Element_no = 20;
elem_len = [0.08 , 0.02] ;  % (Element 2, Element 1)
Nodes = Element_no + 1;
DOF = Nodes*2;
f = zeros(4,1,Element_no);
k = zeros(4,4,Element_no);
Mass = zeros(4,4,Element_no);
%Element stiffness matrix and load vector
for element = 1:1:Element_no
    elem_s_no = mod(element,2) + 1;
    k(:,:,element) = Beam_stiff_matrix(E(elem_s_no),I(elem_s_no),[(element-1)*elem_len(elem_s_no),element*elem_len(elem_s_no)]);
    Mass(:,:,element) = Beam_Mass_matrix(b(elem_s_no)*h(elem_s_no),density(elem_s_no),elem_len(elem_s_no));
    f(:,:,element) = Load_vec_matrix(q,[(element-1)*elem_len(elem_s_no),element*elem_len(elem_s_no)]);
end

K_final = zeros(DOF,DOF);
M_final = zeros(DOF,DOF);
K_final([1:4],[1:4]) = k(:,:,1);
M_final([1:4],[1:4]) = Mass(:,:,1);
F_final = zeros(DOF,1);
F_final(1:4) = f(:,:,1);
%Assembling
for element = 1:1:(Element_no-1)
    K_final([(element*2)+1:(element*2)+4],[(element*2)+1:(element*2)+4]) = K_final([(element*2)+1:(element*2)+4],[(element*2)+1:(element*2)+4]) + k(:,:,element+1);
    M_final([(element*2)+1:(element*2)+4],[(element*2)+1:(element*2)+4]) = M_final([(element*2)+1:(element*2)+4],[(element*2)+1:(element*2)+4]) + Mass(:,:,element+1);
    F_final([(element*2)+1:(element*2)+4]) = F_final([(element*2)+1:(element*2)+4]) + f(:,:,element+1);
end

% k1 = Beam_stiff_matrix(E,I,[0,le]);
% f = Load_vec_matrix(q,[0,le]);
F_final(DOF-1) = F_final(DOF-1) + P;
F_final(DOF) = F_final(DOF) + M;
% f(3) = f(3) + P;
% f(4) = f(4) + M;

%Adding Spring mass system 

k_spring = 20000;
m_spring = 20;
spring_node = 13;
DOF = DOF + 1;

%Altering stiffness matrix
K_final(DOF,DOF) = 0;
M_final(DOF,DOF) = 0;
F_final(DOF,1) = 0;
K_final(DOF,spring_node) = K_final(DOF,spring_node) - k_spring;
K_final(DOF,DOF) = K_final(DOF,DOF) + k_spring;
K_final(spring_node,spring_node) = K_final(spring_node,spring_node) + k_spring;
K_final(spring_node,DOF) = K_final(spring_node,DOF) - k_spring;
%Altering Mass matrix
M_final(DOF,DOF) = M_final(DOF,DOF) + m_spring;

%Done adding spring mass system


%Solving equations for deflection ( Essential Boundary conditions : w1 = 0
%and (dw/dx)1 = 0)
W = zeros(DOF,1);

Kreduced = K_final([3:DOF],[3:DOF]);
Freduced = F_final([3:DOF]);
Mreduced = M_final([3:DOF],[3:DOF]);
K_inv = inv(Kreduced);
wreduced = K_inv*Freduced;

W(3:DOF) = wreduced;
Reaction = K_final*W;
disp("SUCCESS");

%Natural Frequency
F = zeros(DOF-2,1);
F(39,1) = 100;
freq_sqr = eig(Kreduced,Mreduced);
freq_root = sqrt(freq_sqr);
[V,D] = eig(Kreduced,Mreduced);
%disp(transpose(V(:,1))*Mreduced*V(:,1));

%Transforming from Physical to Modal Coordinates y(t)
Mod_F = transpose(V) * F;
Mod_M = transpose(V) * Mreduced * V;
Mod_K = transpose(V) * Kreduced * V;
Solut = cell(DOF-2,1);
Y_real = cell(DOF-2,1);
Y_imag = cell(DOF-2,1);
Y_mag = cell(DOF-2,1);

%Solving Equations in Algebric form
for i = 1:1:DOF-2
    syms omega_applied
    m = freq_sqr(i,1) - (omega_applied)^2;
    n = 2*damp_ratio*freq_root(i,1)*omega_applied;
    Y_real{i} = (Mod_F(i,1)*m)/((m*m)+(n*n));
    Y_imag{i} = (Mod_F(i,1)*n)/((m*m)+(n*n));
    Y_mag{i} = sqrt((Y_real{i}^2) + (Y_imag{i}^2));
end

%Modal superposition of solution from algebric method
Y_val = zeros(DOF-2,1);
Amp_Y_real = zeros(length(init_freq:1:Final_freq),1);
for f = init_freq:1:Final_freq
    omeg = 2*pi*f;
    for i = 1:1:DOF-2
        Y_val(i,1) = subs(Y_mag{i},omega_applied,omeg);
    end
    temp_amp_per_mode = 0;
    for j = 1:1:DOF-2
        temp_amp_per_mode = temp_amp_per_mode + ((V(Point_Study-2,j)*Y_val(j))^2);
    end
   Amp_Y_real(f-init_freq+1,1) = sqrt(temp_amp_per_mode)/W(Point_Study);
end

amp_list = Amp_Y_real;

%Solving Second Order equations to find y(t)

%Method - 1

% for i = 1:1:40
%     C = 2*damp_ratio*freq_root(i);
%     syms y(t) omega
%     eqn = diff(y,t,2) == -Mod_K(i,i)*y + Mod_F(i,1)*cos(omega*t) - C*diff(y,t);
%     Dy = diff(y,t);
%     cond = [y(0)==0, Dy(0)==0];
%     ySol = dsolve(eqn,cond);
%     ySo_simp = simplify(ySol);
%     Solut{i} = ySo_simp;
% end
% 
% %disp(subs(Solut{1},omega,52));
% %disp(Func_Max(Solut{1}));
% % temp = zeros(length(1:1:20),1);
% % 
% % for tim = 1:1:4000 
% %     temp(tim,1) = subs(Solut{1},t,tim/1000);
% % end
% % 
% % plot(linspace(0.001,4,4000),temp(:,1));
% 
% %RMS Amplitude
% 
% Amp_rms = zeros(length(5:1:100),1);
% Normalized_amp = zeros(length(5:1:100),1);
% %Uncomment quick
% 
% for f = 5:1:100
%     omeg = 2*pi*f;
%     inter_soln = cell(40,1);
%     for j = 1:1:40
%         inter_soln{j} = subs(Solut{j},omega,omeg);
%     end
% %     inter_soln = subs(Solut,omega,f);
%     temp_amp_sqr = 0;
%     for i = 1:1:40
%         inter_ans = V(39,i)*inter_soln{i};
%         a = subs(inter_ans,t,pi/(2*omeg))/cos((freq_root(i)/omeg)*pi/2);
%         b = subs(inter_ans,t,pi/(2*freq_root(i)))/cos((omeg/freq_root(i))*pi/2);
%         temp_amp_sqr = temp_amp_sqr + (a^2) + (b^2);
%     end
%     Amp_rms(f-4,1) = sqrt(temp_amp_sqr)/1.414;
%     Normalized_amp(f-4,1) = Amp_rms(f-4,1)/(W(DOF-1));
% end

%Plotting

% plot(linspace(init_freq,Final_freq,Final_freq-init_freq+1),Amp_Y_real);
% xlabel('Frequency Applied(Hz)'); ylabel('Normalized Amplitude(Xk/F)');

%Plotting ends

% amp_var = zeros(length(1:1:100),1);
% for tim = 1:1:100
%     amp_var(i) = subs(Solut{1},[t,omega],[tim/100,2*pi*8]);
%     
% end
% plot(linspace(1,100,100),amp_var);

%Modal Superposition

% Final_ans = zeros(40,1);
% for i = 1:1:40
%     Final_ans = Final_ans + V(:,i)*Solut{i};
% end


%Amplitude vs Omega Calculation and Plot

% temp = zeros(length(1:1:300),1);
% y = Final_ans(39);
% % Amp = zeros(length(5:1:100),1);
% % for f = 5:1:100
% y_2 = subs(y,omega,2*pi*8);
%     for tim = 1:1:300 
%         temp(tim,1) = subs(y_2,t,tim/50);
%     end
% %Amp(f-4,1) = max(temp(:,1));
% %end
% 
% % plot(linspace(5,100,96),Amp);
% % xlabel('Frequency Applied(Hz)'); ylabel('Amplitude(mm)');
% plot(linspace(1,1,300),temp(:,1));

% %Free Vibration Modes
% d = zeros(40,length([1,0.1,10]));
% for t = 1:1:20
%     d(:,t) = V(:,1)*cos(sqrt(D(1,1))*t/10);     %Mode no
% end
% W_freq_1 = zeros(1,21);
% W_freq_1(2:21) = d(1:2:39,1);
% %plot(linspace(0,1,21),W_freq_1,linspace(0,1,21),zeros(length(linspace(0,1,21))),'r');
% %plot(linspace(1,40,20),d(1:39:20,:));
% freq = (1/(2*pi))*sqrt(freq_sqr);
%Mass_mat = Beam_Mass_matrix(1,density,le);


%Plotting the solution

%  w_final = zeros(Element_no*length([-1,0.1,1]),1);
%  w_actual = zeros(Element_no*length([-1,0.1,1]),1);
%  X_vals = zeros(Element_no*length([-1,0.1,1]),1);
%  i = 1;
%  for element = 1:1:Element_no
%      elem_s_no = mod(element,2) + 1;
%  for e = -1 : 0.1 : 1
%      N1 = 0.25*(1-e)*(1-e)*(2+e);
%      N2 = 0.25*(1-e)*(1-e)*(1+e);
%      N3 = 0.25*(1+e)*(1+e)*(2-e);
%      N4 = 0.25*(1+e)*(1+e)*(e-1);
%      w_final(i) = N1*W(2*element - 1) + elem_len*0.5*N2*W(2*element) + N3*W(2*element + 1) + elem_len*0.5*N4*W(2*element +2);
%      X_fin = (1+e)*0.5*element*elem_len + (1-e)*0.5*(element-1)*elem_len;
%      X_vals(i) = X_fin;
%      w_actual(i) = (P*X_fin*X_fin/(6*E(elem_s_no)*I(elem_s_no)))*(3*le - X_fin);
%      i = i+1;
%  end
%  end

% plot(X_vals,-1*w_final,'*',X_vals,-1*w_actual);
% legend('FEM results','Analytical Solution');
% title('20 Element Solution');
% xlabel('X(m)');ylabel('Deflection(m)');
% plot(linspace(1,DOF-2,40),freq(1:DOF-2),'*',linspace(1,DOF-2,40),freq(1:DOF-2));
% xlabel('DOF');ylabel('Natural Frequency(Hz)');
% x = linspace(0,10,21);
% plot(x,-1*w_final,'blue');
% plot(x,-1*w_actual,'red');
% title('Analytical Solution'); xlabel('x(m)');ylabel('Deflection(m)');



