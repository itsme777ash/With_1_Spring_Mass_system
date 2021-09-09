init_freq = 1;
Final_freq = 350;

% plot_Steel = zeros(length(init_freq:1:Final_freq),1);
% plot_PVC = zeros(length(init_freq:1:Final_freq),1);
% plot_Steel_PVC = zeros(length(init_freq:1:Final_freq),1);
plot_Steel = BEAM_MAIN_program(210e9,7850,210e9,7850);
plot_Steel_PVC = BEAM_MAIN_program(3.4e9,1380,210e9,7850);
plot_PVC = BEAM_MAIN_program(3.4e9,1380,3.4e9,1380);
semilogy(linspace(init_freq,Final_freq,Final_freq-init_freq+1),plot_Steel,'blue',linspace(init_freq,Final_freq,Final_freq-init_freq+1),plot_PVC,'green',linspace(init_freq,Final_freq,Final_freq-init_freq+1),plot_Steel_PVC,'magenta');
xlabel('Frequency Applied(Hz)'); ylabel('Normalized Amplitude(Xk/F)');
