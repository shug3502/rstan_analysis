%call output from mrna_transport_sims_and_plot
% and plot this on grid of NCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distn = csvread('mrna_sims.csv');
ts = NaN*ones(1,15); %may want to include times on plots, can work out way to do so
plot_coarse_grained_NCs(distn(1,:),NaN,0);
print('plots/simsv005','-deps');
plot_coarse_grained_NCs(distn(15,:),NaN,0);
print('plots/simsv006','-deps');
exit;
