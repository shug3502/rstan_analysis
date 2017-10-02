%call output from mrna_transport_sims_and_plot
% and plot this on grid of NCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distn = csvread('mrna_sims.csv');
ts = NaN; %may want to include times on plots, can work out way to do so
plot_coarse_grained_NCs(distn,ts,1);
print('plots/simsv005','-deps');

