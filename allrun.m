% Run all simulations
function allrun
monte_sim_he_all(10,10000);
monte_sim_fft_all(10,1000);
monte_sim_uc_all(10,10000);
monte_sim_condat_all(10,1000);
end