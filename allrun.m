% Run all simulations
function allrun(M,N)
monte_sim_he_all(M,N);
monte_sim_fft_all(M,N);
monte_sim_uc_all(M,N);
monte_sim_condat_all(M,N);
end