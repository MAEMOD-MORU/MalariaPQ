source('function.R')

Tanzania_k13_Allele_frequency <- read.csv("data/Tanzania_K13_Allele_frequency.csv")[,2]

out <- run_malariaPQ_ode(prob_sym_s = 0,prob_sym_r=0)

plot(out$time, out$inc_sym, type='l', xlab='Time (months)', ylab='Incidence of symptomatic infections', main='Incidence over time')

#plot inc_sym and inc_asym
plot(out$time, out$inc_asym, type='l', xlab='Time (months)', ylab='Incidence of asymptomatic infections', main='Incidence over time')
lines(out$time, out$inc_sym, col='red')

out  <- run_malariaPQ_ode(prob_sym_s = 0.25, prob_sym_r = 0.25)

prep <- prep_malaria_outputs(out, start_year = 2000, end_year = 2050)

plot_malaria_outputs(prep,
                     plot_type = "Sym and Asym",
                     year_plot = 2000,
                     Tanzania_Incidence = Tanzania_Incidence,
                     Tanzania_k13_Allele_frequency = Tanzania_k13_Allele_frequency)
