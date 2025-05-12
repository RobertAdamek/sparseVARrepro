FROM rocker/r-ver:4.4.2

RUN mkdir /home/output

ENV SCRIPT_NAME=all.R

RUN R -e "options(repos= list(CRAN = 'https://packagemanager.posit.co/cran/2025-01-01')); install.packages(c('Rcpp', 'RcppParallel', 'parallelly', 'RcppArmadillo','RcppThread'))" 

# package related files
COPY DGP.R /home/DGP.R
COPY Simulation.R /home/Simulation.R
COPY sparseVARboot_0.4.4.tar.gz /home/sparseVARboot_0.4.4.tar.gz

RUN R -e "install.packages('/home/sparseVARboot_0.4.4.tar.gz', repos = NULL, type = 'source')" 

# scripts
COPY all.R /home/all.R
# DGP 2
COPY dgp2_size.R /home/dgp2_size.R
COPY dgp2_power_mu0035_prop05.R /home/dgp2_power_mu0035_prop05.R
COPY dgp2_power_mu00175_prop05.R /home/dgp2_power_mu00175_prop05.R
COPY dgp2_power_mu00175_prop09.R /home/dgp2_power_mu00175_prop09.R
# DGP 1
COPY dgp1_size.R /home/dgp1_size.R
COPY dgp1_power_mu05_prop05.R /home/dgp1_power_mu05_prop05.R
# DGP 9 
COPY dgp9_size.R /home/dgp9_size.R
COPY dgp9_power_mu0013_prop05.R /home/dgp9_power_mu0013_prop05.R
# DGP 0
COPY dgp0_size.R /home/dgp0_size.R
COPY dgp0_power_mu025_prop05.R /home/dgp0_power_mu025_prop05.R

ENTRYPOINT ["sh", "-c", "Rscript /home/${SCRIPT_NAME}"]