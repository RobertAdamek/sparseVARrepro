# install package
setwd(this.path::here())
install.packages("../sparseVARboot_0.5.0.tar.gz", repos = NULL, type = "source")
library(sparseVARboot)

# But for us, if we want to change things, the following way is easier
# Note, this does not work for parallel computing. Then we really need to install the package. 
# setwd(this.path::here())
# setwd("../sparseVARboot")
# devtools::load_all()
# setwd(this.path::here())

# DGP 2
source("dgp2_size.R")
source("dgp2_power_mu00175_prop05.R")
source("dgp2_power_mu0035_prop05.R")
source("dgp2_power_mu00175_prop09.R")

# DGP 1
source("dgp1_size.R")
source("dgp1_power_mu05_prop05.R")

# DGP 9
source("dgp9_size.R")
source("dgp9_power_mu0013_prop05.R")

# DGP 0
source("dgp0_size.R")
source("dgp0_power_mu025_prop05.R")
