# load requirements to list reqlibs
reqlibs <- scan("../REQUIREMENTS.txt", what="", sep="\n")

# point to default CRAN mirror
cran_url <- "http://cran.us.r-project.org"
# install all required packages
for(lib in reqlibs) {
  install.packages(lib, repos = cran_url)
}

# test env
source("install/testenv.R")