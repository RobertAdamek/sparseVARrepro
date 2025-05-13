# sparseVARrepro
This repository contains everything necessary to reproduce the simulation results of the paper "Sparse high-dimensional vector autoregressive bootstrap" by Robert Adamek, Stephan Smeekes, and Ines Wilms. 
If you only want to re-run the simulations yourself, we provide a Docker image and instructions how to run it. This method is reproducible regardless of future updates to R or its packages. For transparency, and to make the code useable for other purposes, we also provide the source code used to build the Docker image. 

# How to run the simulations via Docker
1. Download the Docker image from https://drive.google.com/file/d/12LdcaW7fsZSH39hU8c-CCfhjWSOZ-3Ge/view?usp=sharing. This file is called "sparsevar-final.tar", and is quite large (~976MB) 
The code we provide here is only for transparency, and you do not need to run 
https://www.docker.com/

# Contents
- The simulation code uses function from the bespoke package "sparseVARboot". We provide the bundled package in the file "sparseVARboot_0.4.4.tar.gz", which can be installed directly, and also its source files in the "sparseVARboot" folder.

- The folder "Docker image setup" contains all the files used to build the Docker image, including 10 individual scripts that can also be run directly (with minor alterations).

- The "Results" folder contains the files the simulation produces, as well as a script "plotting results.R" which produces plots from these files. The plots themselves are also included in pdf format.


