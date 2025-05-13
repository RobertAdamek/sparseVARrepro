# sparseVARrepro
This repository contains everything necessary to reproduce the simulation results of the paper "Sparse high-dimensional vector autoregressive bootstrap" by Robert Adamek, Stephan Smeekes, and Ines Wilms. 
If you only want to re-run the simulations yourself, we provide a Docker image and instructions how to run it. This method is reproducible regardless of future updates to R or its packages. For transparency, and to make the code useable for other purposes, we also provide the source code used to build the Docker image. 

# How to run the simulations via Docker
*Warning: the simulations can take a significant amount of time to run. We ran them on an Intel® Core™ i9-10900K processor for approximately 360 hours. It also requires at least 16GB of RAM.*
1. Download our Docker image from https://drive.google.com/file/d/12LdcaW7fsZSH39hU8c-CCfhjWSOZ-3Ge/view?usp=sharing. This file is called "sparsevar-final.tar", and is quite large (~976MB). 
2. Download and install Docker Desktop from https://www.docker.com/.
    - Docker Desktop has a user-friendly interface which can also be used to do the following steps, but here we provide precise commands you can use in your computer's Terminal/PowerShell.
3. Start your Terminal/Powershell and navigate to the folder where you saved "sparsevar-final.tar".
    - On Windows, you can open PowerShell by opening the start menu (pressing the Windows key) and typing `powershell`. On Mac, use Command+Space to search for `terminal`.
    - Your starting working directory (on Windows) is probably "C:\Users\yourname". You can tell where you are by looking left of the ">" sign next to your flashing input line (it's "%" in the Mac Terminal).
    - If you're ever lost, the command `pwd` (print working directory) will also tell you where you are.
    - To navigate to a specific folder, use `cd path\to\folder` to set that as your working directory. For example, if the file is in your default Windows downloads folder, you can get there directly with `cd C:\Users\yourname\Downloads`.
    - If you want to quickly return to your home folder, type `cd ~`. "~" is shorthand for "C:\Users\yourname".
    - Insted of giving the full address of the folder (the absolute path) you can also give the relative path: If you're currently in "C:\Users\yourname", you can just type `cd Downloads`, because "Downloads" is a folder in your working directory.
    - To see the contents of the current folder you're in, "ls" will list all the files and folders, along with some metadata. If you're in the right place, "sparsevar-final.tar" should be on the list.
4. Load the Docker image with `docker load -i .\sparsevar-final.tar`. This could take a few minutes.
    - The "." is shorthand for your current working directory, the same way that "~" is shorthand for your home directory. 
5. Check that the image loaded correctly: `docker images` should list all the images you have access to. "sparsevar-final" should be one of them.
6. Create a container that runs the image: `docker run -v .:/home/output sparsevar-final`
    - The default command to create a container is `docker run image-name`; but the optional argument "-v .:/home/output" will copy the simulation results into your current working directory.
    - This method runs the entire simulation. If you want to run the simulation in smaller chunks, you can also add another optional arugment, e.g.: `docker run -e SCRIPT_NAME=dgp0_size.R -v .:/home/output sparsevar-final`. This runs only the "dgp0_size.R" script. There are 10 scripts in total, see the "Docker image setup" folder. Each takes approximately 36 hours to run.
7. When complete, your folder should contain ten .RData files, which contain the simulation results. To generate the plots, run the "plotting results.R" script in the "Results" folder.

# Contents
- The simulation code uses function from the bespoke package "sparseVARboot". We provide the bundled package in the file "sparseVARboot_0.4.4.tar.gz", which can be installed directly, and also its source files in the sparseVARboot folder.

- The folder "Docker image setup" contains all the files used to build the Docker image, including 10 individual scripts that can also be run directly (with minor alterations).

- The "Results" folder contains the files the simulation produces, as well as a script "plotting results.R" which produces plots from these files. The plots themselves are also included in pdf format.


