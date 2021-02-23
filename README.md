# Work in progress
This is (going to be) a rather disorganised pile of utilities for drawing samples from bosonic system probability amplitude. In near future I hope to have first stable release ready, including :
* a simple paralel implementation of diffusion monte carlo sampling algorithm (based of [this article](https://arxiv.org/pdf/physics/9702023.pdf))
* an optimized paralel implemenntation of Zdravko Botev's [adaptive kernel density estimation algorithm](https://www.mathworks.com/matlabcentral/fileexchange/58312-kernel-density-estimator-for-high-dimensions) + some compiled data regarding its performance
* a simple paralel implementation of Metropolis sampling algorithm


Everything is written in Fortran, using some Lapack, Blas and MPI


(note for collaborators : ) různé další soubory (grafy atp.) [zde](https://drive.google.com/drive/folders/1GCk214wVyhbHRdO6Afs1KrC-Fhj9AN-_?usp=sharing).
