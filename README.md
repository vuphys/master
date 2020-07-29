# Master Project
Master Project in Medical Physics at the University of Aberdeen

"Applying and optimising Total Generalised Variation filtering for FFC-MRI images"

Hoang Vu Nguyen

Supervisor: Dr. Lionel Broche

For lastest simulation, please check "monte_sim_x_all" with "x" is the name of the algorithms (e.g., condat)

The "auto TGV" program, which can be used for denoising FFC MRI images automatically, can be found in auto_TGV folder.

All data are stored in DATA

I would be glad to receive any comment or suggestion.

For further information and details of this project. Please send me an email at: harrynguyen.hv@gmail.com

ABSTRACT:

Fast Field Cycling Magnetic Resonance Imaging (FFC MRI) is a new technique that promises to gain more valuable information than conventional MRI. FFC MRI has to deal with many types of noise in its image, which degrades the image quality and reduces clinical information in practice. A filtering method called Total Generalised Variational (TGV) offers high performance in denoising, and preserving the edges of the structure in the image can be a solution for the FFC MRI problem. The project considered four open-source TGV algorithms and conducted several tests on FFC MRI data using a Ground Truth (GT) image to determine and validate TGV behaviour on FFC MRI images. Tests were designed to use the Monte Carlo method in the simulation of the TGV process. Six Full Reference Image Quality Assessment (FR IQA) metrics were also used to quantify the image quality. The research successfully determined the relationship between TGV parameters and the amplitude of the noise, which led to constructing a method to choose TGV parameters for the best outcome. A suitable metric with the highest correlation to the noise was determined, and best-performing algorithms were also revealed. Finally, an auto TGV denoising procedure was developed to produce a process of denoising FFC MRI images automatically. An adapted computer program of the procedure was also built for practical use and later research.
