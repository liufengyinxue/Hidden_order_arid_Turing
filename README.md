# Hidden_order_arid_Turing

"Raster_DF_Convolution.m" : calculate density fluctuation by the convolution method.

"Raster_DF_ThrowBox.m" : calculate density fluctuation by the box-throwing method.

"Arid_***_GPU.m" : simulate the water-vegetation models in MATLAB 2022a, ***={Rietkerk, Klausmeier, Getzin, Hardenberg}

"Arid_***_GPU.ipynb" and "SpatialFunction_iPy.cl": simulate the water-vegetation models in pyopencl, ***={Klausmeier, Getzin, Hardenberg}

“Rietkerk_SPDE_GPU.m” and "get_twod_dW_GPU.m": simulatie the stochastic (white noise) Rietkerk model.

"spectral.m" and "Spectral_density.m": do spectral analysis for matrix data. The only difference between them is whether or not the fast Fourier transform is used. The results from them are qualitatively consistent.

"BG_UPL": calculate unchanneled path length of bare ground.

"LOC_dis.m" and "GL_dis.m": excute two kinds of disturbance.

“SpotStat.m”: count the number of spot vegetation patches.

"Niger_image.zip": the aerial photographs cite from Barbier et al., Journal of Ecology, 2006, 94, 537–547.

"Figure1_classfied.zip": the natural patterns in Figure1A1-A4.
