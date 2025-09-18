# Li-Yang WRF
Based on the Rossby wave ray tracing theory in horizontally nonuniform basic flows, Yang and Li (2025) introduced the concept of the horizontal Rossby wave ray flux (WRF), hereafter referred to as the Li–Yang WRF, which is developed to diagnose the propagation pathways, local activity, and propagation directions of Rossby waves. For a comprehensive account of the theoretical framework and methodological details, the reader is referred to the original publication.

The computational package comprises three functions (Fun1_threshold, Fun2_region_threshold, Fun3_WRF_calculate) and a main script (WRF_universal), all provided as .py files, with the following functionalities:

Fun1_threshold defines a velocity threshold to exclude anomalous wave numbers (optional);

Fun2_region_threshold identifies wave rays that pass through a specified target region;

Fun3_WRF_calculate evaluates the wave ray flux (Li–Yang WRF).

All procedures are implemented within the main script (WRF_universal).


References
Yang, Y. N., and J. P. Li*, 2025: Novel monsoon indices based on vector projection and directed angle for measuring the East Asian summer monsoon. Clim. Dyn., 63, 210. https://doi.org/10.1007/s00382-025-07696-7
