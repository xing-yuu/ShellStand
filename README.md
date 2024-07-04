# Shell stand: Stable thin shell models for 3D fabrication

This repository offers an implementation of *Shell stand* and a library for its algorithm. It is about Open surface model stability optimization; see our CVM 2024 paper for details, DOI: https://doi.org/10.1007/s41095-024-0402-8.

### Pipeline
![image](pipeline.jpg)
Shell optimization overview. Given an input shell model and its initial orientation, we perform global deformation and local thickening modulations to produce a balanced shell model that can be fabricated and stand stably on the ground.
### Dependencies
- CinoLib (ARAP) https://github.com/mlivesu/cinolib
- Sequential Line Search (Bayesian optimization) https://github.com/yuki-koyama/sequential-line-search.git
- MOSEK https://www.mosek.com/
- VTK 8.2.0
- CGAL
- Eigen3
- Qt5
- Yaml-cpp 0.3.0 
