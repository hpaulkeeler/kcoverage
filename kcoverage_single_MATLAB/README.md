# kcoverage_single_MATLAB

If you use this code for published work, please cite paper[1] listed below.

These scripts calculate the k-coverage probability (based on SINR* values) in a single-tier cellular network using a method based on a homogeneous Poisson process model. More details are found in the (submitted) work [1], which presents the model that these scripts are based on.

The script funProbCov.m uses a inclusion-exclusion-like formula and two types of integral to calculate the k-coverage probability in a network with log-normal shadowing (though the shadowing distribution can be somewhat arbitrary [1]) and without fading.

The simpler integral I_n uses a quadrature method or a simple analytic formula (for the zero-noise or “interference limited” case). The more complex high-dimensional integral J_n uses quadrature methods for low dimensions and quasi-random (Sobol) integration for higher (n>2) dimensions.

The script funProbCovFade.m calculates 1-coverage probability for a network with Rayleigh fading (exponentially distributed with unit mean) and log-normal shadowing. Numerical integration of hypergeometric function 2F1 is used when the model has noise. A close-form solution with 2F1 is used in the no noise case.

Simulation scripts are also included for comparison purposes. All network base stations are sampled on a disk region. The disk region needs to be large enough to reduce "edge effects", which become more prominent when fading is included.

Running the file TestSimVsInt.m is a good place to start. The fading-incorporated model is demonstrated in TestSimVsIntFade.m. Default parameters are based on a Walfisch-Ikegami model for a urban environment. The (submitted) associated work [1] has more details on the equations and the model.

*SINR = signal-to-interference-and-noise-ratio

[1] H.P Keeler, B. Błaszczyszyn and M. Karray, 'SINR-based k-coverage probability in cellular networks with arbitrary shadowing', 2013. Online at http://arxiv.org/abs/1301.6491 or https://ieeexplore.ieee.org/document/6620410/.

