# kcoverage_single_R

REQUIREMENTS: The R packages R2Cuba and randtoolbox for functions "cuhre" and "sobol", which are called by funJn.R

If you use this code for published work, please cite paper[1] listed below.

The scripts calculate the k-coverage probability (based on SINR* values) in a single-tier cellular network using a method based on a homogeneous Poisson process model. More details are found in the (submitted) work [1], which presents the model that these scripts are based on.

The script funProbCov.R uses a inclusion-exclusion-like formula and two types of integral to calculate the k-coverage probability in a network with log-normal shadowing (though the shadowing distribution can be somewhat arbitrary [1]) and without fading.

The simpler integral I_n uses a quadrature method or a simple analytic formula (for the zero-noise or “interference limited” case). The more complex high-dimensional integral J_n uses quadrature methods for low dimensions and quasi-random (Sobol) integration for higher (n>2) dimensions.

Running the file TestSimVsInt.R is a good place to start. It will compare the expressions to simulation results (based on the same model from [1]). Default parameters are based on a Walfisch-Ikegami model for a urban environment. The published work [1] has more details on the equations and the model.

*SINR = signal-to-interference-and-noise-ratio

[1] H.P Keeler, B. Błaszczyszyn and M. Karray, 'SINR-based k-coverage probability in cellular networks with arbitrary shadowing', 2013. Online at http://arxiv.org/abs/1301.6491
