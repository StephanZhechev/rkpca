# rkpca
Randomized kernel PCA

We implement our own version of randomized kernel PCA that we use for non-linear dimensionality reduction. 

# Kernel PCA

Kernel PCA is generalization of principal component analysis (PCA), where we first embed our data into a high dimensional space (potentially infinite dimensional), where we perform regular PCA. The embedding is in general non-linear and this procedure lets us capture complicated non-linear relations between our initial features. For more information on PCA and kernel PCA, we refer to [wikipedia](https://en.wikipedia.org/wiki/Kernel_principal_component_analysis).

The problem with KPCA is that it has complexity O(n^3), which scales horribly with the size of the data set. However, there are several approximations of KPCA that are substantially faster. Here we implement two methods, namely random Fourier features and the Nyström low-rank approximation. The main difference between the two is that the former generates random features in a data-independent way, while the latter makes use of random sampling from the data. Unfortunately, in R there is no implementation of either of these randomized kernel PCA algorithms. There are implementations in python (scikit-learn)(needs reference) and a partial [implementation](https://github.com/lopezpaz/randomized_nonlinear_component_analysis/blob/master/code/rca.r) by Lopez-Paz, first author of [5].

# Nyström low-rank approximation
This method is based on the idea that we can randomly subsample our data points and obtain an approximation of the full kernel matrix using the subsampled kernel matrix. The magic is that there are provable bounds on the errors of the approximation and also in most cases a fairly small subsampling size would work very well. In more details, let **X** be our data matrix, consisting of **n** observations (rows) and **d** features (columns). Let **k(x,y)** be the kernel function that we want to use. For the vanilla kernel PCA, we would need to compute the eigendecomposition of the kernel matrix **K**, which has entries

![1](https://user-images.githubusercontent.com/69961386/131489031-b6b0fcc8-a52d-4ec7-be17-cdbaa13e9032.jpg)


As we mentioned above, in practice this matrix will usually be very large. Let **Y** be a submatrix of **X** obtained by randomly sampling **q** rows, i.e. data points. Following [1], we can obtain the following optimal rank-q decomposition of **K** :

![2](https://user-images.githubusercontent.com/69961386/131489192-a22f9618-f078-4be4-bcb2-40306f099442.jpg)


In the decomposition, the first matrix is the submatrix of **K**, only containing the **q** columns, corresponding to our subsampling, and **L** is the kernel matrix of **Y**. Using the formula, we can obtain approximations for the first **q** eigenvalues and eigenvectors of **K** from the eigenvalues and eigenvectors of **L**. Once again, following [1], we can then obtain the following formula for projecting **X** onto the approximated  eigenvectors, thus computing the new feature for each data point:

![3](https://user-images.githubusercontent.com/69961386/131489802-082661cd-14df-470a-a94f-0f09caacc2cc.jpg)

where **A** is the diagonal matrix containing the eigenvalues of **L** and **V** is the matrix with columns the corresponding eigenvectors of **L**. We multiply by the inverse square roots of the eigenvalues in order to normalize the eigenvectors. For details, see [1] and [6].

# Random Fourier features
Here the mathematics is a bit less intuitive on a first glance. Informally, given a kernel **k**, there is a theorem by Bochner, that implies a particular form for the Fourier transform of **k**, namely it happens to be a positive finite Borel measure. That gives us the possibility to sample from it, i.e. to approximate our kernel by a Monte Carlo random sampling. The key part is to find what is the exact formula for the Fourier transform. As it turns out, for the Gaussian kernel, the Fourier transform is given by a multivariate normal distribution with mean 0 and a particular variance. The Fourier transforms of several other kernels have been computed but our implementation only includes the radial basis functions kernel. Our implementation is partially basd on this [code](https://github.com/lopezpaz/randomized_nonlinear_component_analysis/blob/master/code/rca.r) by Lopez-Paz, first author of [5]. For further details, we refer to [3], [4], [5].![4](https://user-images.githubusercontent.com/69961386/131490755-8bba1286-c949-42ed-b263-c1b626d3375b.jpg)


# Our implementation

**Nyström:**
1. Randomly sample a subset of the data points.
1. Compute the partial kernel matrix with one row for each data point in the original set and one column for each data point in the subset.
1. Modify the kernel matrix as if it were representing centered data in feature space.

![4](https://user-images.githubusercontent.com/69961386/131490935-5aad4d04-819e-4726-a759-d675b98433da.jpg)

where ```1_M``` is a matrix with all elements equal to ```1/m```.

1. Compute a truncated SVD for the kernel matrix of the subset of the samples.
1. Output the final transformation.

**Random Fourier features:**
1. Randomly sample from the appropriate distribution (the one corresponding to the kernel in Fourier space).
1. Generate the random features.
1. Center the kernel matrix (same as for the Nyström method).
1. Compute a truncated SVD for the kernel matrix of the subset of the samples.
1. Output the final transformation.


# References (in no particular order)

1. Williams, C. K. I, Seeger, M, *Using the Nyström method to speed up kernel machines*, NIPS 2000.
1. Drineas, P., Mahoney, M. W., *On the Nyström method for approximating a Gram matrix for improved kernel-based learning*, Journal of Machine learning research 6 (2005).
1. Yang, T., Yu-Feng, L., Mahdabi, M., Jing, R., Zhou, Z.- H., *Nyström method vs random Fourier features: a theoritica and empirical comparison*, Advances in Neural Information Processing Systems 25 (NIPS 2012).
1. Rahimi, A., Recht, B, *Random features for large-scale kernel machines*, NIPS'07: Proceedings of the 20th International Conference on Neural Information Processing.
1. Lopez-Paz, D., Sra, S., Smola, A. J., Ghahramani, Z., Schölkopf, B., *Randomized nonlinear component analysis*, Proceedings of the 31st International Conference on Machine Learning, PMLR 32(2):1359-1367, 2014.
1. Schölkopf, B., Smola, A. J., *Nonlinear component analysis as a kernel eigenvalue problem*, July 1998, Neural Computation 10(5):1299-1319.

