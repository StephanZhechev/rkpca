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

# Usage

The main code is packaged into an ```R6``` class.

1. Initialization with the desired attributes. The full documentation can be found in docstrings in the code. Here we provide brief explanations.

    * ```subsample_size``` is self-explanatory.
    * ```randomization_type``` can be ```nystrom``` or ```fourier```.
    * ```num_components``` controls for the number of principal components to be produced. The implementation makes use of an algorithm
      for truncated SVD, so setting this to the desired number rather than computing many components and subsetting, will definitely cut on computational costs.
    * ```kernel``` can be any kernel from the ```kernlab``` library. If ```randomization_type == "fourier"```, only ```rbfdot``` is supported, while
       ```randomization_type == "nystrom"``` supports all the kernels implemented in ```kernlab```.
    * The other arguments control the selected kernel. Only the necessary will be used by the functions later on. For more information, see the [kernlab documentation](https://cran.r-project.org/web/packages/kernlab/kernlab.pdf).
    * If not ```kernel_sigma``` is passed to the constructor, the standard deviation will later be automatically estimated from the data used to fit a randomized kernel PCA using the ```kernlab::sigest``` function.


```
rkpca <- Randomized_kpca$new(
  subsample_size = 50,
  randomization_type = "nystrom",
  num_components = 2,
  kernel = "rbfdot",
  kernel_degree = 3,
  kernel_sigma = NA,
  kernel_offset = 1,
  kernel_scale = 1,
  kernel_order = 1
)
```

1. Fit a randomized kernel PCA on a given data frame or matrix. The ```rpca_fit()``` method works does not take a formula as an argument. That means that the kernel PCA will be estimated for all the columns that are passed.

    * ```data``` should ideally be a data frame but any object coercible to a matrix should do.
    * ```labels``` is an optional argument that can take a vector of labels of length ```nrow(data)```, that will be used for balanced subsetting if the Nyström method is selected. Passing labels when using random Fourier features does not affect the model, which will also be announced in a warning.

```
fitted <- rkpca$rpca_fit(
    data = input_data,
    labels = label_vector
  )
```

1. Predict from a fitted randomized kernel PCA object and new data.
 
    * ```rpca``` should be an object produced by the ```rpca_fit()``` method.
    * ```newdata``` is new data to be projected onto the principal components. It must have the same number of columns as the original data, ideally with the same names.


```
  predicted <- rkpca$rpca_predict(
    rpca = fitted,
    newdata = some_new_data
  )
```


# Performance

This code was run on a laptop with the i5-8265U CPU @ 1.60GHz, 1800 Mhz, 4 Core(s), 8 Logical Processor(s) and 16GB of RAM.

![image](https://user-images.githubusercontent.com/69961386/131563172-e9de4258-5a68-4af7-b4f9-19e8dcc1eac4.png)

Running full kernel PCA on 100k rows would probably take more than 200GB of RAM and would be completely impossible on this laptop.

# License...or something

I wrote this code as part of my work for [DataScience Service GmbH](https://www.datascience-service.at/). However, the code is **completely open source** and comes with no waranty or anything. Feel free to use it for any purpose you like, make copies etc.

# Bugs

Please feel free to submit pull requests with bug fixes and improvements.

# References (in no particular order)

1. Williams, C. K. I, Seeger, M, *Using the Nyström method to speed up kernel machines*, NIPS 2000.
1. Drineas, P., Mahoney, M. W., *On the Nyström method for approximating a Gram matrix for improved kernel-based learning*, Journal of Machine learning research 6 (2005).
1. Yang, T., Yu-Feng, L., Mahdabi, M., Jing, R., Zhou, Z.- H., *Nyström method vs random Fourier features: a theoritica and empirical comparison*, Advances in Neural Information Processing Systems 25 (NIPS 2012).
1. Rahimi, A., Recht, B, *Random features for large-scale kernel machines*, NIPS'07: Proceedings of the 20th International Conference on Neural Information Processing.
1. Lopez-Paz, D., Sra, S., Smola, A. J., Ghahramani, Z., Schölkopf, B., *Randomized nonlinear component analysis*, Proceedings of the 31st International Conference on Machine Learning, PMLR 32(2):1359-1367, 2014.
1. Schölkopf, B., Smola, A. J., *Nonlinear component analysis as a kernel eigenvalue problem*, July 1998, Neural Computation 10(5):1299-1319.

