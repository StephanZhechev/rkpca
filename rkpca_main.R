#' Dependencies
#'
#' library(R6)
#' library(kernlab)
#' library(svd)
#' library(magrittr)
#' library(dplyr)


'%>%' <- magrittr::'%>%'

Randomized_kpca <- R6::R6Class(
  
  "Randomized_kpca",
  
  private = list(
    
    compute_principal_components = function(matrix){
      #' Compute truncated principal components normalized by
      #' the singular values.
      #'
      #' @param matrix (matrix / kernlab::kernel_matrix) must be
      #' an object that can be coerced to a matrix.
      #'
      #' @returns list of normalized truncated principal components
      #' and corresponding singular values.
      
      svd <- svd::propack.svd(as.matrix(matrix), neig = self$num_components)
      
      result <- list(
        "principal_components" = svd$v,
        "singular_values" = svd$d
      )
      
      return(result)
      
    },
    
    
    compute_kernel_matrix = function(
      data,
      labels = NULL
    ) {
      #' Compute kernel matrix from the input data.
      #'
      #' @param data (data frame)
      #'
      #' @param labels (vector) labels for the data points.
      #' Will be used if the randomization type is "nystrom"
      #' when drawing a sample from the data.
      
      private$set_sigma(data)
      
      if(self$randomization_type == "nystrom") {
        
        if (!is.null(labels)) {
          
          result <- private$compute_kernel_matrix_nystrom(
            data,
            labels
          )
          
        } else {
          
          result <- private$compute_kernel_matrix_nystrom(data)
        }
      }
      
      
      if(self$randomization_type == "fourier") {
        
        if (!is.null(labels)) {
          
          warning("Labels are not used for Fourier features. Ignoring!")
        }
        
        result <- private$compute_kernel_matrix_fourier(data)
      }
      
      return(result)
    },
    
    
    compute_kernel_matrix_nystrom = function(
      data,
      labels = c()
    ) {
      #' Compute a Nystroem approximation of a kernel matrix.
      #'
      #' @param data (data frame)
      #'
      #' @param labels (vector) labels for the data points.
      #' Will be used if the when drawing a sample from the data
      #' in order to make the sample balanced with respect to the given classes.
      
      data <- as.matrix(data)
      
      if ("sigma" %in% names(self$kernel_params)) {
        
        if (is.na(self$kernel_params$sigma)) {
          
          self$kernel_params$sigma <- unname(kernlab::sigest(data, scaled = NULL)[2])
        }
      }
      
      indices <- get_sample_indices(
        data = data,
        labels = labels,
        sample_size = self$subsample_size
      )
      
      data_sample <- data[indices, ]
      
      # Compute the prediction partial kernel matrix.
      kernel_matrix <- kernlab::kernelMatrix(
        do.call(
          eval(parse(text = self$kernel)),
          self$kernel_params
        ),
        data_sample
      )
      
      result <- list(
        "indices" = indices,
        "kernel_matrix" = kernel_matrix,
        "standard_deviation" = self$kernel_params$sigma
      )
      
      return(result)
    },
    
    
    compute_kernel_matrix_fourier = function(data) {
      #' Approximation of the kernel matrix using random Fourier features.
      #'
      #' @param data (data frame matrix)
      
      fourier_sample <- matrix(
        rnorm(
          ncol(data) * self$subsample_size,
          sd = 2 * self$kernel_params$sigma
        ),
        ncol(data)
      )
      
      b <- runif(self$subsample_size, 0, 2*pi)
      
      f <- function(new_data) {
        new_data <- as.matrix(new_data)
        return(
          new_data %*% fourier_sample + t(matrix(b, self$subsample_size, nrow(new_data)))
        )
      }
      
      result <- list()
      
      result$kernel_matrix <- cos(f(data))
      
      result$standard_deviation <- self$kernel_params$sigma
      
      return(result)
    },
    
    
    center_kernel_matrix = function(mat){
      #' Adjust the kernel matrix as if the underlying
      #' data is column centered.
      #'
      #' @details While it is not computationally feasible
      #' to center the columns of the embedding of our data
      #' in feature space, we can adjust the corresponding
      #' kernel matrix as if the data were column centered:
      #'
      #' K_1 = K - 1_M * K - K %*% 1_M + 1_M * K * 1_M =
      #' = t(t(K - colSums(K)/m) - rowSums(K)/m) + sum(K)/m^2,
      #'
      #' where 1_M is the matrix of the corresponding dimension
      #' having all its entries equal to 1/ncol(K) and I is the
      #' identity matrix, potentially truncated accordingly.
      #'
      #' @param mat (matrix / kernlab::kernel_matrix)
      #' A matrix or an object that is coercible to a matrix.
      #' It doesn't need to be a square matrix.
      
      m <- ncol(mat)
      
      centered_mat <- t(t(mat - colSums(mat)/m) -  rowSums(mat)/m) + sum(mat)/m^2
      
      return(centered_mat)
    },
    
    set_sigma = function(data) {
      #' Estimate the kernel standard deviation from the data.
      #'
      #' @details Setter method that sets the corresponding parameter
      #' into the internal state of the class instance.
      #'
      #' @param data (data frame / matrix) must be coercible to a matrix.
      
      data <- as.matrix(data)
      
      if ("sigma" %in% names(self$kernel_params)) {
        
        if (is.na(self$kernel_params$sigma)) {
          
          self$kernel_params$sigma <- unname(kernlab::sigest(data, scaled = NULL)[2])
        }
      }
    }
    
  ),
  
  public = list(
    
    # Can be "nystrom" or "fourier".
    randomization_type = "nystrom",
    
    subsample_size = 20,
    
    num_components = 2,
    
    kernel = "rbfdot",
    
    kernel_params = list(),
    
    initialize = function(
      
      subsample_size = -1L,
      num_components = 3L,
      randomization_type = "nystrom",
      kernel = NA,
      kernel_degree = 1L,
      kernel_sigma = NA,
      kernel_offset = 1,
      kernel_scale = 1,
      kernel_order = 1
    ) {
      #' @details
      #' 1. If "randomization_type" == "fourier", the "kernel" and
      #' "kernel_params" will be ignored.
      #' 2. If "randomization_type" == "nystrom", the relevant
      #' kernel parameters for the selected kernel will be taken
      #' and the rest will be ignored.
      #'
      #' @param subsample_size (int = -1L) The subsample size.
      #' If set to -1L, standard KPCA will be run without any
      #' subsampling.
      #'
      #' @param num_components (int = 3L) How many principal
      #' components to return. The class uses truncated SVD,
      #' therefore setting the desired number might cut on
      #' computational cost.
      #'
      #' @param randomization_type (string) in c("fourier", "nystrom").
      #' What type of randomization to be used.
      #'
      #' @param kernel (string) in c("rbfdot", "polydot", "tanhdot",
      #' "vanilladot", "laplacedot", "besseldot", "anovadot", "splinedot").
      #' The class allows for all the kernels available in "kernlab".
      #' For information about the kernel parameters, see section "Kernel Functions"
      #' in the "kernlab" documentation.
      #'
      #' @param kernel_degree (int = 1L)
      #'
      #' @param kernel_sigma (float = NA) will be inferred automatically if NA.
      #'
      #' @param kernel_offset (float = 1)
      #'
      #' @param kernel_scale (float = 1)
      #'
      #' @param kernel_order (complex = 1) order of the Bessel function
      #' to be used as a kernel when "kernel" == "besseldot".
      
      if (subsample_size <= 0) {
        
        stop("Please input a valid subsample_size!")
      }
      
      self$subsample_size <- subsample_size
      
      self$randomization_type <- match.arg(randomization_type, c("nystrom", "fourier"))
      
      if (self$randomization_type == "fourier") {
        
        kernel <- "rbfdot"
        
        self$kernel_params <- list(
          
          "sigma" = kernel_sigma
        )
        
      } else {
        
        kernel <- match.arg(
          kernel,
          c("rbfdot", "polydot", "tanhdot", "vanilladot",
            "laplacedot", "besseldot", "anovadot", "splinedot")
        )
        
        self$kernel <- paste0("kernlab::", kernel)
        
        if (num_components > subsample_size) {
          
          warning("Number of components cannot exceed subsample size!")
          
          num_components <- subsample_size
        }
        
        self$num_components <- num_components
        
        self$kernel_params <- switch(
          kernel,
          "rbfdot" = list("sigma" = kernel_sigma),
          "polydot" = list(
            "degree" = kernel_degree,
            "scale" = kernel_scale,
            "offset" = kernel_offset
          ),
          "tanhdot" = list(
            "scale" = kernel_sigma,
            "offset" = kernel_offset
          ),
          "vanilladot" = list(),
          "laplacedot" = list(sigma = kernel_sigma),
          "besseldot" = list(
            "sigma" = kernel_sigma,
            "order" = kernel_order,
            "degree" = kernel_degree
          ),
          "anovadot" = list(
            "degree" = kernel_degree,
            "sigma" = kernel_sigma
          ),
          "splinedot" = list()
        )
      }
    },
    
    
    rpca_fit = function(
      data,
      labels = NULL
    ) {
      #' Fit a randomized kernel PCA.
      #'
      #' @param data (matrix / data frame)
      #'
      #' @param labels (vector) will be used for balanced
      #' subsetting of grouped data in the Nystroem method.
      
      kernel_matrix_elements <- private$compute_kernel_matrix(data, labels)
      
      if (self$randomization_type == "nystrom") {
        
        indices <- kernel_matrix_elements$indices
      }
      
      standard_deviation <- kernel_matrix_elements$standard_deviation
      
      kernel_matrix <- kernel_matrix_elements$kernel_matrix
      
      kernel_matrix <- private$center_kernel_matrix(kernel_matrix)
      
      pc <- private$compute_principal_components(kernel_matrix)
      
      result <- list(
        
        "transformation_matrix" = pc$principal_components %*%
          diag(1/sqrt(pc$singular_values)),
        
        "singular_values" = pc$singular_values,
        
        "randomization_type" = self$randomization_type,
        
        "kernel" = self$kernel,
        
        "kernel_params" = self$kernel_params
      )
      
      if (self$randomization_type == "nystrom"){
        
        result$kernel_prediction_matrix <- data[indices, ]
        
        result$sampled_indices <- indices
      }
      
      return(result)
    },
    
    
    rpca_predict = function(rpca, newdata) {
      #' Predict from an RPCA model.
      #'
      #' @details The prediction will be computed using
      #' the randomization type and kernel parameters (if any)
      #' from the input "rpca" object and not from the class
      #' attributes!
      #'
      #' @param rpca (list) fitted with rpca_fit().
      #'
      #' @param newdata (matrix / data frame).
      
      newdata <- as.matrix(newdata)
      
      if (rpca$randomization_type == "nystrom") {
        
        if (ncol(newdata) != ncol(rpca$kernel_prediction_matrix)) {
          
          stop("newdata must have the same columns as the data rpca was fitted on!")
        }
        
        # Compute the prediction partial kernel matrix.
        prediction_kernel_matrix <- kernlab::kernelMatrix(
          do.call(
            eval(parse(text = rpca$kernel)),
            rpca$kernel_params
          ),
          newdata,
          as.matrix(rpca$kernel_prediction_matrix)
        )
        
        prediction_kernel_matrix <- private$center_kernel_matrix(prediction_kernel_matrix)
        
        new_features <- prediction_kernel_matrix %*% rpca$transformation_matrix
        
      }
      
      if (rpca$randomization_type == "fourier") {
        
        subsample_size <- dim(rpca$transformation_matrix)[1]
        
        w <- matrix(
          rnorm(
            ncol(newdata) * subsample_size,
            sd = 2 * rpca$kernel_params$sigma
          ),
          ncol(newdata)
        )
        
        b <- runif(subsample_size, 0, 2*pi)
        
        f <- function(x0) x0 %*% w + t(matrix(b, subsample_size, nrow(x0)))
        
        mc_simulated_features <- cos(f(newdata))
        
        mc_simulated_features <- private$center_kernel_matrix(mc_simulated_features)
        
        new_features <- mc_simulated_features %*% rpca$transformation_matrix
        
      }
      
      new_features <- data.frame(new_features)
      
      names(new_features) <- paste0("feature_", 1:length(rpca$singular_values))
      
      return(new_features)
    }
    
  )
)





get_sample_indices = function(
  data,
  labels = c(),
  sample_size = -1L
) {
  #' Get sample indices based on input data
  #' and a vector of labels.
  #'
  #' @details The function mostly handles special cases. The indices are
  #' sampled using via calling get_sample_indices_from_labels().
  #'
  #' @param data (data frame / matrix)
  #'
  #' @param labels (vector) must have
  #' length(labels) == nrow(data).
  #'
  #' @param sample_size (int = -1L) if non-positive, will
  #' return the full set of indices.
  #'
  #' @returns vector of indices to be used for sampling.
  
  if (length(labels) == 0) {
    
    if (sample_size > 0 & sample_size < nrow(data)) {
      
      indices <- sample(1:nrow(data), sample_size)
      
    } else{
      
      indices <- 1:nrow(data)
    }
    
  } else {
    
    if (nrow(data) != length(labels)) {
      
      stop("number of labels is different from the number of data points!")
      
    } else {
      
      indices <- get_sample_indices_from_labels(labels, sample_size)
    }
    
  }
  
  return(indices)
}



get_sample_indices_from_labels = function(
  labels,
  sample_size = -1L
) {
  #' Balanced sampling of indices from a vector with categorical
  #' labels.
  #'
  #' @details This is a helper function that will be used for
  #' Nystrom sampling when the data is split into categories.
  #' Rather than sampling the data itself, we simply sample indices
  #' from the vector of labels.
  #' If the self$subsample_size is smaller than the unique
  #' labels, the sampling will be uniform.
  #'
  #' @param labels (vector)
  #'
  #' @param sample_size (int = -1L) if non-positive, will
  #' return the full set of indices.
  #'
  #' @returns vector of indices to be used for sampling.
  
  if (length(labels) == 0) {
    
    stop("Empty labels vector!")
  }
  
  if (sample_size <= 0 | sample_size >= length(labels)) {
    
    warning("Sample size is non-positive or has larger length than labels!")
    return(1:length(labels))
  }
  
  if (length(unique(labels)) > sample_size) {
    
    warning(
      paste("Sample size is smaller than the number of unique groups! \n",
            "Ignoring the labels when sampling!")
    )
    
    return(sample(1:length(labels), size = sample_size))
  }
  
  df <- data.frame(
    "label" = labels,
    "index" = 1:length(labels)
  )
  
  sample_frac <- sample_size / length(labels)
  
  sample <- df %>% dplyr::group_by(label) %>%
    dplyr::sample_frac(size = sample_frac)
  
  return(sample[, "index", drop = TRUE])
}

