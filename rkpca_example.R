
source("./rpca_package/rkpca_main.R")


# Computational performance
random_vector <- rnorm(10000000)

data <- data.frame(matrix(random_vector, nrow = 100000, ncol = 100))

names(data) <- paste0("col_", 1:100)

rpca <- Randomized_kpca$new(
  subsample_size = 1000,
  num_components = 3,
  kernel = "anovadot",
  kernel_degree = 5
)


system.time(
  model <- rpca$rpca_fit(data)
)

microbenchmark::microbenchmark(
  blep <- rpca$rpca_fit(data),
  times = 10
)

# Visual examples
data(iris)
test <- sample(1:150, 20)

##
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



# A vector of labels can be passed for balanced sampling when using the Nyström method.
sampling_labels <- TRUE
if (sampling_labels) {
  
  label_vector <- iris[-test, 5]
  
} else {
  
  label_vector <- c()
}


{
  
  fitted <- rkpca$rpca_fit(
    data = iris[-test,-5],
    labels = label_vector
  )
  
  to_plot <- rkpca$rpca_predict(rpca = fitted, iris[test, -5])
  
  to_plot$label <- factor(iris[test, 5])
  
  plot(to_plot[, 1:2], col = as.integer(to_plot$label))
  
  
  predicted_test <- rkpca$rpca_predict(
    fitted,
    iris[-test, -5]
  )
  
  predicted_test$label <- factor(iris[-test, 5])
  
  points(predicted_test[, 1:2], col = as.numeric(predicted_test$label))
}


###
# Kernel PCA from "kernlab" without randomization. Example taken from
# the "kernlab" documentation.

kpc <- kernlab::kpca(~., data=iris[-test,-5], kernel="rbfdot",
                     kpar=list(sigma=0.2),features=2)

#print the principal component vectors
kernlab::pcv(kpc)

#plot the data projection on the components
plot(kernlab::rotated(kpc),col=as.integer(iris[-test,5]),
     xlab="1st Principal Component", ylab="2nd Principal Component")

#embed remaining points
emb <- kernlab::predict(kpc,iris[test,-5])
points(emb,col=as.integer(iris[test,5]))


