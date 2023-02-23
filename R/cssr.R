# TODO(gregfaletto): implement protolasso and clusterRepLasso (located in 
# toy_ex_slide_funcs.R, in /Users/gregfaletto/Google Drive/Data Science/LaTeX/Generalized Stability Selection Presentation)

# TODO(gregfaletto): make sure behavior of functions makes sense if selection
# indicator matrix ends up containing all 0s or all 1s (and maybe throw a
# warning or error in this case)

### BELOW IS DONE AND IN RMD FILE

#' Cluster Stability Selection
#'
#' Executes cluster stability selection algorithm. Takes subsamples of data,
#' executes feature selection algorithm on each subsample, and returns matrices
#' of feature selection indicators as well as cluster selection indicators.
#'
#' @param X An n x p numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' p >= 2 features/predictors.
#' @param y The response; can be anything that takes the form of an
#' n-dimensional vector, with the ith entry corresponding to the ith row of X.
#' Typically (and for default fitfun = cssLasso), y should be an n-dimensional
#' numeric vector.
#' @param lambda A tuning parameter or set of tuning parameters that may be used
#' by the feature selection method fitfun. In the default case when
#' fitfun = cssLasso, lambda should be a numeric: the penalty to use for each
#' lasso fit. (css does not require lambda to be any particular object because
#' for a user-specified feature selection method fitfun, lambda can be an
#' arbitrary object. See the description of fitfun below.)
#' @param clusters A list of integer vectors; each vector should contain the 
#' indices of a cluster of features (a subset of 1:p). (If there is only one
#' cluster, clusters can either be a list of length 1 or an integer vector.)
#' All of the provided clusters must be non-overlapping. Every feature not
#' appearing in any cluster will be assumed to be unclustered (that is, they
#' will be treated as if they are in a "cluster" containing only themselves). If
#' clusters is a list of length 0 (or a list only containing clusters of length
#' 1), then css() returns the same results as stability selection (so the
#' returned feat_sel_mat will be identical to clus_sel_mat). Names for the
#' clusters will be needed later; any clusters that are not given names in the
#' provided list will be given names automatically by css. Default is list() (so
#' no clusters are specified).
#' @param fitfun A function; the feature selection function used on each
#' subsample by cluster stability selection. This can be any feature selection
#' method; the only requirement is that it accepts the arguments (and only the
#' arguments) X, y, and lambda and returns an integer vector that is a subset of
#' 1:p. For example, fitfun could be best subset selection or forward stepwise
#' selection or LARS and lambda could be the desired model size; or fitfun could be the
#' elastic net and lambda could be a length-two vector specifying lambda and
#' alpha. Default is cssLasso, an implementation of lasso (relying on the R
#' package glmnet), where lambda must be a positive numeric specifying the L1
#' penalty for the lasso.
#' @param sampling_type A character vector; either "SS" or "MB". For "MB",
#' all B subsamples are drawn randomly (as proposed by Meinshausen and Bühlmann
#' 2010). For "SS", in addition to these B subsamples, the B complementary pair
#' subsamples will be drawn as well (see Faletto and Bien 2022 or Shah and
#' Samworth 2013 for details). Default is "SS", and "MB" is not supported yet.
#' @param B Integer or numeric; the number of subsamples. Note: For
#' sampling.type=="MB" the total number of subsamples will be `B`; for
#' sampling_type="SS" the number of subsamples will be `2*B`. Default is 100
#' for sampling_type="MB" and 50 for sampling_type="SS".
#' @param prop_feats_remove Numeric; if prop_feats_remove is greater than 0,
#' then on each subsample, each feature is randomly dropped from the design
#' matrix that is provided to fitfun with probability prop_feats_remove
#' (independently across features). That is, in a typical subsample,
#' prop_feats_remove*p features will be dropped (though this number will vary).
#' This is similar in spirit (but distinct from) extended stability selection
#' (Beinrucker et. al. 2016); see their paper for some of the benefits of
#' dropping features (besides increasing computational speed and decreasing
#' memory requirements). For sampling_type="SS", the features dropped in
#' each complementary pair of subsamples are identical in order to ensure that
#' the theoretical guarantees of Faletto and Bien (2022) are retained within
#' each individual pair of subsamples. (Note that this feature is not
#' investigated either theoretically or in simulations by Faletto and Bien
#' 2022). Must be between 0 and 1. Default is 0.
#' @param train_inds Optional; an integer or numeric vector containing the
#' indices of observations in X and y to set aside for model training by the
#' function getCssPreds after feature selection. (This will only work if y is
#' real-valued, because getCssPreds using ordinary least squares regression to
#' generate predictions.) If train_inds is not provided, all of the observations
#' in the provided data set will be used for feature selection.
#' @param num_cores Optional; an integer. If using parallel processing, the
#' number of cores to use for parallel processing (num_cores will be supplied
#' internally as the mc.cores argument of parallel::mclapply).
#' @return A list containing the following items: \item{feat_sel_mat}{A B (or
#' `2*B` for sampling.method "SS") x p numeric (binary) matrix.
#' `feat_sel_mat[i, j] = 1` if feature j was selected by the base feature
#' selection method on subsample i, and 0 otherwise.} \item{clus_sel_mat}{A B
#' (or 2*B for SS sampling) x length(clusters) numeric (binary) matrix.
#' `clus_sel_mat[i, j] = 1` if at least one feature from cluster j was selected
#' by the base feature selection method on subsample i, and 0 otherwise.}
#' \item{X}{The X matrix provided to css, coerced from a data.frame to a matrix
#' if needed.} \item{y}{The y vector provided to css.} \item{clusters}{A named
#' list of integer vectors containing all of the clusters provided to css, as
#' well as size 1 clusters of any features not listed in any of the clusters
#' provided to css. All clusters will have names; any clusters not provided with
#' a name in the input to css will be given names automatically by css (of the
#' form c1, etc.).} \item{train_inds}{Identical to the train_inds provided to
#' css.}
#' @author Gregory Faletto, Jacob Bien
#' @references Faletto, G., & Bien, J. (2022). Cluster Stability Selection.
#' \emph{arXiv preprint arXiv:2201.00494}.
#' \url{https://arxiv.org/abs/2201.00494}. \cr Shah, R. D., & Samworth, R. J.
#' (2013). Variable selection with error control: Another look at stability
#' selection. \emph{Journal of the Royal Statistical Society. Series B:
#' Statistical Methodology}, 75(1), 55–80.
#' \url{https://doi.org/10.1109/RITA.2014.2302071}. \cr Meinshausen, N., &
#' Bühlmann, P. (2010). Stability Selection. \emph{Journal of the Royal
#' Statistical Society. Series B: Statistical Methodology}, 72(4), 417–473.
#' \url{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2010.00740.x}.
#' \cr Beinrucker, A., Dogan, Ü., &
#' Blanchard, G. (2016). Extensions of stability selection using subsamples of
#' observations and covariates. \emph{Statistics and Computing}, 26(5), 1059-
#' 1077. \url{https://doi.org/10.1007/s11222-015-9589-y}. \cr Jerome Friedman,
#' Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized
#' Linear Models via Coordinate Descent. \emph{Journal of Statistical Software},
#' 33(1), 1-22. URL \url{https://www.jstatsoft.org/v33/i01/}.
#' @export
css <- function(X, y, lambda, clusters = list(), fitfun = cssLasso,
    sampling_type = "SS", B = ifelse(sampling_type == "MB", 100L, 50L),
    prop_feats_remove = 0, train_inds = integer(), num_cores = 1L
    ){

    # Check inputs

    check_list <- checkCssInputs(X, y, lambda, clusters, fitfun, sampling_type,
        B, prop_feats_remove, train_inds, num_cores)

    feat_names <- check_list$feat_names
    X <- check_list$X
    clusters <- check_list$clusters

    rm(check_list)

    n <- nrow(X)
    p <- ncol(X)

    train_inds <- as.integer(train_inds)

    ### Create subsamples

    sel_inds <- setdiff(1:n, train_inds)
    n_sel <- length(sel_inds)
    if(n_sel < 4){
        stop("Too many training indices provided (must be at least 4 observations left for feature selection, and ideally many more)")
    }

    subsamps_object <- createSubsamples(n_sel, p, B, sampling_type,
        prop_feats_remove)

    ### Get matrix of selected feature sets from subsamples

    stopifnot(!is.matrix(y))

    feat_sel_mat <- getSelMatrix(X[sel_inds, ], y[sel_inds], lambda, B,
        sampling_type, subsamps_object, num_cores, fitfun)

    if(any(!is.na(feat_names))){
        colnames(feat_sel_mat) <- feat_names
        colnames(X) <- feat_names
    }

    ### Get selection proportions for clusters corresponding to each feature

    clus_sel_mat <- getClusterSelMatrix(clusters, feat_sel_mat)

    # Check outputs
    stopifnot(!is.null(colnames(clus_sel_mat)))
    stopifnot(all(colnames(clus_sel_mat) == names(clusters)))

    ret <- list(feat_sel_mat = feat_sel_mat,
        clus_sel_mat = clus_sel_mat,
        X = X,
        y = y,
        clusters = clusters,
        train_inds = train_inds
        )

    class(ret) <- "cssr"

    return(ret)
}

### BELOW IS DONE AND IN RMD FILE

# TODO(gregfaletto) change cluster_size into a vector of sizes (maybe also
# deprecate n_clusters as an input, since this would be inferred by the length
# of cluster_sizes?)

#' Generate randomly sampled data including noisy observations of latent
#' variables
#'
#' Generate a data set including latent features Z, observed features X (which
#' may include noisy or noiseless observations of the latent features in Z),
#' an observed response y which is a linear model of features from Z and X as
#' well as independent mean zero noise, and mu (the responses from y without
#' the added noise). Data is generated in the same way as in the simulations
#' from Faletto and Bien (2022).
#' @param n Integer or numeric; the number of observations to generate. (The
#' generated X and Z will have n rows, and the generated y and mu will have
#' length n.)
#' @param p Integer or numeric; the number of features to generate. The
#' generated X will have p columns.
#' @param k_unclustered Integer or numeric; the number of features in X that
#' will have nonzero coefficients in the true model for y among those features 
#' not generated from the n_clusters latent variables (called "weak signal" 
#' features in the simulations from Faletto and Bien 2022). The coefficients on
#' these features will be determined by beta_unclustered. Must be at least 1.
#' @param cluster_size Integer or numeric; for each of the n_clusters latent
#' variables, X will contain cluster_size noisy proxies that are correlated with
#' the latent variable. Must be at least 2.
#' @param n_clusters Integer or numeric; the number of latent variables to
#' generate, each of which will be associated with an observed cluster in X.
#' Must be at least 1. Default is 1.
#' @param sig_clusters Integer or numeric; the number of generated latent
#' features that will have nonzero coefficients in the true model for y (all of
#' them will have coefficient beta_latent). Must be less than or equal to
#' n_clusters. Default is 1.
#' @param rho Integer or numeric; the correlation of the proxies in each cluster
#' with the latent variable. Must be greater than 0. Default is 0.9.
#' @param beta_latent Integer or numeric; the coefficient used for all
#' sig_clusters latent variables that have nonzero coefficients in the true
#' model for y. Can't equal 0. Default is 1.5.
#' @param beta_unclustered Integer or numeric; the maximum coefficient in the
#' model for y among the k_unclustered features in X not generated from the
#' latent variables. The coefficients of the features will be
#' beta_unclustered/sqrt(1:k_unclustered). Can't equal 0. Default is 1.
#' @param snr Integer or numeric; the signal-to-noise ratio of the response
#' y. If sigma_eps_sq is not specified, the variance of the noise in y will be
#' calculated using the formula sigma_eps_sq = sum(mu^2)/(n * snr). Only one of
#' snr and sigma_eps_sq must be specified. Default is NA.
#' @param sigma_eps_sq Integer or numeric; the variance on the noise added
#' to y. Only one of snr and sigma_eps_sq must be specified. Default is NA.
#' @return A list of the following elements. \item{X}{An n x p numeric matrix of
#' n observations from a p-dimensional multivariate normal distribution
#' generated using the specified parameters. The first n_clusters times
#' cluster_size features will be the clusters of features correlated with the
#' n_clusters latent variables. The next k_unclustered features will be the
#' "weak signal" features, and the remaining p - n_clusters*cluster_size -
#' k_unclustered features will be the unclustered noise features.} \item{y}{A
#' length n numeric vector; the response generated from X, the latent features
#' from Z, and the coefficient vector, along with additive noise.} \item{Z}{The
#' latent features; either a numeric vector (if n_clusters > 1) or a numeric
#' matrix (if n_clusters > 1). Note that (X, Z) is multivariate Gaussian.}
#' item{mu}{A length `n` numeric vector; the expected response given X, Z, and
#' the true coefficient vector (equal to y minus the added noise).}
#' @author Gregory Faletto, Jacob Bien
#' @references
<<faletto2022>>
#' @export
genClusteredData <- function(n, p, k_unclustered, cluster_size, n_clusters=1,
    sig_clusters=1, rho=0.9, beta_latent=1.5, beta_unclustered=1,
    snr=as.numeric(NA), sigma_eps_sq=as.numeric(NA)){

    # Check inputs
    checkgenClusteredDataInputs(p, k_unclustered, cluster_size, n_clusters,
        sig_clusters, rho, beta_latent, beta_unclustered, snr,
        sigma_eps_sq)

    ret <- genZmuY(n=n, p=p, k_unclustered=k_unclustered,
        cluster_size=cluster_size, n_clusters=n_clusters,
        sig_clusters=sig_clusters, beta_latent=beta_latent,
        beta_unclustered=beta_unclustered, snr=snr, sigma_eps_sq=sigma_eps_sq)

    Z <- ret$Z
    y <- ret$y
    mu <- ret$mu
    other_X <- ret$other_X

    # Finally, generate clusters of proxies to complete X. First, get needed
    # variances of noise to add
    noise_var <- getNoiseVar(rho)

    # Generate these noise features
    noise_mat <- matrix(stats::rnorm(n*n_clusters*cluster_size, mean=0,
        sd=sqrt(noise_var)), n, n_clusters*cluster_size)

    # Create matrix of proxies
    proxy_mat <- matrix(as.numeric(NA), n, n_clusters*cluster_size)
    if(n_clusters > 1){
        for(i in 1:n_clusters){
            first_ind <- (i - 1)*cluster_size + 1
            last_ind <- i*cluster_size
            proxy_mat[, first_ind:last_ind] <- Z[, i] +
                noise_mat[, first_ind:last_ind]
        }
    } else{
        stopifnot(ncol(noise_mat) == cluster_size)
        proxy_mat[, 1:cluster_size] <- Z + noise_mat
    }

    X <- cbind(proxy_mat, other_X)
    Z <- as.matrix(Z)

    # Check output
    stopifnot(length(mu) == n)

    stopifnot(nrow(X) == n)
    stopifnot(ncol(X) == p)

    if(any(!is.na(Z))){
        stopifnot(nrow(Z) == n)
        stopifnot(ncol(Z) == n_clusters)
    }

    return(list(X=X, y=y, Z=Z, mu=mu))
}

#' Get variance of noise to add to Z in order to yield proxies X with desired
#' correlations with Z
#'
#' @param cor A numeric vector of desired correlations for each proxy to have
#' with Z. Note: correlations must be positive.
#' @return A vector of variances of independent Gaussian random variables to add
#' to Z in order to yield proxies with the desired correlations with Z.
#' @author Gregory Faletto, Jacob Bien
#' @export
getNoiseVar <- function(cor){
    # Correlation between standard normal Z and X = Z + epsilon where epsilon
    # is normal, independent of Z, and has mean 0 and variance sig_eps_sq:
    # 
    #   E[Z X]/sqrt{Var(Z) Var(X)}
    # = (E[Z^2] + E[Z*epsilon])/sqrt{1*(1 + sig_eps_sq)}
    # = (1 + 0)/sqrt{1 + sig_eps_sq}
    #
    # So we have
    #                 cor = 1/sqrt{1 + sig_eps_sq}
    # \iff 1 + sig_eps_sq = 1/cor^2
    # \iff     sig_eps_sq = 1/cor^2 - 1
    stopifnot(is.numeric(cor) | is.integer(cor))
    stopifnot(all(!is.na(cor)))
    stopifnot(length(cor) >= 1)
    stopifnot(all(cor > 0))
    stopifnot(all(cor <= 1))
    return(1/cor^2 - 1)
}

# #' Get variance of noise to add to Z in order to yield proxies X with desired
# #' correlations between the proxies
# #'
# #' @param cor A numeric vector of desired correlations for each proxy to have
# #' with the other proxies in the cluster. Note: correlations must be positive.
# #' @return A vector of variances of independent Gaussian random variables to add
# #' to Z in order to yield proxies with the desired correlations with Z.
# #' @author Gregory Faletto, Jacob Bien
# getNoiseVar <- function(cor){
#     # Correlation between X_1 = Z + epsilon_1 and X_2 = Z + epsilon_2 where
#     # Z is standard normal and epsilon_i are both normal, independent of Z,
#     # and have mean 0 and variance sig_eps_sq:
#     # 
#     #   E[X_1 X_2]/sqrt{Var(X_1) Var(X_2)}
#     # = (E[Z^2] + E[Z*epsilon_1] + E[Z*epsilon_2] + E[epsilon_1*epsilon_2])/
#     #   sqrt{(1 + sig_eps_sq)^2}
#     # = (1 + 0 + 0 + 0)/(1 + sig_eps_sq)
#     #
#     # So we have
#     #             cor = 1/(1 + sig_eps_sq)
#     # \iff sig_eps_sq = 1/cor - 1
#     stopifnot(is.numeric(cor) | is.integer(cor))
#     stopifnot(all(!is.na(cor)))
#     stopifnot(length(cor) >= 1)
#     stopifnot(all(cor > 0))
#     stopifnot(all(cor <= 1))
#     return(1/cor - 1)
# }

### BELOW IS DONE AND IN RMD FILE


#' Check inputs to genClusteredData
#'
#' @param p Integer or numeric; the number of features to generate. The
#' generated X will have p columns.
#' @param k_unclustered Integer or numeric; the number of features in X that
#' will have nonzero coefficients in the true model for y among those features 
#' not generated from the n_clusters latent variables (called "weak signal" 
#' features in the simulations from Faletto and Bien 2022). The coefficients on
#' these features will be determined by beta_unclustered.
#' @param cluster_size Integer or numeric; for each of the n_clusters latent
#' variables, X will contain cluster_size noisy proxies that are correlated with
#' the latent variable.
#' @param n_clusters Integer or numeric; the number of latent variables to
#' generate, each of which will be associated with an observed cluster in X.
#' Must be at least 1. Default is 1.
#' @param sig_clusters Integer or numeric; the number of generated latent
#' features that will have nonzero coefficients in the true model for y (all of
#' them will have coefficient beta_latent). Must be less than or equal to
#' n_clusters. Default is 1.
#' @param rho Integer or numeric; the covariance of the proxies in each cluster
#' with the latent variable (and each other). Note that the correlation between
#' the features in the cluster will be rho/var. Can't equal 0. Default is 0.9.
#' @param beta_latent Integer or numeric; the coefficient used for all
#' sig_clusters latent variables that have nonzero coefficients in the true
#' model for y. Can't equal 0. Default is 1.5.
#' @param beta_unclustered Integer or numeric; the maximum coefficient in the
#' model for y among the k_unclustered features in X not generated from the
#' latent variables. The coefficients of the features will be
#' beta_unclustered/sqrt(1:k_unclustered). Can't equal 0. Default is 1.
#' @param snr Integer or numeric; the signal-to-noise ratio of the response
#' y. If sigma_eps_sq is not specified, the variance of the noise in y will be
#' calculated using the formula sigma_eps_sq = sum(mu^2)/(n * snr). Only one of
#' snr and sigma_eps_sq must be specified. Default is NA.
#' @param sigma_eps_sq Integer or numeric; the variance on the noise added
#' to y. Only one of snr and sigma_eps_sq must be specified. Default is NA.
#' @return A list of the following elements. \item{X}{An n x p numeric matrix of
#' n observations from a p-dimensional multivariate normal distribution
#' generated using the specified parameters. The first n_clusters times
#' cluster_size features will be the clusters of features correlated with the
#' n_clusters latent variables. The next k_unclustered features will be the
#' "weak signal" features, and the remaining p - n_clusters*cluster_size -
#' k_unclustered features will be the unclustered noise features.} \item{y}{A
#' length n numeric vector; the response generated from X, the latent features
#' from Z, and the coefficient vector, along with additive noise.} \item{Z}{The
#' latent features; either a numeric vector (if n_clusters > 1) or a numeric
#' matrix (if n_clusters > 1). Note that (X, Z) is multivariate Gaussian.}
#' item{mu}{A length `n` numeric vector; the expected response given X, Z, and
#' the true coefficient vector (equal to y minus the added noise).}
#' @author Gregory Faletto, Jacob Bien
checkGenClusteredDataInputs <- function(p, k_unclustered, cluster_size,
    n_clusters, sig_clusters, rho, beta_latent, beta_unclustered, snr,
    sigma_eps_sq){

    stopifnot(is.numeric(sig_clusters) | is.integer(sig_clusters))
    stopifnot(sig_clusters <= n_clusters)
    stopifnot(sig_clusters >= 0)
    stopifnot(sig_clusters == round(sig_clusters))
    
    stopifnot(is.numeric(n_clusters) | is.integer(n_clusters))
    stopifnot(n_clusters == round(n_clusters))
    # TODO(gregfaletto): is it easy to remove the requirement that n_clusters is
    # at least 1 (so that it's possible to generate data with no latent 
    # features)? If so, should only check that cluster_size >= 1 if n_clusters
    # >= 1, and in makeCovarianceMatrix function only need block_size >= 1
    # rather than 2.
    stopifnot(n_clusters >= 1)

    stopifnot(cluster_size >= 2)

    stopifnot(rho > 0)

    stopifnot(beta_latent != 0)
    stopifnot(beta_unclustered != 0)

    stopifnot(is.numeric(k_unclustered) | is.integer(k_unclustered))
    stopifnot(k_unclustered >= 2)
    stopifnot(k_unclustered == round(k_unclustered))

    stopifnot(p >= n_clusters*cluster_size + k_unclustered)

    # Same as make_sparse_blocked_linear_model_random, but ith coefficient
    # of weak signal features is beta_unclustered/sqrt(i) in order to have
    # a definitive ranking of weak signal features.
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }

    if(is.na(snr)){
        stopifnot(all(!is.na(sigma_eps_sq)))
        stopifnot(is.numeric(sigma_eps_sq) | is.integer(sigma_eps_sq))
        stopifnot(length(sigma_eps_sq) == 1)
        stopifnot(sigma_eps_sq >= 0)
    } else{
        stopifnot(is.numeric(snr) | is.integer(snr))
        stopifnot(length(snr) == 1)
        stopifnot(snr > 0)
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Get lambda value for lasso
#'
#' Chooses a lambda value for the lasso used on a subsample of size n/2 (as in
#' cluster stability selection) by cross-validation.
#' @param X An n x p numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the p >= 2 features/predictors that will be used by cluster stability
#' selection.
#' @param y The response; an n-dimensional numeric or integer vector. (Unlike
#' in the more general css setup, this response must be real-valued since
#' lambda will be determined using the lasso with cross-validation.)
#' @param lambda_choice Character; either "min" or "1se". If "min", chooses
#' the lambda that minimizes the cross-validated error; if "1se", chooses the
#' largest lambda within one standard error of the minimum error lambda
#' (resulting in a smaller selected set, which may be desirable because the
#' model size corresponding to the minimum error lambda tends to be larger
#' than optimal. See, for example, Bühlmann and Meinshausen 2006, Prop. 1 and
#' Bühlmann and van de Geer 2011, Section 2.5.1.). Default is "1se".
#' @param nfolds Numeric or integer; the number of folds for cross-validation.
#' Must be at least 4 (as specified by cv.glmnet). Default is 10.
#' @return A numeric; the selected value of lambda.
#' @author Gregory Faletto, Jacob Bien
#' @references Bühlmann, P., & Meinshausen, N. (2006). High-Dimensional Graphs
#' and Variable Selection With the Lasso. \emph{The Annals of Statistics},
#' 34(3), 1436–1462. \url{https://doi.org/10.1214/009053606000000281}.
#' \cr Peter Bühlmann and Sara van de Geer. Statistics for High-Dimensional
#' Data. \emph{Springer Series in Statistics}. Springer, Heidelberg, 2011. ISBN
#' 978-3-642-20191-2. \url{http://dx.doi.org/10.1007/978-3-642-20192-9}. \cr
#' Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization
#' Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of
#' Statistical Software}, 33(1), 1-22. URL \url{https://www.jstatsoft.org/v33/i01/}.
#' @export
getLassoLambda <- function(X, y, lambda_choice="1se", nfolds=10){
    stopifnot(is.character(lambda_choice))
    stopifnot(length(lambda_choice) == 1)
    stopifnot(!is.na(lambda_choice))
    stopifnot(lambda_choice %in% c("min", "1se"))

    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X) | is.integer(X))
    n <- nrow(X)

    stopifnot(is.numeric(nfolds) | is.integer(nfolds))
    stopifnot(length(nfolds) == 1)
    stopifnot(nfolds == round(nfolds))
    stopifnot(nfolds > 3)

    # Since we are using the lasso, we require y to be a real-valued response
    # (unlike for the general cluster stability selection procedure, where y
    # can have a more general format as long as a suitable feature selection
    # function is provided by the user)
    stopifnot(is.numeric(y) | is.integer(y))
    stopifnot(n == length(y))

    # Sample size to use: inflate n/2 by a factor of nfolds/(nfolds - 1),
    # so that each individual lasso fit is of size floor(n/2)

    n_sample <- min(round(n/2*nfolds/(nfolds - 1)), n)
    nfolds <- min(nfolds, n_sample)

    inds_size <- sample(1:n, n_sample)
    size_results <- glmnet::cv.glmnet(x=X[inds_size, ], y=y[inds_size],
        family="gaussian", nfolds=nfolds)

    lambda_ret <- size_results[[paste("lambda.", lambda_choice, sep="")]]

    # Check output
    stopifnot(length(lambda_ret) == 1)
    stopifnot(is.numeric(lambda_ret) | is.integer(lambda_ret))
    stopifnot(lambda_ret >= 0)

    return(lambda_ret)
}

### BELOW IS DONE AND IN RMD FILE

#' Fit model and generate predictions from new data
#'
#' Generate predictions on test data using cluster stability-selected model.
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param testX A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to generate predictions. Must contain the same
#' features (in the same number of columns) as the matrix provided to css, and
#' if the columns of testX are labeled, the names must match the variable names
#' provided to css.
#' @param weighting Character; determines how to calculate the weights to
#' combine features from the selected clusters into weighted averages, called
#' cluster representatives. Must be one of "sparse", "weighted_avg", or
#' "simple_avg'. For "sparse", all the weight is put on the most frequently
#' selected individual cluster member (or divided equally among all the clusters
#' that are tied for the top selection proportion if there is a tie). For
#' "weighted_avg", the weight used for each cluster member is calculated in
#' proportion to the individual selection proportions of each feature. For
#' "simple_avg", each cluster member gets equal weight regardless of the
#' individual feature selection proportions (that is, the cluster representative
#' is just a simple average of all the cluster members). See Faletto and Bien
#' (2022) for details. Default is "weighted_avg".
#' @param cutoff Numeric; getCssPreds will make use only of those clusters with
#' selection proportions equal to at least cutoff. Must be between 0 and 1.
#' Default is 0 (in which case either all clusters are used, or max_num_clusts
#' are used, if max_num_clusts is specified).
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.) Default is 1.
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' max_num_clusts is ignored).
#' @param trainX A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to estimate the linear model from the selected
#' clusters. trainX is only necessary to provide if no train_inds were
#' designated in the css function call to set aside observations for model
#' estimation (though even if train_inds was provided, trainX and trianY will be
#' used for model estimation if they are both provided to getCssPreds). Must 
#' contain the same features (in the same number of columns) as the matrix 
#' provided to css, and if the columns of trainX are labeled, the names must
#' match the variable names provided to css. Default is NA (in which case
#' getCssPreds uses the observations from the train_inds that were provided to
#' css to estimate a linear model).
#' @param trainY The response corresponding to trainX. Must be a real-valued
#' response (unlike in the general css setup) because predictions will be
#' generated by an ordinary least squares model. Must have the same length as
#' the number of rows of trainX. Like trainX, only needs to be provided if no
#' observations were set aside for model estimation by the parameter train_inds
#' in the css function call. Default is NA (in which case getCssPreds uses the
#' observations from the train_inds that were provided to css).
#' @return A vector of predictions corresponding to the observations from testX.
#' @author Gregory Faletto, Jacob Bien
#' @references Faletto, G., & Bien, J. (2022). Cluster Stability Selection.
#' \emph{arXiv preprint arXiv:2201.00494}.
#' \url{https://arxiv.org/abs/2201.00494}.
#' @export
getCssPreds <- function(css_results, testX, weighting="weighted_avg", cutoff=0,
    min_num_clusts=1, max_num_clusts=NA, trainX=NA, trainY=NA){
    # TODO(gregfaletto) Consider adding an argument for a user-provided prediction
    # function in order to allow for more general kinds of predictions than
    # OLS.

    # Check inputs
    
    check_list <- checkGetCssPredsInputs(css_results, testX, weighting, cutoff,
        min_num_clusts, max_num_clusts, trainX, trainY)

    trainXProvided <- check_list$trainXProvided
    trainX <- check_list$trainX
    testX <- check_list$testX
    feat_names <- check_list$feat_names
    max_num_clusts <- check_list$max_num_clusts

    rm(check_list)

    n_train <- nrow(trainX)
    n <- nrow(testX)
    p <- ncol(testX)

    # Take provided training design matrix and testX and turn them into
    # matrices of cluster representatives using information from css_results
    if(trainXProvided){
        train_X_clusters <- formCssDesign(css_results, weighting, cutoff,
            min_num_clusts, max_num_clusts, newx=trainX)
        if(!is.numeric(trainY) & !is.integer(trainY)){
            stop("The provided trainY must be real-valued, because predictions will be generated by ordinary least squares regression.")
        }
        y_train <- trainY
    } else{
        train_X_clusters <- formCssDesign(css_results, weighting, cutoff,
            min_num_clusts, max_num_clusts)
        y_train <- css_results$y[css_results$train_inds]
        if(!is.numeric(y_train) & !is.integer(y_train)){
            stop("Can't generated predictions from the data that was provided to css because the provided y was not real-valued (getCssPreds generated predictions using ordinary least squares regression).")
        }
    }

    stopifnot(length(y_train) == nrow(train_X_clusters))

    testX_clusters <- formCssDesign(css_results, weighting, cutoff,
        min_num_clusts, max_num_clusts, newx=testX)

    stopifnot(ncol(testX_clusters) == ncol(train_X_clusters))

    # Get names for clusters
    clust_X_names <- paste("c_fit_", 1:ncol(testX_clusters), sep="")
    if(!is.null(colnames(train_X_clusters))){
        stopifnot(identical(colnames(train_X_clusters), colnames(testX_clusters)))
        clust_X_names <- colnames(train_X_clusters)
    }

    # Fit linear model on training data via OLS
    if(nrow(train_X_clusters) < ncol(train_X_clusters)){
        err_mess <- paste("css not provided with enough indices to fit OLS model for predictions (number of training indices: ",
            nrow(train_X_clusters), ", number of clusters: ",
            ncol(train_X_clusters),
            "). Try reducing number of clusters by increasing cutoff, or re-run css with a larger number of training indices.",
            sep="")
        stop(err_mess)
    }

    df <- data.frame(y=y_train, train_X_clusters)
    colnames(df)[2:ncol(df)] <- clust_X_names
    model <- stats::lm(y ~., data=df)

    # Use fitted model to generate predictions on testX
    df_test <- data.frame(testX_clusters)
    colnames(df_test) <- clust_X_names
    predictions <- stats::predict.lm(model, newdata=df_test)
    names(predictions) <- NULL

    # Check output
    stopifnot(is.numeric(predictions) | is.integer(predictions))
    stopifnot(length(predictions) == n)
    stopifnot(all(!is.na(predictions)))

    return(predictions)
}

### BELOW IS DONE AND IN RMD FILE

#' Obtain a selected set of clusters and features
#'
#' Generate sets of selected clusters and features from cluster stability
#' selection.
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param weighting Character; determines how to calculate the weights for
#' individual features within the selected clusters. Only those features with
#' nonzero weight within the selected clusters will be returned. Must be one of
#' "sparse", "weighted_avg", or "simple_avg'. For "sparse", all the weight is
#' put on the most frequently
#' selected individual cluster member (or divided equally among all the clusters
#' that are tied for the top selection proportion if there is a tie). For
#' "weighted_avg", only the features within a selected cluster that were
#' themselves selected on at least one subsample will have nonzero weight. For
#' "simple_avg", each cluster member gets equal weight regardless of the
#' individual feature selection proportions (that is, all cluster members within
#' each selected cluster will be returned.). See Faletto and Bien (2022) for
#' details. Default is "sparse".
#' @param cutoff Numeric; getCssSelections will select and return only of those
#' clusters with selection proportions equal to at least cutoff. Must be between
#' 0 and 1. Default is 0 (in which case either all clusters are selected, or
#' max_num_clusts are selected, if max_num_clusts is specified).
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.) Default is 1.
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' max_num_clusts is ignored).
#' @return A named list with two items. \item{selected_clusts}{A named list of
#' integer vectors; each vector contains the indices of the features in one of
#' the selected clusters.} \item{selected_feats}{A named integer vector; the
#' indices of the features with nonzero weights from all of the selected
#' clusters.}
#' @author Gregory Faletto, Jacob Bien
#' @references Faletto, G., & Bien, J. (2022). Cluster Stability Selection.
#' \emph{arXiv preprint arXiv:2201.00494}.
#' \url{https://arxiv.org/abs/2201.00494}.
#' @export
getCssSelections <- function(css_results, weighting="sparse", cutoff=0,
    min_num_clusts=1, max_num_clusts=NA){
    # Check inputs
    stopifnot(class(css_results) == "cssr")
    checkCutoff(cutoff)
    checkWeighting(weighting)

    p <- ncol(css_results$feat_sel_mat)

    checkMinNumClusts(min_num_clusts, p, length(css_results$clusters))

    max_num_clusts <- checkMaxNumClusts(max_num_clusts, min_num_clusts, p,
        length(css_results$clusters))

    sel_results <- getSelectedClusters(css_results, weighting, cutoff,
        min_num_clusts, max_num_clusts)

    # sel_results$selected_clusts is guaranteed to have length at least 1 by
    # getSelectedClusters
    sel_clust_names <- names(sel_results$selected_clusts)

    stopifnot(length(sel_clust_names) >= 1)
    stopifnot(all(sel_clust_names %in% names(css_results$clusters)))

    sel_clusts <- list()
    for(i in 1:length(sel_clust_names)){
        sel_clusts[[i]] <- css_results$clusters[[sel_clust_names[i]]]
        names(sel_clusts)[i] <- sel_clust_names[i]
    }

    stopifnot(is.list(sel_clusts))
    stopifnot(length(sel_clusts) == length(sel_clust_names))

    # sel_results$selected_feats is guaranteed to have length at least as long
    # as sel_results$selected_clusts by getSelectedClusters
    return(list(selected_clusts=sel_clusts,
        selected_feats=sel_results$selected_feats))
}



### BELOW IS DONE AND IN RMD FILE

#' Prepares a data.frame summarazing cluster stability selection output to print
#'
#' Print a summary of the information from the css function.
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param cutoff Numeric; the outputted data.frame will display only those
#' clusters with selection proportions equal to at least cutoff. Must be between
#' 0 and 1. Default is 0 (in which case either all clusters are displayed, or
#' max_num_clusts are, if max_num_clusts is specified).
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.) Default is 1.
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' max_num_clusts is ignored).
#' @return A data.frame; each row contains a cluster, arranged in decreasing
#' order of cluster selection proportion from top to bottom. The columns are
#' ClustName (the name of the cluster that was either provided to css or made by
#' css if no name was provided); ClustProtoName (the name of the selection
#' prototype from the cluster, which is the feature with the greatest individual
#' selection proportion among all the cluster members, with ties broken by
#' choosing the feature with the highest correlation with the response if the
#' response is real-valued; only returned if the features are named),
#' ClustProtoNum (the column number of the prototype in the X matrix provided to
#' css), and ClustSize (the size of the cluster).
#' @author Gregory Faletto, Jacob Bien
#' @export
printCssDf <- function(css_results, cutoff=0, min_num_clusts=1,
    max_num_clusts=NA){
    # Check inputs
    stopifnot(class(css_results) == "cssr")
    checkCutoff(cutoff)

    p <- ncol(css_results$feat_sel_mat)

    checkMinNumClusts(min_num_clusts, p, length(css_results$clusters))

    max_num_clusts <- checkMaxNumClusts(max_num_clusts, min_num_clusts, p,
        length(css_results$clusters))

    sel_clusts <- getCssSelections(css_results, cutoff=cutoff,
        min_num_clusts=min_num_clusts,
        max_num_clusts=max_num_clusts)$selected_clusts

    # sel_clusts is guaranteed to have length at least 1 by
    # getCssSelections 

    # Get prototypes (feature from each cluster with highest selection
    # proportion, breaking ties by using marginal correlations of features with
    # y from data provided to css if y is real-valued)
    prototypes <- getSelectionPrototypes(css_results, sel_clusts)
    
    # Cluster selection proportions
    if(length(sel_clusts) > 1){
        sel_clust_sel_props <- colMeans(css_results$clus_sel_mat[,
            names(sel_clusts)])
    } else{
        sel_clust_sel_props <- mean(css_results$clus_sel_mat[,
            names(sel_clusts)])
    }
    

    # Data.frame: name of cluster, cluster prototype, selection proportion,
    # cluster size

    if(!is.null(names(prototypes))){
        print_df <- data.frame(ClustName=names(sel_clusts),
            ClustProtoName=names(prototypes), ClustProtoNum=unname(prototypes),
            ClustSelProp=sel_clust_sel_props, ClustSize=lengths(sel_clusts))
    } else{
        print_df <- data.frame(ClustName=names(sel_clusts),
            ClustProtoNum=unname(prototypes), ClustSelProp=sel_clust_sel_props,
            ClustSize=lengths(sel_clusts))
    }

    print_df <- print_df[order(print_df$ClustSelProp, decreasing=TRUE), ]

    rownames(print_df) <- NULL

    stopifnot(is.data.frame(print_df))
    stopifnot(nrow(print_df) >= 1)

    return(print_df)
}

### BELOW IS DONE AND IN RMD FILE

#' Print cluster stability selection output
#'
#' Print a summary of the information from the css function (using output from
#' printCssDf function).
#' @param x An object of class "cssr" (the output of the function css).
#' @param cutoff Numeric; print.cssr will display only those
#' clusters with selection proportions equal to at least cutoff. Must be between
#' 0 and 1. Default is 0 (in which case either all clusters are displayed, or
#' max_num_clusts are, if max_num_clusts is specified).
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.) Default is 1.
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' max_num_clusts is ignored).
#' @param ... Additional arguments to generic print.data.frame function
#' @return A data.frame; each row contains a cluster, arranged in decreasing
#' order of cluster selection proportion from top to bottom. The columns are
#' ClustName (the name of the cluster that was either provided to css or made by
#' css if no name was provided); ClustProtoName (the name of the selection
#' prototype from the cluster, which is the feature with the greatest individual
#' selection proportion among all the cluster members, with ties broken by
#' choosing the feature with the highest correlation with the response if the
#' response is real-valued; only returned if the features are named),
#' ClustProtoNum (the column number of the prototype in the X matrix provided to
#' css), and ClustSize (the size of the cluster).
#' @author Gregory Faletto, Jacob Bien
#' @export
print.cssr <- function(x, cutoff=0, min_num_clusts=1, max_num_clusts=NA, ...){
    df <- printCssDf(css_results=x, cutoff=cutoff,
        min_num_clusts=min_num_clusts, max_num_clusts=max_num_clusts)
    print.data.frame(df, ...)
}

### BELOW IS DONE AND IN RMD FILE

#' Obtain a design matrix of cluster representatives
#'
#' Takes a matrix of observations from the original feature space and returns
#' a matrix of representatives from the selected clusters based on the results
#' of cluster stability selection.
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param newX A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to generate the design matrix of cluster
#' representatives. Must contain the same features (in the same
#' number of columns) as the X matrix provided to css, and if the columns of
#' newX are labeled, the names must match the variable names provided to css.
#' newX may be omitted if train_inds were provided to css to set aside
#' observations for model estimation. If this is the case, then when newX is
#' omitted getCssDesign will return a design matrix of cluster representatives
#' formed from the train_inds observations from the matrix X provided to css.
#' (If no train_inds were provided to css, newX must be provided to
#' getCssDesign.) Default is NA.
#' @param weighting Character; determines how to calculate the weights to
#' combine features from the selected clusters into weighted averages, called
#' cluster representatives. Must be one of "sparse", "weighted_avg", or
#' "simple_avg'. For "sparse", all the weight is put on the most frequently
#' selected individual cluster member (or divided equally among all the clusters
#' that are tied for the top selection proportion if there is a tie). For
#' "weighted_avg", the weight used for each cluster member is calculated in
#' proportion to the individual selection proportions of each feature. For
#' "simple_avg", each cluster member gets equal weight regardless of the
#' individual feature selection proportions (that is, the cluster representative
#' is just a simple average of all the cluster members). See Faletto and Bien
#' (2022) for details. Default is "weighted_avg".
#' @param cutoff Numeric; getCssDesign will only include those clusters with
#' selection proportions equal to at least cutoff. Must be between 0 and 1.
#' Default is 0 (in which case either all clusters are used, or max_num_clusts
#' are used, if max_num_clusts is specified).
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.) Default is 1.
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' max_num_clusts is ignored).
#' @return A design matrix with either nrow(newX) (or length(train_inds), if
#' train_inds was provided to css and newX was not provided to getCssDesign)
#' observations and number of columns equal to the number of selected clusters,
#' containing the cluster representatives for each cluster.
#' @author Gregory Faletto, Jacob Bien
#' @export
getCssDesign <- function(css_results, newX=NA, weighting="weighted_avg",
    cutoff=0, min_num_clusts=1, max_num_clusts=NA){
    # Check inputs
    stopifnot(class(css_results) == "cssr")

    check_results <- checkNewXProvided(newX, css_results)

    newX <- check_results$newX
    newXProvided <- check_results$newXProvided

    rm(check_results)

    n_train <- nrow(newX)

    results <- checkXInputResults(newX, css_results$X)

    newX <- results$newx
    feat_names <- results$feat_names

    rm(results)

    n <- nrow(newX)
    p <- ncol(newX)

    checkCutoff(cutoff)
    checkWeighting(weighting)
    checkMinNumClusts(min_num_clusts, p, length(css_results$clusters))

    max_num_clusts <- checkMaxNumClusts(max_num_clusts, min_num_clusts, p,
        length(css_results$clusters))

    # Take provided training design matrix and testX and turn them into
    # matrices of cluster representatives using information from css_results
    if(newXProvided){
        newX_clusters <- formCssDesign(css_results, weighting, cutoff,
            min_num_clusts, max_num_clusts, newx=newX)
    } else{
        newX_clusters <- formCssDesign(css_results, weighting, cutoff,
            min_num_clusts, max_num_clusts)
    }

    return(newX_clusters)
}

### BELOW IS DONE AND IN RMD FILE

#' Obtain a selected set of clusters and features using cluster stability
#' selection
#'
#' Takes in data X and y and returns a set of clusters (and a set of features)
#' that are useful for predicting y from the data in X. This is a wrapper
#' function for css and getCssSelections. Using cssSelect is simpler, but it
#' has fewer options, and it executes the full (computationally expensive)
#' subsampling procedured every time it is called. In contrast, css can be
#' called just once, and then getCssSelections can quickly return results using
#' different values of cutoff, max_num_clusts, etc. from the calculations done
#' in one call to css.
#' @param X An n x p numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the p >= 2 features/predictors.
#' @param y A length-n numeric vector containing the responses; `y[i]` is the
#' response corresponding to observation `X[i, ]`. (Note that for the css
#' function, y does not have to be a numeric response, but for this function,
#' the underlying selection procedure is the lasso, so y must be a real-valued
#' response.)
#' @param clusters Optional; either an integer vector of a list of integer
#' vectors; each vector should contain the indices of a cluster of features (a
#' subset of 1:p). (If there is only one cluster, clusters can either be a list
#' of length 1 or an integer vector.) All of the provided clusters must be
#' non-overlapping. Every feature not appearing in any cluster will be assumed
#' to be unclustered (that is, they  will be treated as if they are in a
#' "cluster" containing only themselves). If clusters is a list of length 0 (or
#' a list only containing clusters of length 1), then css() returns the same
#' results as stability selection (so feat_sel_mat will be identical to
#' clus_sel_mat). Names for the clusters will be needed later; any clusters that
#' are not given names in the list clusters will be given names automatically by
#' css. Default is list() (so no clusters are specified, and every feature is
#' assumed to be in a "cluster" containng only itself).
#' @param lambda Optional; the tuning parameter to be used by the lasso for
#' feature selection in each subsample. If lambda is not provided, cssSelect
#' will choose one automatically by cross-validation. Default is NA.
#' @param cutoff Numeric; cssSelect will only select those clusters with
#' selection proportions equal to at least cutoff. Must be between 0 and 1.
#' Default is NA (in which case max_num_clusts are used).
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' either cutoff is used to choose the number of clusters, or if cutoff was also
#' unspecified, cssSelect chooses max_num_clusts by cross-validation).
#' @param auto_select_size Logical; if TRUE, then max_num_clusts will be
#' automatically estimated using the lasso with cross-validation. Default is
#' TRUE, though his argument is ignored if either cutoff or max_num_clusts is
#' provided. (If desired output is to return all clusters, you should set
#' auto_select_size to FALSE and do not provide cutoff or max_num_clusts.)
#' @return A named list with two items. \item{selected_clusts}{A list of
#' integer vectors; each vector contains the indices of one of the selected
#' clusters.} \item{selected_feats}{An integer vector; the indices of the
#' all of the selected features within all of the selected clusters (typically
#' only one feature is selected from each cluster).}
#' @author Gregory Faletto, Jacob Bien
#' @export
cssSelect <- function(X, y, clusters = list(), lambda=NA, cutoff=NA,
    max_num_clusts=NA, auto_select_size=TRUE){

    # Check inputs (most inputs will be checked by called functions)

    stopifnot(!is.na(auto_select_size))
    stopifnot(length(auto_select_size) == 1)
    stopifnot(is.logical(auto_select_size))

    stopifnot(is.matrix(X) | is.data.frame(X))
    stopifnot(all(!is.na(X)))


    # Check if x is a matrix; if it's a data.frame, convert to matrix.
    if(is.data.frame(X)){
        p <- ncol(X)
        X <- stats::model.matrix(~ ., X)
        X <- X[, colnames(X) != "(Intercept)"]

        if(length(clusters) > 0 & (p != ncol(X))){
            stop("When stats::model.matrix converted the provided data.frame X to a matrix, the number of columns changed (probably because the provided data.frame contained a factor variable with at least three levels). Please convert the data.frame X to a matrix yourself using model.matrix and provide cluster assignments according to the columns of the new matrix.")
        }

    }

    stopifnot(is.matrix(X))
    stopifnot(all(!is.na(X)))

    if(!is.numeric(y) & !is.integer(y)){
        stop("The provided y must be real-valued, because cssSelect uses the lasso for feature selection. (In order to use a different form of response, use the css function and provide your own selection function accommodating your choice of y.)")
    }

    stopifnot(length(lambda) == 1)
    if(is.na(lambda)){
        lambda <- getLassoLambda(X, y)
    }

    css_results <- css(X, y, lambda, clusters)

    # If no indication of how to select model size was provided, choose model
    # size by cross-validation
    if(is.na(cutoff) & is.na(max_num_clusts)){
        if(auto_select_size){
            max_num_clusts <- getModelSize(X, y, css_results$clusters)
        }
    }

    if(is.na(cutoff)){
        cutoff <- 0
    }

    # Get selected features
    getCssSelections(css_results, weighting="sparse", cutoff=cutoff,
        min_num_clusts=1, max_num_clusts=max_num_clusts)
}

### BELOW IS DONE AND IN RMD FILE

#' Wrapper function to generate predictions from cluster stability selected
#' model in one step
#'
#' Select clusters using cluster stability selection, form cluster
#' representatives, fit a linear model, and generate predictions from a matrix
#' of unlabeled data. This is a wrapper function for css and getCssPreds. Using
#' cssPredict is simpler, but it has fewer options, and it executes the full
#' (computationally expensive) subsampling procedured every time it is called.
#' In contrast, css can be called just once, and then cssPredict can quickly
#' return results for different matrices of new data or using different values
#' of cutoff, max_num_clusts, etc. by using the calculations done in one call to
#' css.
#'
#' @param X_train_selec An n x p numeric matrix (preferably) or a data.frame
#' (which will be coerced internally to a matrix by the function model.matrix)
#' containing the p >= 2 features/predictors. The data from X_train_selec and
#' y_train_selec will be split into two parts; half of the data will be used for
#' feature selection by cluster stability selection, and half will be used for
#' estimating a linear model on the selected cluster representatives.
#' @param y_train_selec A length-n numeric vector containing the responses;
#' `y[i]` is the response corresponding to observation `X[i, ]`. Unlke the more
#' general setup of css, y_train_selec must be real-valued because predictions
#' will be generated by ordinary least squares.
#' @param X_test A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to generate predictions. Must contain the same
#' features (in the same number of columns) as X_train_selec, and if the columns
#' of X_test are named, they must match the names of X_train_selec.
#' @param clusters Optional; either an integer vector of a list of integer
#' vectors; each vector should contain the indices of a cluster of features (a
#' subset of 1:p). (If there is only one cluster, clusters can either be a list
#' of length 1 or an integer vector.) All of the provided clusters must be
#' non-overlapping. Every feature not appearing in any cluster will be assumed
#' to be unclustered (that is, they  will be treated as if they are in a
#' "cluster" containing only themselves). If clusters is a list of length 0 (or
#' a list only containing clusters of length 1), then css() returns the same
#' results as stability selection (so feat_sel_mat will be identical to
#' clus_sel_mat). Names for the clusters will be needed later; any clusters that
#' are not given names in the list clusters will be given names automatically by
#' css. Default is list() (so no clusters are specified, and every feature is
#' assumed to be in a "cluster" containng only itself).
#' @param lambda Optional; the tuning parameter to be used by the lasso for
#' feature selection in each subsample. If lambda is not provided, cssPredict
#' will choose one automatically by cross-validation. Default is NA.
#' @param cutoff Numeric; getCssPreds will make use only of those clusters with
#' selection proportions equal to at least cutoff. Must be between 0 and 1.
#' Default is 0 (in which case either all clusters are used, or max_num_clusts
#' are used, if max_num_clusts is specified).
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' max_num_clusts is ignored).
#' @param train_inds Optional; an integer or numeric vector containing the
#' indices of observations in X and y to set aside for model training after
#' feature selection. If train_inds is not provided, half of the data will be
#' used for feature selection and half for model estimation (chosen at random).
#' @param auto_select_size Logical; if TRUE, then max_num_clusts will be
#' automatically estimated using the lasso with cross-validation. Default is
#' TRUE, though his argument is ignored if either cutoff or max_num_clusts is
#' provided. (If desired output is to generate predictions using all clusters,
#' you should set auto_select_size to FALSE and do not provide cutoff or
#' max_num_clusts.)
#' @return A numeric vector of length nrow(X_test) of predictions
#' corresponding to the observations from X_test.
#' @author Gregory Faletto, Jacob Bien
#' @export
cssPredict <- function(X_train_selec, y_train_selec, X_test, clusters=list(),
    lambda=NA, cutoff=NA, max_num_clusts=NA, train_inds=NA,
    auto_select_size=TRUE){

    # Check inputs (most inputs will be checked by called functions)
    if(!is.numeric(y_train_selec) & !is.integer(y_train_selec)){
        stop("The provided y_train_selec must be real-valued, because predictions will be generated by ordinary least squares regression.")
    }

    stopifnot(!is.na(auto_select_size))
    stopifnot(length(auto_select_size) == 1)
    stopifnot(is.logical(auto_select_size))

    stopifnot(is.matrix(X_train_selec) | is.data.frame(X_train_selec))
    stopifnot(all(!is.na(X_train_selec)))

    # Check if x is a matrix; if it's a data.frame, convert to matrix.
    if(is.data.frame(X_train_selec)){
        p <- ncol(X_train_selec)
        X_train_selec <- stats::model.matrix(~ ., X_train_selec)
        X_train_selec <- X_train_selec[, colnames(X_train_selec) !=
            "(Intercept)"]

        if(length(clusters) > 0 & (p != ncol(X_train_selec))){
            stop("When stats::model.matrix converted the provided data.frame X_train_selec to a matrix, the number of columns changed (probably because the provided data.frame contained a factor variable with at least three levels). Please convert X_train_selec to a matrix yourself using model.matrix and provide cluster assignments according to the columns of the new matrix.")
        }
    }

    stopifnot(is.matrix(X_train_selec))
    stopifnot(all(!is.na(X_train_selec)))

    n <- nrow(X_train_selec)

    if(any(is.na(train_inds))){
        train_inds <- sample(n, size=round(n/2))
    }

    stopifnot(length(lambda) == 1)
    if(is.na(lambda)){
        lambda <- getLassoLambda(X_train_selec[setdiff(1:n, train_inds), ],
            y_train_selec[setdiff(1:n, train_inds)]) 
    }

    css_results <- css(X=X_train_selec, y=y_train_selec, lambda=lambda,
        clusters=clusters, train_inds=train_inds)

    # If no indication of how to select model size was provided, choose model
    # size by cross-validation
    if(is.na(cutoff) & is.na(max_num_clusts)){
        if(auto_select_size){
            max_num_clusts <- getModelSize(X_train_selec[train_inds, ],
                y_train_selec[train_inds], css_results$clusters)
        }
    }

    if(is.na(cutoff)){
        cutoff <- 0
    }

    # Get predictions
    getCssPreds(css_results, testX=X_test, weighting="weighted_avg",
        cutoff=cutoff, max_num_clusts=max_num_clusts)
}

#
#
#
#
#
#
# helper functions
#
#
#
#
#
#
#

### BELOW IS DONE AND IN RMD FILE

#' Creates lists of subsamples for stability selection.
#'
#' @param n Integer or numeric; sample size of the data set.
#' @param p Integer or numeric; number of features.
#' @param B Integer or numeric; the number of subsamples. Note: For
#' sampling.type=="MB" the number of lasso fits will be `B`; for
#' sampling_type="SS" the number of lasso fits will be `2*B`.
#' @param sampling_type A character vector (either "SS" or "MB"); the sampling
#' type used for stability selection.
#' @param num_feats_remove Integer; number of features select automatically on
#' every iteration. Determined earlier from input prop_feats_remove to css.
#' @return A list of length `B` (or `2*B` for sampling_type = "SS"). If
#' prop_feats_remove = 0, each list element is an integer vector of length
#' floor(`n/2`) containing the indices of a subsample of 1:`n`. (For
#' sampling_type=="SS", the last `B` subsamples will be complementary pairs of
#' the first `B` subsamples; see Faletto and Bien 2022 or Shah and Samworth 2013
#' for details.) If prop_feats_remove > 0, each element is a named list with
#' members "subsample" (same as above) and "feats_to_keep", a logical vector
#' of length p; feats_to_keep[j] = TRUE if feature j is chosen for this
#' subsample, and false otherwise.
#' @author Gregory Faletto, Jacob Bien
#' @references Faletto, G., & Bien, J. (2022). Cluster Stability Selection.
#' \emph{arXiv preprint arXiv:2201.00494}.
#' \url{https://arxiv.org/abs/2201.00494}. \cr Shah, R. D., & Samworth, R. J.
#' (2013). Variable selection with error control: Another look at stability
#' selection. \emph{Journal of the Royal Statistical Society. Series B:
#' Statistical Methodology}, 75(1), 55–80.
#' \url{https://doi.org/10.1109/RITA.2014.2302071}.
createSubsamples <- function(n, p, B, sampling_type, prop_feats_remove=0){

    # Check inputs

    stopifnot(length(n) == 1)
    stopifnot(is.numeric(n) | is.integer(n))
    stopifnot(n == round(n))
    stopifnot(n > 0)

    stopifnot(length(p) == 1)
    stopifnot(is.numeric(p) | is.integer(p))
    stopifnot(p == round(p))
    stopifnot(p > 0)

    checkSamplingType(sampling_type)
    checkPropFeatsRemove(prop_feats_remove, p)

    if(prop_feats_remove == 0){
        subsamples <- getSubsamps(n, B, sampling_type)
        return(subsamples)
    } else{
        # In this case, we generate subsamples as well as logical vectors
        # of features to keep
        subsamps_and_feats <- list()
        subsamples <- getSubsamps(n, B, sampling_type)
        for(i in 1:B){
            # Logical p-vector, where each entry is TRUE with probability
            # 1 - prop_feats_remove
            feats_to_keep_i <- as.logical(stats::rbinom(n=p, size=1,
                prob=1 - prop_feats_remove))
            # Make sure at least two entries are equal to TRUE (so that at
            # least two features are present for every subsample)--if not,
            # randomly choose features to add
            while(sum(feats_to_keep_i) < 2){
                false_inds <- which(!feats_to_keep_i)
                sel_feat <- sample(false_inds, size=1)
                feats_to_keep_i[sel_feat] <- TRUE
            }
            subsamps_and_feats[[i]] <- list(subsample=subsamples[[i]],
                feats_to_keep=feats_to_keep_i)
        }

        if(sampling_type=="SS"){
            stopifnot(length(subsamples) == 2*B)
            for(i in 1:B){
                # Keep the same features as in the other subsample (this
                # ensures that the theoretical guarantee of Shah and Samworth
                # 2013 remains valid on every individual pair of subsamples)
                subsamps_and_feats[[B + i]] <- list(subsample=subsamples[[B + i]],
                    feats_to_keep=subsamps_and_feats[[i]]$feats_to_keep)
            }
        }

        # Check output
        stopifnot(all(names(subsamps_and_feats) == c("subsample",
            "feats_to_keep")))

        return(subsamps_and_feats)
    }
    # Shouldn't be possible to reach this part of function
    stop("createSubsamples failed to return anything")
}

### BELOW IS DONE AND IN RMD FILE

#' Generate list of subsamples
#'
#` Generate list of `B` (or `2*B` for sampling_type="SS") subsamples of size
#` `n/2`
#' @param n Integer or numeric; sample size of the data set.
#' @param B Integer or numeric; the number of subsamples. Note: For
#' sampling.type=="MB" the number of lasso fits will be `B`; for
#' sampling_type="SS" the number of lasso fits will be `2*B`.
#' @param sampling_type A character vector (either "SS" or "MB"); the sampling
#' type used for stability selection.
#' @return A list of length `B` (or `2*B` for sampling_type="SS"), where each
#' element is an integer vector of length floor(`n/2`) containing the indices
#' of a subsample of 1:`n`. For sampling_type=="SS", the last `B` subsamples
#' will be complementary pairs of the first `B` subsamples (see Faletto and
#' Bien 2022 or Shah and Samworth 2013 for details).
#' @author Gregory Faletto, Jacob Bien
#' @references Faletto, G., & Bien, J. (2022). Cluster Stability Selection.
#' \emph{arXiv preprint arXiv:2201.00494}.
#' \url{https://arxiv.org/abs/2201.00494}. \cr Shah, R. D., & Samworth, R. J.
#' (2013). Variable selection with error control: Another look at stability
#' selection. \emph{Journal of the Royal Statistical Society. Series B:
#' Statistical Methodology}, 75(1), 55–80.
#' \url{https://doi.org/10.1109/RITA.2014.2302071}.
getSubsamps <- function(n, B, sampling_type){
    subsamples <- list()
    for(i in 1:B){
        subsamples[[i]] <- sort(sample.int(n=n, size=floor(n/2), replace=FALSE))
    }
    stopifnot(length(subsamples) == B)
    # TODO(gregfaletto): add support for sampling_type="MS"
    if(sampling_type=="SS"){
        for(i in 1:B){
            # For the ith entry, take a subsample of size floor(n/2) from the
            # remaining n - floor(n/2) observations. (Only necessary to actually
            # take the subsample if n is odd; if not, the set we want is
            # setdiff(1:n, subsamples[[i]]), so skip the sample function.)
            if(n/2 == floor(n/2)){
                subsamples[[B + i]] <- sort(setdiff(1:n, subsamples[[i]]))
            } else{
                subsamples[[B + i]] <- sort(sample(x=setdiff(1:n,
                    subsamples[[i]]), size=floor(n/2), replace=FALSE))
            }

            # Check output

            stopifnot(is.integer(subsamples[[B + i]]))
            stopifnot(all(subsamples[[B + i]] ==
                round(subsamples[[B + i]])))
            stopifnot(floor(n/2) == length(subsamples[[B + i]]))
            stopifnot(length(subsamples[[B + i]]) ==
                length(unique(subsamples[[B + i]])))
        }
        stopifnot(length(subsamples) == 2*B)
    }
    return(subsamples)
}

### BELOW IS DONE AND IN RMD FILE

#' Generates matrix of selection indicators from stability selection.
#'
#' @param x an n x p numeric matrix or a data.frame containing the predictors.
#' @param y A response vector; can be any response that takes the form of a
#' length n vector and is used (or not used) by fitfun. Typically (and for
#' default fitfun = cssLasso), y should be an n-dimensional numeric vector
#' containing the response.
#' @param lambda A tuning parameter or set of tuning parameters that may be used
#' by the feature selection method. For example, in the default case when
#' fitfun = cssLasso, lambda is a numeric: the penalty to use for each lasso
#' fit.
#' @param B Integer or numeric; the number of subsamples. Note: For
#' sampling.type=="MB" the number of lasso fits will be `B`; for
#' sampling_type="SS" the number of lasso fits will be `2*B`.
#' @param sampling_type A character vector (either "SS" or "MB"); the sampling
#' type used for stability selection.
#' @param subsamps_object A list of length `B` (or `2*B` if sampling_type="SS"),
#' where each element is one of the following: \item{subsample}{An integer
#' vector of size `n/2` containing the indices of the observations in the
#' subsample.} \item{drop_var_input}{A named list containing two elements: one
#' named "subsample", matching the previous description, and a logical vector
#' named "feats_to_keep" containing the indices of the features to be
#' automatically selected.} (The first object is the output of the function
#' createSubsamples when the provided prop_feats_remove is 0, the default, and
#' the second object is the output of createSubsamples when prop_feats_remove >
#' 0.)
#' @param num_cores Optional; an integer. If using parallel processing, the
#' number of cores to use for parallel processing (num_cores will be supplied
#' internally as the mc.cores argument of parallel::mclapply).
#' @param fitfun A function that takes in arguments X, y, and lambda and returns
#' a vector of indices of the columns of X (selected features). Default is
#' cssLasso, an implementation of the lasso using glmnet.
#' @return A binary integer matrix of dimension `B` x `p` (if sampling_type ==
#' "MB") or `2*B` x `p` (if sampling_type == "SS"). res[i, j] = 1 if feature j
#' was selected on subsample i and equals 0 otherwise. (That is, each row is a
#' selected set.)
#' @author Gregory Faletto, Jacob Bien
getSelMatrix <- function(x, y, lambda, B, sampling_type, subsamps_object,
    num_cores, fitfun=cssLasso){

    # Check inputs

    stopifnot(is.matrix(x))
    stopifnot(all(!is.na(x)))

    n <- nrow(x)
    p <- ncol(x)

    stopifnot(length(y) == n)
    stopifnot(!is.matrix(y))
    # Intentionally don't check y or lambda further to allow for flexibility--these
    # inputs should be checked within fitfun.

    checkSamplingType(sampling_type)

    # Get list of selected feature sets from subsamples

    res_list <- parallel::mclapply(X=subsamps_object, FUN=cssLoop, x=x, y=y,
        lambda=lambda, fitfun=fitfun, mc.cores=num_cores)

    # Store selected sets in B x p (or `2*B` x p for "SS") binary matrix
    if(sampling_type=="SS"){
        res <- matrix(0L, 2*B, p)
    } else if(sampling_type=="MB"){
        res <- matrix(0L, B, p)
    } else{
        stop("!(sampling_type %in% c(SS, MB)) (don't know how to set dimensions of res")
    }

    stopifnot(length(res_list) == nrow(res))

    # Get selected sets into matrix res

    for(i in 1:nrow(res)){
        if(length(res_list[[i]]) == 0){
            # If no features are selected, don't fill in anything in this row
            next
        }

        if(!is.integer(res_list[[i]])){
            print(paste("failed on iteration", i))
            print(res_list[[i]])
            stop("Something seems to be wrong with the feature selection method (fitfun failed to return an integer vector)")
        }
        stopifnot(length(res_list[[i]]) <= p & length(res_list[[i]]) > 0)
        stopifnot(all(!is.na(res_list[[i]])))
        stopifnot(max(res_list[[i]]) <= p)
        stopifnot(min(res_list[[i]]) >= 1)
        stopifnot(length(res_list[[i]]) == length(unique(res_list[[i]])))
        stopifnot(!("try-error" %in% class(res_list[[i]]) |
            "error" %in% class(res_list[[i]]) |
            "simpleError" %in% class(res_list[[i]]) |
            "condition" %in% class(res_list[[i]])))

        # Store selected variables in the ith row of res
        res[i, res_list[[i]]] <- 1L
    }

    # Check output

    stopifnot(is.matrix(res))
    if(sampling_type=="SS"){
        stopifnot(nrow(res) == 2*B)
    } else{
        stopifnot(nrow(res) == B)
    }
    stopifnot(ncol(res) == p)
    stopifnot(all(res %in% c(0, 1)))

    return(res)
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function run on each subsample
#' 
#' Runs provided feature selection method fitfun on each subsample for cluster
#' stability selection (this function is called within mclapply).
#' @param input Could be one of two things: \item{subsample}{An integer vector
#' of size `n/2` containing the indices of the observations in the subsample.}
#' \item{drop_var_input}{A named list containing two elements: one named
#' "subsample" and the same as the previous description, and a logical vector
#' named "feats_to_keep" containing the indices of the features to be
#' automatically selected.} (The first object is the output of the function
#' createSubsamples when the provided prop_feats_remove is 0, the default, and
#' the second object is the output of createSubsamples when prop_feats_remove >
#' 0.)
#' @param x an n x p numeric matrix containing the predictors. (This should be
#' the full design matrix provided to css.)
#' @param y A response; can be any response that takes the form of a length n
#' vector and is used (or not used) by fitfun. Typically (and for default fitfun
#' = cssLasso), y should be an n-dimensional numeric vector containing the
#' response. This should be the full response provided to css.
#' @param lambda A tuning parameter or set of tuning parameters that may be used
#' by the feature selection method. For example, in the default case when
#' fitfun = cssLasso, lambda is a numeric: the penalty to use for each lasso
#' fit.
#' @param fitfun A function that takes in arguments X, y, and lambda and returns
#' a vector of indices of the columns of X (selected features).
#' @return An integer vector; the indices of the features selected by fitfun.
#' @author Gregory Faletto, Jacob Bien
cssLoop <- function(input, x, y, lambda, fitfun){
    # Check inputs
    stopifnot(is.matrix(x))
    stopifnot(all(!is.na(x)))

    colnames(x) <- character()
    n <- nrow(x)
    p <- ncol(x)

    stopifnot(length(y) == n)
    stopifnot(!is.matrix(y))
    # Intentionally don't check y or lambda further to allow for flexibility--these
    # inputs should be checked within fitfun.

    if(!is.list(input)){
        subsample <- input
        feats_to_keep <- rep(TRUE, p)
    } else{
        stopifnot(all(names(input) == c("subsample", "feats_to_keep")))
        subsample <- input$subsample
        feats_to_keep <- input$feats_to_keep
    }

    stopifnot(is.integer(subsample))
    stopifnot(all(subsample == round(subsample)))
    stopifnot(floor(n/2) == length(subsample))
    stopifnot(length(subsample) == length(unique(subsample)))

    stopifnot(is.logical(feats_to_keep))
    stopifnot(length(feats_to_keep) == p)

    selected <- do.call(fitfun, list(X=x[subsample, feats_to_keep],
        y=y[subsample], lambda=lambda))

    selected <- which(feats_to_keep)[selected]

    # Check output
    checkCssLoopOutput(selected, p, as.integer(which(feats_to_keep)))

    return(as.integer(selected))
}

### BELOW IS DONE AND IN RMD FILE

#' Provided fitfun implementing the lasso
#'
#' Function used to select features with the lasso on each subsample in cluster
#' stability selection. Uses glmnet implementation of the lasso.
#' @param X A design matrix containing the predictors. (In practice this will
#' be a subsample of the full design matrix provided to css.)
#' @param y A numeric vector containing the response.
#' @param lambda Numeric; a nonnegative number for the lasso penalty to use
#' on each subsample. (For now, only one lambda value can be provided to
#' cssLasso; in the future, we plan to allow for multiple lambda values to be
#' provided to cssLasso, as described in Faletto and Bien 2022.)
#' @return An integer vector; the indices of the features selected by the lasso.
#' @author Gregory Faletto, Jacob Bien
#' @references Faletto, G., & Bien, J. (2022). Cluster Stability Selection.
#' \emph{arXiv preprint arXiv:2201.00494}.
#' \url{https://arxiv.org/abs/2201.00494}. \cr Jerome Friedman, Trevor Hastie,
#' Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models
#' via Coordinate Descent. \emph{Journal of Statistical Software}, 33(1), 1-22.
#' URL \url{https://www.jstatsoft.org/v33/i01/}.
#' @export
cssLasso <- function(X, y, lambda){
    # Check inputs

    # TODO(gregfaletto) allow cssLasso to accept a vector of lambda values rather
    # than just a single one.
    checkCssLassoInputs(X, y, lambda)

    n <- nrow(X)
    p <- ncol(X)

    # Fit a lasso path (full path for speed, per glmnet documentation)

    lasso_model <- glmnet::glmnet(X, y, family="gaussian")
    stopifnot(all.equal(class(lasso_model), c("elnet", "glmnet")))

    # Get coefficients at desired lambda

    pred <- glmnet::predict.glmnet(lasso_model, type="nonzero",
        s=lambda, exact=TRUE, newx=X, x=X, y=y)

    if(is.null(pred[[1]])){return(integer())}

    stopifnot(is.data.frame(pred))
    stopifnot(!("try-error" %in% class(pred) | "error" %in% class(pred) |
        "simpleError" %in% class(pred) | "condition" %in% class(pred)))

    if(length(dim(pred)) == 2){
        selected_glmnet <- pred[, 1]
    } else if(length(dim(pred)) == 3){
        selected_glmnet <- pred[, 1, 1]
    } else if(length(dim(pred)) == 1){
        selected_glmnet <- pred
    } else{
        stop("length(dim(pred)) not 1, 2, or 3")
    }

    stopifnot(length(selected_glmnet) >= 1)
    stopifnot(length(selected_glmnet) <= ncol(X))
    stopifnot(all(selected_glmnet == round(selected_glmnet)))
    stopifnot(length(selected_glmnet) == length(unique(selected_glmnet)))
    selected_glmnet <- as.integer(selected_glmnet)

    selected <- sort(selected_glmnet)

    return(selected)
}

### BELOW IS DONE AND IN RMD FILE

#' Formats clusters in standardized way, optionally estimating cluster
#' prototypes
#'
#' @param clusters Either an integer vector of a list of integer vectors; each
#' vector should contain the indices of a cluster of features. (If there is only
#' one cluster, clusters can either be a list of length 1 or simply an integer
#' vector.) If clusters is specified then R is ignored.
#' @param p integer or numeric; the numbe of features in x (should match 
#' ncol(x), if x is provided)
#' @param clust_names A character vector of the names of the clusters in
#' clusters.
#' @param get_prototypes Logical: if TRUE, will identify prototype from each
#' cluster (the feature from each cluster that is most correlated with the
#' response) for the protolasso. In this case, x and y must be provided.
#' @param x n x p numeric matrix; design matrix. Only needs to be provided if
#' get_prototypes is TRUE.
#' @param y Numeric response vector; only needs to be provided if get_prototypes
#' is TRUE. Note: in general, the css function does not require y to be a
#' numeric vector, because the provided fitfun could use a different form of y
#' (for example, a categorical response variable). However, y must be numeric in
#' order to provide prototypes because the prototypes are determined using the
#' correlation between cluster members (columns of x) and y.
#' @param R Numeric p x p matrix; not currently used. Entry ij contains the 
#' "substitutive value" of feature i for feature j (diagonal must consist of
#' ones, all entries must be between 0 and 1, and matrix must be symmetric)
#' @return A named list with the following elements: \item{clusters}{A named
#' list where each entry is an integer vector of indices of features that are in
#' a common cluster. (The length of list clusters is equal to the number of
#' clusters.) All identified clusters are non-overlapping. All features appear
#' in exactly one cluster (any unclustered features will be put in their own
#' "cluster" of size 1).} \item{multiple}{Logical; TRUE if there is more than
#' one cluster of size greater than 1, FALSE otherwise.} \item{prototypes}{only
#' returned if get_prototypes=TRUE. An integer vector whose length is equal to
#' the number of clusters. Entry i is the index of the feature belonging to
#' cluster i that is most highly correlated with y (that is, the prototype for
#' the cluster, as in the protolasso; see Reid and Tibshirani 2016).}
#' @author Gregory Faletto, Jacob Bien
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}.
formatClusters <- function(clusters=NA, p=-1, clust_names=NA, 
    get_prototypes=FALSE, x=NA, y=NA, R=NA){

    # Check inputs
    clusters <- checkFormatClustersInput(clusters, p, clust_names,
        get_prototypes, x, y, R)

    n <- nrow(x)

    multiple <- FALSE

    if(any(lengths(clusters) > 1)){ # & length(clusters) > 1
        # Only care about clusters with more than one element (only ones that
        # need to be treated differently)
        # keep track of whether there's more than one cluster or not
        multiple <- sum(lengths(clusters) > 1) > 1
    }

    # For any features not already in a cluster, add a cluster containing only
    # that feature
    orig_length_clusters <- length(clusters)

    stopifnot(p >= 1)
    for(i in 1:p){
        feat_i_found <- FALSE
        if(orig_length_clusters > 0){
            for(j in 1:orig_length_clusters){
                # If i is in cluster j, break out of this loop and consider the
                # next i
                if(i %in% clusters[[j]]){
                    feat_i_found <- TRUE
                    break
                }
            }
        }

        # If feature i wasn't found in any cluster, add a cluster containing
        # just feature i
        if(!feat_i_found){
            clusters[[length(clusters) + 1]] <- i
        }
    }

    n_clusters <- length(clusters)

    # Add names to clusters
    if(is.null(names(clusters))){
        names(clusters) <- paste("c", 1:n_clusters, sep="")
    } else{
        # What clusters need names?
        unnamed_clusts <- which(is.na(names(clusters)) | names(clusters) == "")
        if(length(unnamed_clusts) > 0){
            proposed_clust_names <- paste("c", unnamed_clusts, sep="")
            # Replace any proposed cluster names that are already in use
            if(any(proposed_clust_names %in% names(clusters))){
                proposed_clust_names[proposed_clust_names %in% names(clusters)] <- paste("c",
                    unnamed_clusts[proposed_clust_names %in% names(clusters)] + p,
                    sep="")
            }
            while_counter <- 0
            while(any(proposed_clust_names %in% names(clusters))){
                proposed_clust_names[proposed_clust_names %in% names(clusters)] <- paste(proposed_clust_names[proposed_clust_names %in% names(clusters)],
                    "_1", sep="")
                while_counter <- while_counter + 1
                if(while_counter >= 100){
                    stop("Function formatClusters stuck in an infinite while loop")
                }
            }
            stopifnot(length(unnamed_clusts) == length(proposed_clust_names))
            names(clusters)[unnamed_clusts] <- proposed_clust_names
        }
    }

    # Check output

    checkClusters(clusters, p)
    stopifnot(is.logical(multiple))
    stopifnot(length(multiple) == 1)
    stopifnot(!is.na(multiple))

    if(get_prototypes){
        prototypes <- getPrototypes(clusters, x, y)

        return(list(clusters=clusters, multiple=multiple,
            prototypes=prototypes))
    } else{
        return(list(clusters=clusters, multiple=multiple))
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Estimate prototypes from a list of clusters
#'
#' Takes in list of clusters, x, and y and returns an integer vector (of length
#' equal to the number of clusters) of the indices of the feature prototypes
#' (the features from each cluster most correlated with the response).
#'
#' @param clusters A list where each entry is an integer vector of indices of
#' features that are in a common cluster. (The length of list clusters must be
#' equal to the number of clusters.) All identified clusters must be
#' non-overlapping. Must only include clusters of size 2 or larger.
#' @param x n x p numeric matrix; design matrix.
#' @param y Numeric response vector. Note: in general, the css function does not
#' require y to be a numeric vector, because the provided fitfun could use a
#' different form of y (for example, a categorical response variable). However,
#' y must be numeric in order to provide prototypes because the prototypes are
#' determined using the correlation between cluster members (columns of x) and
#' y.
#' @return An integer vector of the same length as clusters. Entry j is the
#' index of the feature identified as the prototype for cluster j.
#' @author Gregory Faletto, Jacob Bien
getPrototypes <- function(clusters, x, y){
    # Check inputs

    stopifnot(!is.list(clusters) | all(lengths(clusters) >= 1))
    stopifnot(is.list(clusters) | length(clusters) >= 1)

    stopifnot(all(!is.na(x)))
    stopifnot(is.matrix(x))

    n <- nrow(x)
    p <- ncol(x)

    checkY(y, n)

    # Identify prototypes
    if(length(clusters) > 0){
        if(is.list(clusters)){
            prototypes <- rep(as.integer(NA), length(clusters))
            for(i in 1:length(clusters)){
                prototypes[i] <- identifyPrototype(clusters[[i]], x, y)
            }
        } else{
            # If clusters is not a list, then there is only one cluster.
            prototypes <- identifyPrototype(clusters, x, y)
        }
    }  else{
        prototypes <- integer()
    }

    # Check output

    stopifnot(is.integer(prototypes))
    if(length(clusters) > 0){
        if(is.list(clusters)){
            stopifnot(length(prototypes) == length(clusters))
        } else {
            stopifnot(length(prototypes) == 1)
        }
    } else{
        stopifnot(length(prototypes) == 0)
    }

    stopifnot(all(!is.na(prototypes)))
    stopifnot(length(prototypes) == length(unique(prototypes)))

    return(prototypes)
}

### BELOW IS DONE AND IN RMD FILE

#' Estimate prototypes from a single cluster
#'
#' Takes in a single cluster, x, and y and returns an integer of the index of
#' the feature prototype (the feature from the cluster most correlated with the
#' response).
#'
#' @param cluster_members_i An integer vector of indices of features that are in
#' a common cluster. Must have length at least 2.
#' @param x n x p numeric matrix; design matrix.
#' @param y Numeric response vector. Note: in general, the css function does not
#' require y to be a numeric vector, because the provided fitfun could use a
#' different form of y (for example, a categorical response variable). However,
#' y must be numeric in order to provide prototypes because the prototypes are
#' determined using the correlation between cluster members (columns of x) and
#' y.
#' @return integer; the index of the feature identified as the prototype for
#' the cluster.
#' @author Gregory Faletto, Jacob Bien
identifyPrototype <- function(cluster_members_i, x, y){
    # Check input
    stopifnot(is.integer(cluster_members_i))
    # If cluster only has one member, that member is the prototype
    if(length(cluster_members_i) == 1){
        return(cluster_members_i)
    }

    # Choose which cluster member to represent cluster for stability
    # metric purposes by choosing the one most highly correlated
    # with y

    cors_i <- apply(x[, cluster_members_i], 2, corFunction, y=y)
    max_index_i <- which.max(cors_i)[1]

    stopifnot(length(max_index_i) == 1)
    stopifnot(max_index_i %in% 1:length(cluster_members_i))

    ret <- cluster_members_i[max_index_i]

    # Check output

    stopifnot(is.integer(ret))
    stopifnot(length(ret) == 1)
    stopifnot(ret %in% cluster_members_i)
    stopifnot(identical(x[, cluster_members_i][, max_index_i],
        x[, ret]))

    return(ret)
}

### BELOW IS DONE AND IN RMD FILE

#' Get cluster selection matrix
#'
#' Given a matrix of feature selection indicator variables and a list of
#' clusters of features, returns a matrix of cluster indicator variables.
#'
#' @param clusters A named list where each entry is an integer vector of indices
#' of features that are in a common cluster, as in the output of formatClusters.
#' (The length of list clusters is equal to the number of clusters.) All
#' identified clusters must be non-overlapping, and all features must appear in
#' exactly one cluster (any unclustered features should be in their own
#' "cluster" of size 1).
#' @param res A binary integer matrix. es[i, j] = 1 if feature j was selected
#' on subsample i and equals 0 otherwise, as in the output of getSelMatrix.
#' (That is, each row is a selected set.)
#' @return A binary integer matrix with the same number of rows
#' as res and length(clusters) columns. Entry i, j is 1 if at least
#' one member of cluster j was selected on subsample i, and 0 otherwise.
#' @author Gregory Faletto, Jacob Bien
getClusterSelMatrix <- function(clusters, res){

    # Check input
    checkGetClusterSelMatrixInput(clusters, res)

    p <- ncol(res)

    n_clusters <- length(clusters)

    # Matrix of cluster selection proportions (with n_clusters columns)
    res_n_clusters <- matrix(rep(0L, nrow(res)*n_clusters), nrow=nrow(res),
        ncol=n_clusters)
    colnames(res_n_clusters) <- names(clusters)

    for(j in 1:n_clusters){
        # Identify rows of res where at least one feature from this cluster
        # was selected
        rows_j_sel <- apply(res, 1, function(x){any(x[clusters[[j]]] == 1)})

        # Put ones in these rows of res_n_clusters[, j]
        res_n_clusters[rows_j_sel, j] <- 1L

        if(length(clusters[[j]]) <= 1){
            next
        }
    }

    # Check output
    stopifnot(is.matrix(res_n_clusters))
    stopifnot(identical(colnames(res_n_clusters), names(clusters)))
    stopifnot(all(res_n_clusters %in% c(0, 1)))
    stopifnot(ncol(res_n_clusters) == length(clusters))

    return(res_n_clusters)
}

### BELOW IS DONE AND IN RMD FILE

#' Generate covariance matrix for simulated clustered data
#'
#' @param p Integer or numeric; the total number of features in the covariance
#' matrix to be created, including latent features, the associated noisy proxies
#' with each latent feature, and other (weak signal and noise) features.
#' @param n_blocks Integer or numeric; the number of latent variables in the
#' data, each of which is associated with an observed cluster in X. Must be at
#' least 1.
#' @param block_size Integer or numeric; for each of the n_blocks latent
#' variables, the covariance matrix will include the original latent feature
#' plus block_size - 1 noisy proxies that are correlated with the latent
#' variable.
#' @param rho Integer or numeric; the covariance of the proxies in each cluster
#' with the latent variable (and each other). Note that the correlation between
#' the features in the cluster will be rho/var. rho cannot equal 0.
#' @param var Integer or numeric; the variance of all of the observed features
#' in X (both the proxies for the latent variables and the k_unclustered other
#' features). var cannot equal 0.
#' @return A `p` x `p` numeric matrix representing the covariance matrix for
#' the latent features, the associated proxies, and the remaining features. All
#' features not in a block will be independent from each other and the blocks
#' and have variance var.
#' @author Gregory Faletto, Jacob Bien
makeCovarianceMatrix <- function(p, nblocks, block_size, rho, var) {
    # Check inputs

    stopifnot(nblocks >= 1)
    stopifnot(rho != 0)
    stopifnot(var != 0)
    stopifnot(abs(rho) <= abs(var))
    stopifnot(block_size >= 2)
    stopifnot(p >= nblocks*block_size)

    # start with p x p identity matrix
    Sigma <- var*diag(p)

    # create matrix with nblocks rows, each containing a vector of
    # indices of highly correlated features
    block_feats <- matrix(seq(nblocks*block_size), nrow=nblocks, byrow=TRUE)

    stopifnot(length(unique(block_feats)) == length(block_feats))

    # add covariances of highly correlated features to sigma
    for(i in 1:nblocks){
        for(j in 1:(block_size - 1)){
            for(k in (j+1):block_size){
                feat_1 <- block_feats[i, j]
                feat_2 <- block_feats[i, k]
                Sigma[feat_1, feat_2] <- rho
                Sigma[feat_2, feat_1] <- rho
            }
        }
    }
    stopifnot(is.numeric(Sigma))
    stopifnot(is.matrix(Sigma))
    stopifnot(nrow(Sigma) == p & ncol(Sigma) == p)
    stopifnot(all(Sigma == t(Sigma)))

    return(Sigma)
}

### BELOW IS DONE AND IN RMD FILE

#' Generated coefficients for y in latent variable model
#'
#' @param p Integer or numeric; the number of features that will be observed in
#' x plus the number of latent variables (each corresponding to a cluster).
#' @param k_unblocked Integer or numeric; the number of features in X that
#' will have nonzero coefficients in the true model for y among those features 
#' not generated from the n_clusters latent variables (called "weak signal" 
#' features in the simulations from Faletto and Bien 2022). The coefficients on
#' these features will be determined by beta_low.
#' @param beta_low Integer or numeric; the maximum coefficient in the
#' model for y among the k_unblocked features in X not generated from the
#' latent variables. The coefficients of the features will be
#' beta_low/sqrt(1:k_unblocked).
#' @param beta_high Integer or numeric; the coefficient used for all
#' sig_blocks latent variables that have nonzero coefficients in the true
#' model for y.
#' @param nblocks Integer or numeric; the number of latent variables that were
#' generated, each of which will be associated with an observed cluster in X.
#' @param sig_blocks Integer or numeric; the number of generated latent
#' features that will have nonzero coefficients in the true model for y (all of
#' them will have coefficient beta_latent). In particular, the first sig_blocks
#' latent variables will have coefficient beta_latent, and the remaining nblocks
#' - sig_blocks features will have coefficient 0. Must be less than or equal to
#' n_clusters.
#' @param block_size Integer or numeric; for each of the n_blocks latent
#' variables, the covariance matrix will include the original latent feature
#' plus block_size - 1 noisy proxies that are correlated with the latent
#' variable.
#' @return A named list with the following elements: \item{beta}{A vector of
#' length `p` containing the coefficients for the true model for y. All entries
#' will equal 0 except for the sig_blocks latent variables that will have
#' coefficient beta_high and the k_unblocked independent features with
#' coefficient determined by beta_low.} \item{blocked_dgp_vars}{An integer
#' vector of length sig_blocks containing the indices of the features
#' corresponding to the latent features that will have nonzero coefficient
#' beta_high in the true model for y.} \item{sig_unblocked_vars}{An integer
#' vector of length k_unblocked containing the indices of the observed features
#' that are independent of the blocked features and have coefficient beta_low in
#' the true model for y. If k_unblocked = 0, this will just be NA.}
#' \item{insig_blocked_vars}{An integer vector containing the indices of the
#' features corresponding to the latent features that will have coefficient 0 in
#' the true model for y. If nblocks=0, this will just be NA.}
#' \item{latent_vars}{An integer vector of length nblocks containing the indices
#' of all of the latent features.}
#' @author Gregory Faletto, Jacob Bien
#' @references Faletto, G., & Bien, J. (2022). Cluster Stability Selection.
#' \emph{arXiv preprint arXiv:2201.00494}.
#' \url{https://arxiv.org/abs/2201.00494}.
makeCoefficients <- function(p, k_unblocked, beta_low, beta_high, nblocks,
    sig_blocks, block_size){

    # Check inputs
    stopifnot(k_unblocked >= 0)
    stopifnot(sig_blocks <= nblocks)
    stopifnot(p >= nblocks*block_size + k_unblocked)
    stopifnot(sig_blocks >= 0)

    # Initialize beta
    beta <- numeric(p)

    # identify indices of first coefficient in each significant block (these
    # features will have coefficient beta_high)
    latent_vars <- NA
    if(nblocks >= 1){
        latent_vars <- as.integer(((0:(nblocks - 1))*block_size + 1))

        stopifnot(all(latent_vars) %in% 1:p)
        stopifnot(all(latent_vars) %in% 1:(block_size*nblocks))
        stopifnot(length(unique(latent_vars)) == nblocks)
        stopifnot(length(latent_vars) == nblocks)
    }

    blocked_dgp_vars <- latent_vars[1:sig_blocks]

    stopifnot(sig_blocks == length(blocked_dgp_vars))
    
    beta[blocked_dgp_vars] <- beta_high

    # identify remaining coefficients in blocks (which ought to be set to 0)
    insig_blocked_vars <- NA

    if(nblocks >= 1){
        insig_blocked_vars <- setdiff(1:(block_size*nblocks), blocked_dgp_vars)
        stopifnot(all(beta[insig_blocked_vars] == 0))
    }
    # find significant unblocked variables (if applicable) and fill in
    # coefficients
    sig_unblocked_vars <- NA

    if(k_unblocked > 0){
        # Range of weak signal coefficients
        beta_lows <- beta_low/sqrt(1:k_unblocked)
        sig_unblocked_vars <- (nblocks*block_size + 1):
            (nblocks*block_size + k_unblocked)
        sig_unblocked_vars <- as.integer(sig_unblocked_vars)

        stopifnot(length(sig_unblocked_vars) == k_unblocked)
        stopifnot(length(unique(sig_unblocked_vars)) == k_unblocked)
        stopifnot(all(sig_unblocked_vars) %in% 1:p)

        beta[sig_unblocked_vars] <- beta_lows
    }

    stopifnot(length(intersect(blocked_dgp_vars, sig_unblocked_vars)) == 0)
    stopifnot(length(intersect(sig_unblocked_vars, insig_blocked_vars)) == 0)
    stopifnot(length(intersect(blocked_dgp_vars, insig_blocked_vars)) == 0)

    stopifnot(length(insig_blocked_vars) + length(blocked_dgp_vars) ==
        nblocks*block_size)

    stopifnot(sig_blocks + length(insig_blocked_vars) + k_unblocked <= p)

    stopifnot(sum(beta != 0) == sig_blocks + k_unblocked)
    stopifnot(is.numeric(beta) | is.integer(beta))

    return(list(beta=beta, blocked_dgp_vars=blocked_dgp_vars,
        sig_unblocked_vars=sig_unblocked_vars,
        insig_blocked_vars=insig_blocked_vars, latent_vars=latent_vars))
}

### BELOW IS DONE AND IN RMD FILE

#' Generate observed and latent variables along with conditional mean
#'
#' @param n Integer or numeric; the number of observations to generate. (The
#' generated X and Z will have n rows, and the generated y and mu will have
#' length n.)
#' @param p Integer or numeric; the number of observed features (the generated X
#' will have p columns).
#' @param beta A numeric or integer vector of length `p` + sig_blocks containing
#' the coefficients for the true model for y.
#' @param Sigma A (`p` + n_blocks) x (`p` + n_blocks) numeric matrix
#' representing the covariance matrix for the latent features, the associated
#' proxies, and the remaining features.
#' @param blocked_dgp_vars An integer vector of length sig_blocks containing the
#' indices of the features corresponding to the latent features that have
#' nonzero coefficient beta_high in the true model for y.
#' @param latent_vars An integer vector of length n_blocks containing the
#' indices of all of the latent features.
#' @param block_size Integer or numeric; for each of the n_blocks latent
#' variables, X will contain block_size noisy proxies that are correlated with
#' the latent variable.
#' @param n_blocks Integer or numeric; the number of latent variables to
#' generate, each of which will be associated with an observed cluster in X.
#' Must be at least 1. Default is 1.
#' @param snr Integer or numeric; the signal-to-noise ratio of the response
#' y. If sigma_eps_sq is not specified, the variance of the noise in y will be
#' calculated using the formula sigma_eps_sq = sum(mu^2)/(n * snr). Only one of
#' snr and sigma_eps_sq must be specified. Default is NA.
#' @param sigma_eps_sq Integer or numeric; the variance on the noise added
#' to y. Only one of snr and sigma_eps_sq must be specified. Default is NA.
#' @return A named list with the following elements: \item{X}{An `n` x `p`
#' numeric matrix containing the observed proxies for the latent variables as
#' well as the observed unblocked (iid) variables.} \item{mu}{A length `n`
#' numeric vector; the expected response given X, Z, and the true
#' coefficient vector (equal to y minus the added noise).} \item{z}{An `n` x
#' n_blocks numeric matrix containing the n_blocks latent variables. Note that
#' (X, z) is multivariate Gaussian.} \item{sd}{Numeric; the standard deviation
#' of the noise added to mu to get y (calculated either from snr or
#' sigma_eps_sq).}
#' @author Gregory Faletto, Jacob Bien
genMuXZSd <- function(n, p, beta, Sigma, blocked_dgp_vars,
    latent_vars, block_size, n_blocks=1, snr=NA, sigma_eps_sq=NA){
    # Check inputs

    stopifnot(length(blocked_dgp_vars) <= n_blocks)
    stopifnot(nrow(Sigma) == p + n_blocks)
    stopifnot(ncol(Sigma) == p + n_blocks)
    
    if(any(!is.na(sigma_eps_sq))){
        stopifnot(is.numeric(sigma_eps_sq) | is.integer(sigma_eps_sq))
        stopifnot(length(sigma_eps_sq) == 1)
        stopifnot(sigma_eps_sq >= 0)
    } else{
        if(any(is.na(snr))){
            stop("Must provide one of snr or sigma_eps_sq")
        }
        stopifnot(is.numeric(snr) | is.integer(snr))
        stopifnot(length(snr) == 1)
        stopifnot(snr > 0)
    }

    stopifnot(length(beta) == p + n_blocks)
    stopifnot(all(beta[blocked_dgp_vars] != 0))

    stopifnot(length(latent_vars) == n_blocks)

    x <- MASS::mvrnorm(n=n, mu=rep(0, p + n_blocks), Sigma=Sigma)

    stopifnot(length(beta) == ncol(x))

    mu <- as.numeric(x %*% beta)

    # Remove true blocked signal feature from each block from x now that I've
    # generated mu
    if(n_blocks > 0){
        z <- matrix(as.numeric(NA), nrow=n, ncol=n_blocks)
        stopifnot(length(latent_vars) > 0)
    } else{
        z <- NA
        stopifnot(length(latent_vars) == 0)
    }
    
    if(length(latent_vars) > 0){
        z[, 1:n_blocks] <- x[, latent_vars]
    }
    
    x <- x[, setdiff(1:(p + n_blocks), latent_vars)]

    # If SNR is null, use sigma_eps_sq
    if(!is.na(sigma_eps_sq)){
        sd <- sqrt(sigma_eps_sq)
    }else{
        sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 /(n * sigma^2)
    }

    # Check output

    stopifnot(length(mu) == n)

    stopifnot(nrow(x) == n)
    stopifnot(ncol(x) == p)

    if(any(!is.na(z))){
        stopifnot(nrow(z) == n)
        stopifnot(ncol(z) == n_blocks)
    }

    stopifnot(is.numeric(sd) | is.integer(sd))
    stopifnot(length(sd) == 1)
    stopifnot(!is.na(sd))
    stopifnot(sd >= 0)

    return(list(X=x, mu=mu, z=z, sd=sd))
}

### BELOW IS DONE AND IN RMD FILE

#' Create design matrix of cluster representatives from matrix of raw features
#' using results of css function
#'
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param weighting Character; determines how to calculate the weights to
#' combine features from the selected clusters into weighted averages, called
#' cluster representatives. Must be one of "sparse", "weighted_avg", or
#' "simple_avg'. For "sparse", all the weight is put on the most frequently
#' selected individual cluster member (or divided equally among all the clusters
#' that are tied for the top selection proportion if there is a tie). For
#' "weighted_avg", the weight used for each cluster member is calculated in
#' proportion to the individual selection proportions of each feature. For
#' "simple_avg", each cluster member gets equal weight regardless of the
#' individual feature selection proportions (that is, the cluster representative
#' is just a simple average of all the cluster members). See Faletto and Bien
#' (2022) for details. Default is "weighted_avg".
#' @param cutoff Numeric; css will return only those clusters with selection
#' proportions equal to at least cutoff. Must be between 0 and 1. Default is 0
#' (in which case all clusters are returned in decreasing order of selection
#' proportion).
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.) Default is 1.
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' max_num_clusts is ignored).
#' @param newx A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to generate the design matrix of cluster
#' representatives. Must contain the same features (in the same
#' number of columns) as the X matrix provided to css, and if the columns of
#' newx are labeled, the names must match the variable names provided to css.
#' newx may be omitted if train_inds were provided to css to set aside
#' observations for model estimation. If this is the case, then when newx is
#' omitted formCssDesign will return a design matrix of cluster representatives
#' formed from the train_inds observations from the matrix X provided to css.
#' (If no train_inds were provided to css, newX must be provided to
#' formCssDesign.) Default is NA.
#' @return A design matrix with the same number of rows as newx (or the 
#' train_inds provided to css) where the columns are the constructed cluster
#' representatives.
#' @author Gregory Faletto, Jacob Bien
#' @references Faletto, G., & Bien, J. (2022). Cluster Stability Selection.
#' \emph{arXiv preprint arXiv:2201.00494}.
#' \url{https://arxiv.org/abs/2201.00494}.
formCssDesign <- function(css_results, weighting="weighted_avg", cutoff=0,
    min_num_clusts=1, max_num_clusts=NA, newx=NA){

    # Check inputs
    ret <- checkFormCssDesignInputs(css_results, weighting, cutoff,
        min_num_clusts, max_num_clusts, newx)

    newx <- ret$newx
    max_num_clusts <- ret$max_num_clusts

    rm(ret)

    n <- nrow(newx)
    p <- ncol(newx)

    # Get the names of the selected clusters and the weights for the features
    # within each cluster, according to the provided weighting rule
    weights <- getSelectedClusters(css_results, weighting, cutoff,
        min_num_clusts, max_num_clusts)$weights

    n_sel_clusts <- length(weights)

    # Form matrix of cluster representatives of selected clusters
    X_clus_reps <- matrix(rep(as.numeric(NA), n*n_sel_clusts), nrow=n,
        ncol=n_sel_clusts)
    colnames(X_clus_reps) <- rep(as.character(NA), n_sel_clusts)

    for(i in 1:n_sel_clusts){
        clust_i_name <- names(weights)[i]

        stopifnot(length(clust_i_name) == 1)
        stopifnot(clust_i_name %in% names(weights))

        colnames(X_clus_reps)[i] <- clust_i_name

        clust_i <- css_results$clusters[[clust_i_name]]

        stopifnot(length(clust_i) >= 1)
        stopifnot(all(clust_i) %in% 1:p)

        weights_i <- weights[[clust_i_name]]

        stopifnot(length(clust_i) == length(weights_i))

        if(length(weights_i) > 1){
            X_clus_reps[, i] <- newx[, clust_i] %*% weights_i
        } else{
            X_clus_reps[, i] <- newx[, clust_i]*weights_i
        }
    }

    # Check output
    stopifnot(all(!is.na(X_clus_reps)))
    stopifnot(ncol(X_clus_reps) == n_sel_clusts)
    stopifnot(nrow(X_clus_reps) == n)

    return(X_clus_reps)
}

### BELOW IS DONE AND IN RMD FILE

#' From css output, obtain names of selected clusters and selection proportions,
#' indices of all selected features, and weights of individual cluster members
#'
#' If cutoff is too high for at least min_num_clusts clusters to be selected,
#' then it will be lowered until min_num_clusts can be selected. After that, if
#' the cutoff is too low such that more than max_num_clusts are selected, then
#' the cutoff will be increased until no more than max_num_clusts are selected.
#' Note that because clusters can have tied selection proportions, it is
#' possible that the number of selected clusters will be strictly lower than
#' max_num_clusts or strictly greater than min_num_clusts. In fact, it is
#' possible that both cutoffs won't be able to be satisfied simulteaneously,
#' even if there is a strictly positive difference between max_num_clusts and
#' min_num_clusts. If this occurs, max_num_clusts will take precedence over
#' min_num_clusts. getSelectedClusters will throw an error if the provided
#' inputs don't allow it to select any clusters. 
#' 
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param weighting Character; determines how to calculate the weights for
#' individual features within the selected clusters. Only those features with
#' nonzero weight within the selected clusters will be returned. Must be one of
#' "sparse", "weighted_avg", or "simple_avg'. For "sparse", all the weight is
#' put on the most frequently
#' selected individual cluster member (or divided equally among all the clusters
#' that are tied for the top selection proportion if there is a tie). For
#' "weighted_avg", only the features within a selected cluster that were
#' themselves selected on at least one subsample will have nonzero weight. For
#' "simple_avg", each cluster member gets equal weight regardless of the
#' individual feature selection proportions (that is, all cluster members within
#' each selected cluster will be returned.). See Faletto and Bien (2022) for
#' details.
#' @param cutoff Numeric; getCssSelections will select and return only of those
#' clusters with selection proportions equal to at least cutoff. Must be between
#' 0 and 1.
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.)
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) If NA, max_num_clusts is ignored.
#' @return A named list with the following elements: \item{selected_clusts}{A
#' named numeric vector containing the selection proportions for the selected
#' clusters. The name of each entry is the name of the corresponding cluster.}
#' \item{selected_feats}{A named integer vector; the indices of the features
#' with nonzero weights from all of the selected clusters.} \item{weights}{A
#' named list of the same length as the number of selected clusters. Each list
#' element weights[[j]] is a numeric vector of the weights to use for the jth
#' selected cluster, and it has the same name as the cluster it corresponds
#' to.}
#' @author Gregory Faletto, Jacob Bien
getSelectedClusters <- function(css_results, weighting, cutoff, min_num_clusts,
    max_num_clusts){
    # Check input
    stopifnot(class(css_results) == "cssr")

    # Eliminate clusters with selection proportions below cutoff
    clus_sel_props <- colMeans(css_results$clus_sel_mat)

    # Get selected clusters
    selected_clusts <- clus_sel_props[clus_sel_props >= cutoff]
    B <- nrow(css_results$feat_sel_mat)

    # Check that selected_clusts has length at least min_num_clusts
    while(length(selected_clusts) < min_num_clusts){
        cutoff <- cutoff - 1/B
        selected_clusts <- clus_sel_props[clus_sel_props >= cutoff]
    }

    # Check that selected_clusts has length at most max_num_clusts
    if(!is.na(max_num_clusts)){
        n_clusters <- ncol(css_results$clus_sel_mat)
        while(length(selected_clusts) > max_num_clusts){
            cutoff <- cutoff + 1/B
            if(cutoff > 1){
                break
            }
            # Make sure we don't reduce to a selected set of size 0
            if(any(clus_sel_props >= cutoff)){
                selected_clusts <- clus_sel_props[clus_sel_props >= cutoff]
            } else{
                break
            }
        }
    }

    stopifnot(length(selected_clusts) >= 1)

    clust_names <- names(selected_clusts)

    n_sel_clusts <- length(selected_clusts)

    # Check that n_sel_clusts is as expected, and throw warnings or an error if
    # not
    checkSelectedClusters(n_sel_clusts, min_num_clusts, max_num_clusts,
        max(clus_sel_props))
    
    ### Get selected features from selected clusters
    clusters <- css_results$clusters
    stopifnot(all(clust_names %in% names(clusters)))

    # Get a list of weights for all of the selected clusters
    weights <- getAllClustWeights(css_results, selected_clusts, weighting)

    # Get selected features from each cluster (those features with nonzero
    # weights)
    selected_feats <- integer()
    for(i in 1:n_sel_clusts){
        clus_i_name <- clust_names[i]
        clust_i <- clusters[[clus_i_name]]
        weights_i <- weights[[i]]
        selected_feats <- c(selected_feats, clust_i[weights_i != 0])
    }

    feat_names <- colnames(css_results$feat_sel_mat)

    names(selected_feats) <- feat_names[selected_feats]

    # Check output (already checked weights wihin getAllClustWeights)

    checkGetSelectedClustersOutput(selected_clusts, selected_feats,
        weights, n_clusters=length(clusters), p=ncol(css_results$feat_sel_mat))

    return(list(selected_clusts=selected_clusts,
        selected_feats=selected_feats, weights=weights))
}

### BELOW IS DONE AND IN RMD FILE

#' Identify selection prototypes from selected clusters
#'
#' Takes in list of selected clusters and returns an integer vector of the
#' indices of the features that were most frequently selected from each cluster
#'
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param selected_clusts A list of integer vectors; each vector must contain
#' the indices of features in a cluster.
#' @return An integer vector (of length equal to the number of clusters) of the
#' indices of the feature prototypes (the features from each cluster that were
#' selected the most often individually by the base method in cluster stability
#' selection). In the case of a tie, the tie is broken by choosing the feature
#' most correlated with the response in the full data set provided to css.
#' @author Gregory Faletto, Jacob Bien
getSelectionPrototypes <- function(css_results, selected_clusts){
    
    # Check inputs
    stopifnot(class(css_results) == "cssr")

    stopifnot(is.list(selected_clusts))
    n_selected_clusts <- length(selected_clusts)
    stopifnot(n_selected_clusts >= 1)
    stopifnot(all(lengths(selected_clusts) >= 1))

    prototypes <- rep(as.integer(NA), n_selected_clusts)
    for(i in 1:n_selected_clusts){
        clust_i <- selected_clusts[[i]]
        sel_props_i <- colMeans(css_results$feat_sel_mat)[clust_i]
        proto_i <- clust_i[sel_props_i == max(sel_props_i)]
        stopifnot(length(proto_i) >= 1)
        if(length(proto_i) > 1){
            if(is.numeric(css_results$y) | is.integer(css_results$y)){
                # Break tie by looking at marginal correlations
                corrs_i <- stats::cor(css_results$X[, proto_i], css_results$y)[, 1]
                corrs_i <- abs(corrs_i)
                proto_i <- proto_i[corrs_i == max(corrs_i)]
            }
        }
        # If there is still a tie, break by choosing the first feature of those
        # remaining
        prototypes[i] <- proto_i[1]
        names(prototypes)[i] <- colnames(css_results$X)[proto_i]
    }

    # Check output

    stopifnot(is.integer(prototypes))
    stopifnot(all(!is.na(prototypes)))
    stopifnot(length(prototypes) == length(unique(prototypes)))

    return(prototypes)
}

### BELOW IS DONE AND IN RMD FILE

#' Automated estimation of model size
#'
#' This function is uses the lasso with cross-validation to estimate the
#' model size. Before using the lasso, in each cluster all features will be
#' dropped from X except the one feature with the highest marginal correlation
#' with y, as in the protolasso (Reid and Tibshirani 2016).
#' 
#' @param X An n x p numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the p >= 2 features/predictors.
#' @param y A length-n numeric vector containing the responses; `y[i]` is the
#' response corresponding to observation `X[i, ]`. (Note that for the css
#' function, y does not have to be a numeric response, but for this function,
#' the underlying selection procedure is the lasso, so y must be a real-valued
#' response.)
#' @param clusters A named list where each entry is an integer vector of indices
#' of features that are in a common cluster, as in the output of css.
#' (The length of list clusters is equal to the number of clusters.) All
#' identified clusters must be non-overlapping, and all features must appear in
#' exactly one cluster (any unclustered features should be in their own
#' "cluster" of size 1).
#' @return An integer; the estimated size of the model. The minimum returned
#' value is 1, even if the lasso with cross-validation chose no features.
#' @author Gregory Faletto, Jacob Bien
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}.
#' @export
getModelSize <- function(X, y, clusters){

    stopifnot(is.matrix(X) | is.data.frame(X))

    # Check if x is a matrix; if it's a data.frame, convert to matrix.
    if(is.data.frame(X)){
        p <- ncol(X)

        X <- stats::model.matrix(~ ., X)
        X <- X[, colnames(X) != "(Intercept)"]

        if(length(clusters) > 0 & (p != ncol(X))){
            stop("When stats::model.matrix converted the provided data.frame X to a matrix, the number of columns changed (probably because the provided data.frame contained a factor variable with at least three levels). Please convert X to a matrix yourself using model.matrix and provide cluster assignments according to the columns of the new matrix.")
        }
    }

    stopifnot(is.matrix(X))
    stopifnot(all(!is.na(X)))
    stopifnot(is.numeric(X) | is.integer(X))
    n <- nrow(X)

    # Since the model size will be determined by cross-validation, the provided
    # y must be real-valued (this should be checked internally in other
    # functions before getModelSize is called, but this check is here just in
    # case).
    if(!is.numeric(y) & !is.integer(y)){
        stop("getModelSize is trying to determine max_num_clusts using the lasso with cross-validation, but the y provided to getModelSize was not real-valued.")
    }
    stopifnot(length(y) == n)

    # Check clusters argument
    clusters <- checkCssClustersInput(clusters)

    ### Format clusters into a list where all features are in exactly one
    # cluster (any unclustered features are put in their own "cluster" of size
    # 1).
    clust_names <- as.character(NA)
    if(!is.null(names(clusters)) & is.list(clusters)){
        clust_names <- names(clusters)
    }

    clusters <- formatClusters(clusters, p=ncol(X),
        clust_names=clust_names)$clusters

    X_size <- X

    if(length(clusters) > 0){
        # Create modified design matrix by dropping all but one feature from
        # each cluster
        feats_to_drop <- integer()
        for(i in 1:length(clusters)){
            cluster_i <- clusters[[i]]
            if(length(cluster_i) > 1){
                feat_to_keep <- identifyPrototype(cluster_i, X, y)
                feats_to_drop <- c(feats_to_drop, setdiff(cluster_i,
                    feat_to_keep))
            }
        }
        if(length(feats_to_drop) > 0){
            X_size <- X_size[, -1*feats_to_drop]
        }
    }

    size_results <- glmnet::cv.glmnet(x=X_size, y=y, family="gaussian")
    coefs <- as.numeric(glmnet::coef.glmnet(size_results, s="lambda.1se"))

    # Number of nonzero coefficients (subtract one in order to ignore intercept)
    size <- length(coefs[coefs != 0]) - 1

    # Check output
    stopifnot(is.numeric(size) | is.integer(size))
    stopifnot(!is.na(size))
    stopifnot(length(size) == 1)
    stopifnot(size == round(size))

    return(as.integer(max(size, 1)))
}

### BELOW IS DONE AND IN RMD FILE

#' Calculate weights for each cluster member of all of the selected clusters.
#' 
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param sel_clusters A named numeric vector containing the selection
#' proportions for the selected clusters. The name of each entry is the name
#' of the corresponding cluster.
#' @param weighting Character; determines how to calculate the weights for
#' individual features within the selected clusters. Only those features with
#' nonzero weight within the selected clusters will be returned. Must be one of
#' "sparse", "weighted_avg", or "simple_avg'. For "sparse", all the weight is
#' put on the most frequently
#' selected individual cluster member (or divided equally among all the clusters
#' that are tied for the top selection proportion if there is a tie). For
#' "weighted_avg", only the features within a selected cluster that were
#' themselves selected on at least one subsample will have nonzero weight. For
#' "simple_avg", each cluster member gets equal weight regardless of the
#' individual feature selection proportions (that is, all cluster members within
#' each selected cluster will be returned.). See Faletto and Bien (2022) for
#' details.
#' @return A named list of the same length as sel_clusters of numeric vectors.
#' weights[[j]] is the weights to use for the jth selected cluster, and it has
#' the same name as the cluster it corresponds to.
#' @author Gregory Faletto, Jacob Bien
getAllClustWeights <- function(css_results, sel_clusters, weighting){

    # Check inputs
    stopifnot(class(css_results) == "cssr")

    stopifnot(is.numeric(sel_clusters))
    p_ret <- length(sel_clusters)
    stopifnot(length(unique(names(sel_clusters))) == p_ret)
    stopifnot(p_ret > 0)

    checkWeighting(weighting)

    # Get selection proportions and clusters
    feat_sel_props <- colMeans(css_results$feat_sel_mat)

    p <- length(feat_sel_props)
    stopifnot(p >= p_ret)

    clusters <- css_results$clusters
    stopifnot(all(names(sel_clusters) %in% names(clusters)))

    # Identify weights
    weights <- list()

    for(j in 1:p_ret){
        # Find the members of the cluster feature j is a member of
        cluster_j <- clusters[[names(sel_clusters)[j]]]
        # Get the weights for this cluster and add them to the list
        weights[[j]] <- getClustWeights(cluster_j, weighting, feat_sel_props)
    }

    # Add names to weights
    names(weights) <- names(sel_clusters)

    # Check output

    stopifnot(length(weights) == p_ret)
    stopifnot(is.list(weights))

    for(i in 1:p_ret){
        stopifnot(length(clusters[[names(sel_clusters)[i]]]) ==
            length(weights[[i]]))
        stopifnot(all(weights[[i]] >= 0))
        stopifnot(all(weights[[i]] <= 1))
        stopifnot(abs(sum(weights[[i]]) - 1) < 10^(-6))
    }
    return(weights)
}

### BELOW IS DONE AND IN RMD FILE

#' Calculate weights for members of a cluster using selection proportions
#'
#' Given a cluster of features, the selection proportions for each cluster
#' member, and a specified weighting scheme, calculate the appropriate weights
#' for the cluster.
#' @param cluster_i An integer vector containing the indices of the members
#' of a cluster.
#' @param weighting Character; determines how to calculate the weights for
#' individual features within the selected clusters. Only those features with
#' nonzero weight within the selected clusters will be returned. Must be one of
#' "sparse", "weighted_avg", or "simple_avg'. For "sparse", all the weight is
#' put on the most frequently
#' selected individual cluster member (or divided equally among all the clusters
#' that are tied for the top selection proportion if there is a tie). For
#' "weighted_avg", only the features within a selected cluster that were
#' themselves selected on at least one subsample will have nonzero weight. For
#' "simple_avg", each cluster member gets equal weight regardless of the
#' individual feature selection proportions (that is, all cluster members within
#' each selected cluster will be returned.). See Faletto and Bien (2022) for
#' details.
#' @param feat_sel_props A numeric vector of selection proportions corresponding
#' to each of the p features.
#' @return A numeric vector of the same length as cluster_i containing the
#' weights corresponding to each of the features in cluster_i. The weights
#' will all be nonnegative and sum to 1.
#' @author Gregory Faletto, Jacob Bien
getClustWeights <- function(cluster_i, weighting, feat_sel_props){

    stopifnot(is.integer(cluster_i) | is.numeric(cluster_i))
    stopifnot(all(cluster_i == round(cluster_i)))
    n_weights <- length(cluster_i)
    stopifnot(length(unique(cluster_i)) == n_weights)

    p <- length(feat_sel_props)
    stopifnot(all(cluster_i %in% 1:p))

    # Get the selection proportions of each cluster member
    sel_props <- feat_sel_props[cluster_i]

    stopifnot(all(sel_props >= 0))
    stopifnot(all(sel_props <= 1))

    weights_i <- rep(as.numeric(NA), n_weights)

    # Weighted or simple average?
    if(weighting == "sparse"){
        # Sparse cluster stability selection: All features in cluster with
        # selection proportion equal to the max
        # for the cluster get equal weight; rest of cluster gets 0 weight
        if(sum(sel_props) == 0){
            weights_i <- rep(1/n_weights, n_weights)
        } else{
            maxes <- sel_props==max(sel_props)

            stopifnot(sum(maxes) > 0)
            stopifnot(sum(maxes) <= n_weights)

            weights_i <- rep(0, n_weights)
            weights_i[maxes] <- 1/sum(maxes)
        }
    } else if(weighting == "weighted_avg"){
        # Get weights for weighted average
        if(sum(sel_props) == 0){
            weights_i <- rep(1/n_weights, n_weights)
        } else{
            weights_i <- sel_props/sum(sel_props)
        }
    } else if(weighting == "simple_avg"){
        weights_i <- rep(1/n_weights, n_weights)
    } else{
        stop("weighting must be one of sparse, simple_avg, or weighted_avg")
    }

    stopifnot(abs(sum(weights_i) - 1) < 10^(-6))
    stopifnot(length(weights_i) == n_weights)
    stopifnot(length(weights_i) >= 1)
    stopifnot(all(weights_i >= 0))
    stopifnot(all(weights_i <= 1))

    return(weights_i)
}

### BELOW IS DONE AND IN RMD FILE

#' Absolute value of sample correlation between two vectors
#'
#' Calculates the absolute value of correlation of t and y. If either input has
#' only one unique value, returns 0 by definition.
#' @param t A numeric or integer vector.
#' @param y A numeric or integer vector; must have the same length as t.
#' @return A numeric vector of the same length as cluster_i containing the
#' weights corresponding to each of the features in cluster_i. The weights
#' will all be nonnegative and sum to 1.
#' @author Gregory Faletto, Jacob Bien
corFunction <- function(t, y){
    # Check inputs
    stopifnot(is.numeric(t) | is.integer(t))
    stopifnot(is.numeric(y) | is.integer(y))
    stopifnot(length(t) == length(y))
    if(length(unique(t)) == 1){
        return(0)
    }
    if(length(unique(y)) == 1){
        warning("The second argument to corFunction only had one unique entry")
        return(0)
    }
    return(abs(stats::cor(t, y)))
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that inputs to several functions are as expected,
#' and modify inputs if needed
#'
#' @param newx A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to generate the design matrix of cluster
#' representatives. Must contain the same features (in the same
#' number of columns) as the X matrix provided to css, and if the columns of
#' newX are labeled, the names must match the variable names provided to css.
#' @param css_X The X matrix provided to css, as in the output of the css
#' function (after having been coerced from a data.frame to a matrix by css if
#' needed).
#' @return A named list with the following elements. \item{feat_names}{A 
#' character vector containing the column names of newx (if the provided newx
#' had column names). If the provided newx did not have column names, feat_names
#' will be NA.} \item{newx}{The provided newx matrix, coerced from a data.frame
#' to a matrix if the provided newx was a data.frame.}
#' @author Gregory Faletto, Jacob Bien
checkXInputResults <- function(newx, css_X){

    # Check if x is a matrix; if it's a data.frame, convert to matrix.
    if(is.data.frame(newx)){
        newx <- stats::model.matrix(~ ., newx)
        newx <- newx[, colnames(newx) != "(Intercept)"]
    }

    feat_names <- as.character(NA)
    if(!is.null(colnames(newx))){
        feat_names <- colnames(newx)
        stopifnot(identical(feat_names, colnames(css_X)))
    } else{
        # In this case, newx has no column names, so same better be true of
        # css_X
        if(!is.null(colnames(css_X))){
            warning("New X provided had no variable names (column names) even though the X provided to css did.")
        }
    }

    stopifnot(is.matrix(newx))
    stopifnot(all(!is.na(newx)))

    n <- nrow(newx)
    p <- ncol(newx)
    stopifnot(p >= 2)
    if(length(feat_names) > 1){
        stopifnot(length(feat_names) == p)
        stopifnot(!("(Intercept)" %in% feat_names))
    } else{
        stopifnot(is.na(feat_names))
    }

    colnames(newx) <- character()

    # Confirm that newx matches css_results$X
    if(p != ncol(css_X)){
        err <- paste("Number of columns in newx must match number of columns from matrix provided to css. Number of columns in new provided X: ",
            p, ". Number of columns in matrix provided to css: ", ncol(css_X),
            ".", sep="")
        stop(err)
    }
    if(length(feat_names) != 1 & all(!is.na(feat_names))){
        if(!identical(feat_names, colnames(css_X))){
            stop("Provided feature names for newx do not match feature names provided to css")
        }
    }

    return(list(feat_names=feat_names, newx=newx))
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that inputs to the function css are as expected,
#' and modify inputs if needed
#'
#' @param X An n x p numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' p >= 2 features/predictors.
#' @param y The response; can be anything that takes the form of an
#' n-dimensional vector, with the ith entry corresponding to the ith row of X.
#' Typically (and for default fitfun = cssLasso), y should be an n-dimensional
#' numeric vector.
#' @param lambda A tuning parameter or set of tuning parameters that may be used
#' by the feature selection method fitfun. In the default case when
#' fitfun = cssLasso, lambda should be a numeric: the penalty to use for each
#' lasso fit. (css does not require lambda to be any particular object because
#' for a user-specified feature selection method fitfun, lambda can be an
#' arbitrary object. See the description of fitfun below.)
#' @param clusters A list of integer vectors; each vector should contain the 
#' indices of a cluster of features (a subset of 1:p). (If there is only one
#' cluster, clusters can either be a list of length 1 or an integer vector.)
#' All of the provided clusters must be non-overlapping. Every feature not
#' appearing in any cluster will be assumed to be unclustered (that is, they
#' will be treated as if they are in a "cluster" containing only themselves). If
#' clusters is a list of length 0 (or a list only containing clusters of length
#' 1), then css() returns the same results as stability selection (so the
#' returned feat_sel_mat will be identical to clus_sel_mat). Names for the
#' clusters will be needed later; any clusters that are not given names in the
#' provided list will be given names automatically by css. Default is list() (so
#' no clusters are specified).
#' @param fitfun A function; the feature selection function used on each
#' subsample by cluster stability selection. This can be any feature selection
#' method; the only requirement is that it accepts the arguments (and only the
#' arguments) X, y, and lambda and returns an integer vector that is a subset of
#' 1:p. For example, fitfun could be best subset selection or forward stepwise
#' selection or LARS and lambda could be the desired model size; or fitfun could be the
#' elastic net and lambda could be a length-two vector specifying lambda and
#' alpha. Default is cssLasso, an implementation of lasso (relying on the R
#' package glmnet), where lambda must be a positive numeric specifying the L1
#' penalty for the lasso.
#' @param sampling_type A character vector; either "SS" or "MB". For "MB",
#' all B subsamples are drawn randomly (as proposed by Meinshausen and Bühlmann
#' 2010). For "SS", in addition to these B subsamples, the B complementary pair
#' subsamples will be drawn as well (see Faletto and Bien 2022 or Shah and
#' Samworth 2013 for details). Default is "SS", and "MB" is not supported yet.
#' @param B Integer or numeric; the number of subsamples. Note: For
#' sampling.type=="MB" the total number of subsamples will be `B`; for
#' sampling_type="SS" the number of subsamples will be `2*B`. Default is 100
#' for sampling_type="MB" and 50 for sampling_type="SS".
#' @param prop_feats_remove Numeric; if prop_feats_remove is greater than 0,
#' then on each subsample, each feature is randomly dropped from the design
#' matrix that is provided to fitfun with probability prop_feats_remove
#' (independently across features). That is, in a typical subsample,
#' prop_feats_remove*p features will be dropped (though this number will vary).
#' This is similar in spirit (but distinct from) extended stability selection
#' (Beinrucker et. al. 2016); see their paper for some of the benefits of
#' dropping features (besides increasing computational speed and decreasing
#' memory requirements). For sampling_type="SS", the features dropped in
#' each complementary pair of subsamples are identical in order to ensure that
#' the theoretical guarantees of Faletto and Bien (2022) are retained within
#' each individual pair of subsamples. (Note that this feature is not
#' investigated either theoretically or in simulations by Faletto and Bien
#' 2022). Must be between 0 and 1. Default is 0.
#' @param train_inds Optional; an integer or numeric vector containing the
#' indices of observations in X and y to set aside for model training by the
#' function getCssPreds after feature selection. (This will only work if y is
#' real-valued, because getCssPreds using ordinary least squares regression to
#' generate predictions.) If train_inds is not provided, all of the observations
#' in the provided data set will be used for feature selection.
#' @param num_cores Optional; an integer. If using parallel processing, the
#' number of cores to use for parallel processing (num_cores will be supplied
#' internally as the mc.cores argument of parallel::mclapply).
#' @return A named list with the following elements: \item{feat_names}{A 
#' character vector containing the column names of X (if the provided X
#' had column names). If the provided X did not have column names, feat_names
#' will be NA.} \item{X}{The provided X, converted to a matrix if it was
#' originally provided as a data.frame, and with feature names removed if they
#' had been provided.}\item{clusters}{A list of integer vectors; each vector
#' will contain the indices of a cluster of features. Any duplicated clusters
#' provided in the input will be removed.}
#' @author Gregory Faletto, Jacob Bien
checkCssInputs <- function(X, y, lambda, clusters, fitfun, sampling_type, B,
    prop_feats_remove, train_inds, num_cores){

    stopifnot(is.matrix(X) | is.data.frame(X))

    clust_names <- as.character(NA)
    if(!is.null(names(clusters)) & is.list(clusters)){
        clust_names <- names(clusters)
    }

    # Check if x is a matrix; if it's a data.frame, convert to matrix.
    if(is.data.frame(X)){
        p <- ncol(X)

        X <- stats::model.matrix(~ ., X)
        X <- X[, colnames(X) != "(Intercept)"]

        if(length(clusters) > 0 & (p != ncol(X))){
            stop("When stats::model.matrix converted the provided data.frame X to a matrix, the number of columns changed (probably because the provided data.frame contained a factor variable with at least three levels). Please convert X to a matrix yourself using model.matrix and provide cluster assignments according to the columns of the new matrix.")
        }
    }

    stopifnot(is.matrix(X))
    stopifnot(all(!is.na(X)))

    feat_names <- as.character(NA)
    if(!is.null(colnames(X))){
        feat_names <- colnames(X)
    }

    n <- nrow(X)
    p <- ncol(X)

    if(!is.null(colnames(X))){
        feat_names <- colnames(X)
    }

    stopifnot(p >= 2)
    if(length(feat_names) > 1){
        stopifnot(length(feat_names) == p)
    } else{
        stopifnot(is.na(feat_names))
    }

    colnames(X) <- character()

    stopifnot(length(y) == n)
    # Intentionally don't check y or lambda further to allow for flexibility--these
    # inputs should be checked within fitfun.

    # Check clusters argument
    clusters <- checkCssClustersInput(clusters)

    ### Format clusters into a list where all features are in exactly one
    # cluster (any unclustered features are put in their own "cluster" of size
    # 1).
    clusters <- formatClusters(clusters, p=p, clust_names=clust_names)$clusters

    stopifnot(class(fitfun) == "function")
    stopifnot(length(fitfun) == 1)
    if(!identical(formals(fitfun), formals(cssLasso))){
        err_mess <- paste("fitfun must accept arguments named X, y, and lambda. Detected arguments to fitfun:",
            paste(names(formals(fitfun)), collapse=", "))
        stop(err_mess)
    }

    checkSamplingType(sampling_type)
    checkB(B)
    checkPropFeatsRemove(prop_feats_remove, p)

    stopifnot(is.numeric(train_inds) | is.integer(train_inds))
    if(length(train_inds) > 0){
        stopifnot(all(!is.na(train_inds)))
        stopifnot(all(round(train_inds) == train_inds))
        stopifnot(length(train_inds) == length(unique(train_inds)))
        stopifnot(all(train_inds >= 1))
        stopifnot(all(train_inds) <= n)
        stopifnot(length(train_inds) <= n - 2)
        stopifnot(length(train_inds) >= 1)
    }

    stopifnot(length(num_cores) == 1)
    stopifnot(is.integer(num_cores) | is.numeric(num_cores))
    stopifnot(!is.na(num_cores))
    stopifnot(num_cores == round(num_cores))
    stopifnot(num_cores >= 1)
    stopifnot(num_cores <= parallel::detectCores())

    return(list(feat_names=feat_names, X=X, clusters=clusters))
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that inputs to the function getCssPreds are as
#' expected, and modify inputs if needed.
#'
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param testX A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to generate predictions. Must contain the same
#' features (in the same number of columns) as the matrix provided to css.
#' @param weighting Character; determines how to calculate the weights to
#' combine features from the selected clusters into weighted averages, called
#' cluster representatives. Must be one of "sparse", "weighted_avg", or
#' "simple_avg'. For "sparse", all the weight is put on the most frequently
#' selected individual cluster member (or divided equally among all the clusters
#' that are tied for the top selection proportion if there is a tie). For
#' "weighted_avg", the weight used for each cluster member is calculated in
#' proportion to the individual selection proportions of each feature. For
#' "simple_avg", each cluster member gets equal weight regardless of the
#' individual feature selection proportions (that is, the cluster representative
#' is just a simple average of all the cluster members). See Faletto and Bien
#' (2022) for details. Default is "weighted_avg".
#' @param cutoff Numeric; getCssPreds will make use only of those clusters with
#' selection proportions equal to at least cutoff. Must be between 0 and 1.
#' Default is 0 (in which case either all clusters are used, or max_num_clusts
#' are used, if max_num_clusts is specified).
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.) Default is 1.
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Default is NA (in which case
#' max_num_clusts is ignored).
#' @param trainX A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to estimate the linear model from the selected
#' clusters. trainX is only necessary to provide if no train_inds were
#' designated in the css function call to set aside observations for model
#' estimation (though even if train_inds was provided, trainX and trianY will be
#' used for model estimation if they are both provided to getCssPreds). Must 
#' contain the same features (in the same number of columns) as the matrix 
#' provided to css, and if the columns of trainX are labeled, the names must
#' match the variable names provided to css. Default is NA (in which case
#' getCssPreds uses the observations from the train_inds that were provided to
#' css to estimate a linear model).
#' @param trainY The response corresponding to trainX. Must be a real-valued
#' response (unlike in the general css setup) because predictions will be
#' generated by an ordinary least squares model. Must have the same length as
#' the number of rows of trainX. Like trainX, only needs to be provided if no
#' observations were set aside for model estimation by the parameter train_inds
#' in the css function call. Default is NA (in which case getCssPreds uses the
#' observations from the train_inds that were provided to css).
#' @return A named list with the following elements: \item{trainXProvided}{
#' Logical; indicates whether a valid trainX input was provided.} \item{trainX}{
#' The provided trainX matrix, coerced from a data.frame to a matrix if the
#' provided trainX was a data.frame. (If a valid trainX was not provided, this
#' output simply passes whatever was provided as trainX.)} \item{testX}{The
#' provided testX matrix, coerced from a data.frame to a matrix if the provided
#' testX was a data.frame.} \item{feat_names}{A character vector containing the
#' column names of testX (if the provided testX had column names). If the
#' provided testX did not have column names, feat_names will be NA.}
#' \item{max_num_clusts}{The provided max_num_clusts, coerced to an integer if
#' needed, and coerced to be less than or equal to the total number of clusters
#' from the output of css_results.}
#' @author Gregory Faletto, Jacob Bien
checkGetCssPredsInputs <- function(css_results, testX, weighting, cutoff,
    min_num_clusts, max_num_clusts, trainX, trainY){
    # Check inputs
    stopifnot(class(css_results) == "cssr")

    check_results <- checkNewXProvided(trainX, css_results)

    trainX <- check_results$newX
    trainXProvided <- check_results$newXProvided

    rm(check_results)

    n_train <- nrow(trainX)

    if(trainXProvided){
        if(all(!is.na(trainY)) & length(trainY) > 1){
            stopifnot(is.numeric(trainY))
            stopifnot(n_train == length(trainY))
        } else{
            if(length(css_results$train_inds) == 0){
                stop("css was not provided with indices to set aside for model training (train_inds), so must provide both trainX and trainY in order to generate predictions")
            }
            trainXProvided <- FALSE
            warning("trainX provided but no trainY provided; instead, training model using the train_inds observations provided to css to set aside for model training.")
        }
    } else{
        if(length(css_results$train_inds) == 0){
            stop("css was not provided with indices to set aside for model training (train_inds), so must provide both trainX and trainY in order to generate predictions")
        }
        if(all(!is.na(trainY)) & length(trainY) > 1){
            warning("trainY provided but no trainX provided; instead, training model using the train_inds observations provided to css to set aside for model training.")
        }
    }

    results <- checkXInputResults(testX, css_results$X)

    testX <- results$newx
    feat_names <- results$feat_names

    if(all(!is.na(feat_names))){
        stopifnot(length(feat_names) == ncol(testX))
        stopifnot(!("(Intercept)" %in% feat_names))
        colnames(testX) <- feat_names
    }

    rm(results)

    n <- nrow(testX)
    p <- ncol(testX)

    stopifnot(n >= 1)
    stopifnot(p == ncol(trainX))
    if(!is.null(colnames(trainX)) & is.null(colnames(testX))){
        warning("Column names were provided for trainX but not for testX (are you sure they both contain identical features in the same order?)")
    }
    if(is.null(colnames(trainX)) & !is.null(colnames(testX))){
        warning("Column names were provided for testX but not for trainX (are you sure they both contain identical features in the same order?)")
    }
    if(!is.null(colnames(trainX)) & !is.null(colnames(testX))){
        stopifnot(all(colnames(trainX) == colnames(testX)))
    }

    checkCutoff(cutoff)
    checkWeighting(weighting)
    checkMinNumClusts(min_num_clusts, p, length(css_results$clusters))
    max_num_clusts <- checkMaxNumClusts(max_num_clusts, min_num_clusts, p,
        length(css_results$clusters))

    return(list(trainXProvided=trainXProvided, trainX=trainX, testX=testX,
        feat_names=feat_names, max_num_clusts=max_num_clusts))

}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument cutoff to several functions is
#' as expected
#'
#' @param cutoff Numeric; only those clusters with selection proportions equal
#' to at least cutoff will be selected by cluster stability selection. Must be
#' between 0 and 1.
#' @author Gregory Faletto, Jacob Bien
checkCutoff <- function(cutoff){
    stopifnot(is.numeric(cutoff) | is.integer(cutoff))
    stopifnot(length(cutoff) == 1)
    stopifnot(!is.na(cutoff))
    stopifnot(cutoff >= 0)
    stopifnot(cutoff <= 1)
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument min_num_clusts to several 
#' functions is as expected
#'
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.)
#' @param p The number of features; since this is an upper bound on the number
#' of clusters of features, it is also an upper bound on min_num_clusts.
#' @param n_clusters The number of clusters; note that this is an upper bound
#' on min_num_clusts
#' @author Gregory Faletto, Jacob Bien
checkMinNumClusts <- function(min_num_clusts, p, n_clusters){
    stopifnot(length(min_num_clusts) == 1)
    stopifnot(is.numeric(min_num_clusts) | is.integer(min_num_clusts))
    stopifnot(!is.na(min_num_clusts))
    stopifnot(min_num_clusts == round(min_num_clusts))
    stopifnot(min_num_clusts >= 1)
    stopifnot(min_num_clusts <= p)
    stopifnot(min_num_clusts <= n_clusters)
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument max_num_clusts to several 
#' functions is as expected
#'
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) Can be NA, in which case
#' max_num_clusts will be ignored.
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.) max_num_clusts must be at least as
#' large as min_num_clusts.
#' @param p The number of features; since this is an upper bound on the number
#' of clusters of features, it is also an upper bound on max_num_clusts.
#' @param n_clusters The number of clusters; note that this is an upper bound
#' on max_num_clusts
#' @return The provided max_num_clusts, coerced to an integer if needed, and
#' coerced to be less than or equal to the total number of clusters.
#' @author Gregory Faletto, Jacob Bien
checkMaxNumClusts <- function(max_num_clusts, min_num_clusts, p, n_clusters){
    stopifnot(length(max_num_clusts) == 1)
    if(!is.na(max_num_clusts)){
        stopifnot(is.numeric(max_num_clusts) | is.integer(max_num_clusts))
        stopifnot(max_num_clusts == round(max_num_clusts))
        stopifnot(max_num_clusts >= 1)
        stopifnot(max_num_clusts <= p)
        max_num_clusts <- as.integer(min(n_clusters, max_num_clusts))
        stopifnot(max_num_clusts >= min_num_clusts)
    }
    return(max_num_clusts)
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument weighting to several 
#' functions is as expected
#'
#' @param weighting Character; determines how to calculate the weights to
#' combine features from the selected clusters into weighted averages, called
#' cluster representatives. Must be one of "sparse", "weighted_avg", or
#' "simple_avg'.
#' @author Gregory Faletto, Jacob Bien
checkWeighting <- function(weighting){
    stopifnot(length(weighting)==1)
    stopifnot(!is.na(weighting))
    if(!is.character(weighting)){
        stop("Weighting must be a character")
    }
    if(!(weighting %in% c("sparse", "simple_avg", "weighted_avg"))){
        stop("Weighting must be a character and one of sparse, simple_avg, or weighted_avg")
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument sampling_type to several 
#' functions is as expected
#'
#' @param sampling_type A character vector; either "SS" or "MB". "MB" is not
#' supported yet.
#' @author Gregory Faletto, Jacob Bien
checkSamplingType <- function(sampling_type){
    stopifnot(is.character(sampling_type))
    stopifnot(length(sampling_type) == 1)
    stopifnot(!is.na(sampling_type))
    stopifnot(sampling_type %in% c("SS", "MB"))
    if(sampling_type == "MB"){
        stop("sampling_type MB is not yet supported (and isn't recommended anyway)")
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument prop_feats_remove to several 
#' functions is as expected
#'
#' @param prop_feats_remove Numeric; proportion of features that are dropped on
#' each subsample. Must be between 0 and 1.
#' @param p The number of features; must be greater than 2 if prop_feats_remove
#' is greater than 0.
#' @author Gregory Faletto, Jacob Bien
checkPropFeatsRemove <- function(prop_feats_remove, p){
    stopifnot(length(prop_feats_remove) == 1)
    stopifnot(is.numeric(prop_feats_remove) | is.integer(prop_feats_remove))
    stopifnot(!is.na(prop_feats_remove))
    stopifnot(prop_feats_remove >= 0 & prop_feats_remove < 1)
    if(prop_feats_remove > 0){
        # Make sure p is at least 2 or else this doesn't make sense
        stopifnot(p >= 2)
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument B to several functions is as
#' expected
#'
#' @param B Integer or numeric; the number of subsamples. Note: For
#' sampling.type=="MB" the total number of subsamples will be `B`; for
#' sampling_type="SS" the number of subsamples will be `2*B`.
#' @author Gregory Faletto, Jacob Bien
checkB <- function(B){
    stopifnot(length(B) == 1)
    stopifnot(is.numeric(B) | is.integer(B))
    stopifnot(!is.na(B))
    stopifnot(B == round(B))
    stopifnot(B > 0)
    if(B < 10){
        warning("Small values of B may lead to poor results.")
    } else if (B > 2000){
        warning("Large values of B may require long computation times.")
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument clusters to several functions is
#' as expected
#'
#' @param clusters A named list where each entry is an integer vector of indices
#' of features that are in a common cluster, as in the output of css or
#' formatClusters. (The length of list clusters is equal to the number of
#' clusters.) All identified clusters must be non-overlapping, and all features
#' must appear in exactly one cluster (any unclustered features should be in
#' their own "cluster" of size 1).
#' @param p The number of features; must be at least as large as the number of
#' clusters.
#' @author Gregory Faletto, Jacob Bien
checkClusters <- function(clusters, p){
    stopifnot(is.list(clusters))
    stopifnot(all(lengths(clusters) >= 1))
    stopifnot(all(!is.na(clusters)))

    n_clusters <- length(clusters)

    stopifnot(n_clusters == length(unique(clusters)))
    stopifnot(n_clusters <= p)
    stopifnot(!is.null(names(clusters)))
    stopifnot(is.character(names(clusters)))
    stopifnot(all(!is.na(names(clusters)) & names(clusters) != ""))
    stopifnot(length(unique(names(clusters))) == n_clusters)

    all_clustered_feats <- integer()
    for(i in 1:n_clusters){
        stopifnot(is.integer(clusters[[i]]))
        all_clustered_feats <- c(all_clustered_feats, clusters[[i]])
    }

    stopifnot(length(all_clustered_feats) == p)
    stopifnot(length(unique(all_clustered_feats)) == p)
    stopifnot(all(all_clustered_feats <= p))
    stopifnot(all(all_clustered_feats >= 1))
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the argument y to several functions is
#' as expected
#'
#' @param y Numeric response vector.
#' @param n Number of observations of covariates; should match length of y.
#' @author Gregory Faletto, Jacob Bien
checkY <- function(y, n){
    stopifnot(all(!is.na(y)))
    stopifnot(is.numeric(y) | is.integer(y))
    stopifnot(length(unique(y)) > 1)
    stopifnot(length(n) == 1)
    stopifnot(!is.na(n))
    stopifnot(is.numeric(n) | is.integer(n))
    stopifnot(n == round(n))
    stopifnot(n > 0)
    stopifnot(n == length(y))
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the inputs to cssLasso are as expected. 
#'
#' @param X A design matrix containing the predictors. (Note that we don't need
#' to check X very much, because X will have already been checked by the
#' function checkCssInputs when it was provided to css.)
#' @param y A numeric vector containing the response.
#' @param lambda Numeric; a nonnegative number for the lasso penalty to use
#' on each subsample. (For now, only one lambda value can be provided to
#' cssLasso; in the future, we plan to allow for multiple lambda values to be
#' provided to cssLasso, as described in Faletto and Bien 2022.)
#' @author Gregory Faletto, Jacob Bien
checkCssLassoInputs <- function(X, y, lambda){

    n <- nrow(X)
    p <- ncol(X)

    if(!is.numeric(y)){
        stop("For method cssLasso, y must be a numeric vector.")
    }
    if(is.matrix(y)){
        stop("For method cssLasso, y must be a numeric vector (inputted y was a matrix).")
    }
    if(n != length(y)){
        stop("For method cssLasso, y must be a vector of length equal to nrow(X).")
    }
    if(length(unique(y)) <= 1){
        stop("Subsample with only one unique value of y detected--for method cssLasso, all subsamples of y of size floor(n/2) must have more than one unique value.")
    }
    if(!is.numeric(lambda) & !is.integer(lambda)){
        stop("For method cssLasso, lambda must be a numeric.")
    }
    if(any(is.na(lambda))){
        stop("NA detected in provided lambda input to cssLasso")
    }
    if(length(lambda) != 1){
        stop("For method cssLasso, lambda must be a numeric of length 1.")
    }
    if(lambda < 0){
        stop("For method cssLasso, lambda must be nonnegative.")
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the outputs of the provided feature selection
#' method are as required. 
#'
#' @param selected An integer vector; the indices of the features selected by
#' the lasso.
#' @param p The total number of observed features; all selected features must be
#' in 1:p.
#' @param feats_on_subsamp Integer; the indices of the features considered by
#' the feature selection method. All selected features must be among these
#' features.
#' @author Gregory Faletto, Jacob Bien
checkCssLoopOutput <- function(selected, p, feats_on_subsamp){
    if(!exists("selected")){
        stop("The provided feature selection method fitfun failed to return anything on (at least) one subsample")
    }
    if(!is.integer(selected) & !is.numeric(selected)){
        stop("The provided feature selection method fitfun failed to return an integer or numeric vector on (at least) one subsample")
    }
    if(any(is.na(selected))){
        stop("The provided feature selection method fitfun returned a vector containing NA values on (at least) one subsample")
    }
    if(!all(selected == round(selected))){
        stop("The provided feature selection method fitfun failed to return a vector of valid (integer) indices on (at least) one subsample")
    }
    if(length(selected) != length(unique(selected))){
        stop("The provided feature selection method fitfun returned a vector of selected features containing repeated indices on (at least) one subsample")
    }
    if(length(selected) > p){
        stop("The provided feature selection method fitfun returned a vector of selected features longer than p on (at least) one subsample")
    }
    if(length(selected) > 0){
        if(max(selected) > p){
            stop("The provided feature selection method fitfun returned a vector of selected features containing an index greater than ncol(X) on (at least) one subsample")
        }
        if(min(selected) <= 0){
            stop("The provided feature selection method fitfun returned a vector of selected features containing a non-positive index on (at least) one subsample")
        }
    }
    if("try-error" %in% class(selected) |
        "error" %in% class(selected) | "simpleError" %in% class(selected) |
        "condition" %in% class(selected)){
        stop("The provided feature selection method fitfun returned an error on (at least) one subsample")
    }
    if(!all(selected %in% feats_on_subsamp)){
        stop("The provided feature selection method somehow selected features that were not provided for it to consider.")
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that the new X matrix provided to getCssDesign or
#' getCssPreds matches the characteristics of the X that was provided to css.
#'
#' @param trainX A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix). Must contain
#' the same features (in the same number of columns) as the X matrix provided to
#' css, and if the columns of trainX are labeled, the names must match the
#' variable names provided to css. trainX may be omitted if train_inds were
#' provided to css to set aside observations.
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @return A named list with the following elements: \item{newX}{If trainX was
#' provided, this is the provided trainX matrix, coerced from a data.frame to a
#' matrix if the provided trainX was a data.frame. If trainX was not provided,
#' this is a matrix made up of the training indices provided to css in the
#' train_inds argument.} \item{newXProvided}{Logical; indicates whether a valid
#' trainX input was provided.}
#' @author Gregory Faletto, Jacob Bien
checkNewXProvided <- function(trainX, css_results){
    newXProvided <- FALSE

    if(all(!is.na(trainX)) & length(trainX) > 1){
        newXProvided <- TRUE
        trainX <- checkXInputResults(trainX, css_results$X)$newx
        
        n_train <- nrow(trainX)
        stopifnot(n_train > 1)
    } else{
        if(length(css_results$train_inds) == 0){
            stop("css was not provided with indices to set aside for model training (train_inds), so must provide new X in order to generate a design matrix")
        }
        trainX <- css_results$X[css_results$train_inds, ]
    } 
    stopifnot(is.matrix(trainX))
    stopifnot(is.numeric(trainX) | is.integer(trainX))
    stopifnot(all(!is.na(trainX)))
    stopifnot(ncol(trainX) >= 2)

    return(list(newX=trainX, newXProvided=newXProvided))
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to check operations within getSelectedClusters function
#'
#' @param n_sel_clusts The number of selected clusters; should be constrained
#' by min_num_clusts and max_num_clusts (though it may not be possible to
#' satisfy both constraints simulteneously, in which case a warning will be
#' thrown).
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.)
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.) If NA, max_num_clusts is ignored.
#' @param max_sel_prop Numeric; the maximum selection proportion observed for 
#' any cluster.
#' @author Gregory Faletto, Jacob Bien
checkSelectedClusters <- function(n_sel_clusts, min_num_clusts, max_num_clusts,
    max_sel_prop){
    if(n_sel_clusts == 0){
        err <- paste("No clusters selected with this cutoff (try a cutoff below the maximum cluster selection proportion, ",
            max_sel_prop, ")", sep="")
        stop(err)
    }

    stopifnot(n_sel_clusts >= 1)

    # It may be impossible to get at least min_num_clusts or at most
    # max_num_clusts; if so, give a warning
    if(n_sel_clusts < min_num_clusts){
        warn <- paste("Returning fewer than min_num_clusts = ", min_num_clusts,
            " clusters because decreasing the cutoff any further would require returning more than max_num_clusts = ",
            max_num_clusts, " clusters", sep="")
        warning(warn)
    }
    if(!is.na(max_num_clusts)){
        if(n_sel_clusts > max_num_clusts){
            warn <- paste("Returning more than max_num_clusts = ",
                max_num_clusts,
                " clusters because increasing the cutoff any further would require returning 0 clusters",
                sep="")
            warning(warn)
        }
    }
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to ensure that the inputs to formatClusters are as expected
#'
#' @param clusters Either an integer vector of a list of integer vectors; each
#' vector should contain the indices of a cluster of features. (If there is only
#' one cluster, clusters can either be a list of length 1 or simply an integer
#' vector.) If clusters is specified then R is ignored.
#' @param p integer or numeric; the numbe of features in x (should match 
#' ncol(x), if x is provided)
#' @param clust_names A character vector of the names of the clusters in clusters.
#' @param get_prototypes Logical: if TRUE, will identify prototype from each
#' cluster (the feature from each cluster that is most correlated with the
#' response) for the protolasso. In this case, x and y must be provided.
#' @param x n x p numeric matrix; design matrix. Only needs to be provided if
#' get_prototypes is TRUE.
#' @param y Numeric response vector; only needs to be provided if get_prototypes
#' is TRUE. Note: in general, the css function does not require y to be a
#' numeric vector, because the provided fitfun could use a different form of y
#' (for example, a categorical response variable). However, y must be numeric in
#' order to provide prototypes because the prototypes are determined using the
#' correlation between cluster members (columns of x) and y.
#' @param R Numeric p x p matrix; not currently used. Entry ij contains the 
#' "substitutive value" of feature i for feature j (diagonal must consist of
#' ones, all entries must be between 0 and 1, and matrix must be symmetric)
#' @return A list of integer vectors; each vector will contain the indices of a
#' cluster of features. Any duplicated clusters provided in the input will be
#' removed.
#' @author Gregory Faletto, Jacob Bien
checkFormatClustersInput <- function(clusters, p, clust_names, 
    get_prototypes, x, y, R){

    if(any(is.na(clusters)) & any(is.na(R))){
        stop("Must specify one of clusters or R (or does one of these provided inputs contain NA?)")
    }

    stopifnot(is.integer(p) | is.numeric(p))
    stopifnot(length(p) == 1)
    stopifnot(p == round(p))
    stopifnot(!is.na(p))
    if(p > 0){
        stopifnot(p >= 2)
    }

    use_R <- FALSE
    if(any(is.na(clusters)) | length(clusters) == 0){
        if(all(!is.na(R))){
            stopifnot(is.matrix(R))
            stopifnot(all(dim(R) == p))
            stopifnot(all(diag(R) == 1))
            stopifnot(identical(R, t(R)))
            stopifnot(all(!is.na(R)))
            stopifnot(all(R %in% c(0, 1)))
            use_R <- TRUE
        }
    } else{
        stopifnot(!is.list(clusters) | all(lengths(clusters) >= 1))
        stopifnot(is.list(clusters) | length(clusters) >= 1)
        stopifnot(all(!is.na(clusters)))

        if(is.list(clusters) & length(clusters) > 0){
            for(i in 1:length(clusters)){
                stopifnot(length(clusters[[i]]) == length(unique(clusters[[i]])))
                stopifnot(all(!is.na(clusters[[i]])))
                stopifnot(all(clusters[[i]] >= 1))
                stopifnot(is.integer(clusters[[i]]))
            }

            stopifnot(length(clusters) == length(unique(clusters)))

            clusters <- clusters[!duplicated(clusters)]

            if(length(clusters) >= 2){
                # Check that clusters are non-overlapping
                for(i in 1:(length(clusters) - 1)){
                    for(j in (i+1):length(clusters)){
                        stopifnot(length(intersect(clusters[[i]],
                            clusters[[j]])) == 0)
                    }
                }
            }

            if(any(!is.na(clust_names))){
                stopifnot(length(clust_names) == length(clusters))
            }
        } else if(!is.list(clusters)){
            clusters_temp <- clusters
            clusters <- list()
            clusters[[1]] <- clusters_temp
            rm(clusters_temp)
        }
    }

    stopifnot(length(get_prototypes) == 1)
    stopifnot(is.logical(get_prototypes))

    if(any(!is.na(clust_names))){
        stopifnot(is.character(clust_names))
    }

    if(get_prototypes){
        stopifnot(all(!is.na(x)))
        stopifnot(is.matrix(x))

        n <- nrow(x)

        checkY(y, n)
    }

    if(use_R){
        # Determine clusters from R
        clusters <- list()

        for(i in 1:nrow(R)){
            clusters[[i]] <- as.integer(which(R[i, ] > 0))
            stopifnot(length(clusters[[i]]) == length(unique(clusters[[i]])))
            stopifnot(all(!is.na(clusters[[i]])))
            stopifnot(is.integer(clusters[[i]]))
        }

        clusters <- unique(clusters)
        stopifnot(is.list(clusters))

        if(length(clusters) >= 2){
            # Check that clusters are non-overlapping
            for(i in 1:(length(clusters) - 1)){
                for(j in (i+1):length(clusters)){
                    if(length(intersect(clusters[[i]], clusters[[j]])) != 0){
                        stop("Invalid R matrix with overlapping clusters (clusters must not be overlapping)")
                    }
                }
            }
        }
    }

    stopifnot(is.list(clusters))

    return(clusters)
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to check inputs to getClusterSelMatrix function
#'
#' @param clusters A named list where each entry is an integer vector of indices
#' of features that are in a common cluster, as in the output of formatClusters.
#' (The length of list clusters is equal to the number of clusters.) All
#' identified clusters must be non-overlapping, and all features must appear in
#' exactly one cluster (any unclustered features should be in their own
#' "cluster" of size 1).
#' @param res A binary integer matrix. res[i, j] = 1 if feature j was selected
#' on subsample i and equals 0 otherwise, as in the output of getSelMatrix.
#' (That is, each row is a selected set.)
#' @return The parameter B, corresponding to half of the subsamples for 
#' sampling_type "SS".
#' @author Gregory Faletto, Jacob Bien
checkGetClusterSelMatrixInput <- function(clusters, res){
    stopifnot(is.matrix(res))
    stopifnot(all(res %in% c(0, 1)))
    p <- ncol(res)
    stopifnot(nrow(res) > 0)
    checkClusters(clusters, p)
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to check that the inputs to formCssDesign are as expected
#'
#' @param css_results An object of class "cssr" (the output of the function
#' css).
#' @param weighting Character; determines how to calculate the weights to
#' combine features from the selected clusters into weighted averages, called
#' cluster representatives. Must be one of "sparse", "weighted_avg", or
#' "simple_avg'.
#' @param cutoff Numeric; css will return only those clusters with selection
#' proportions equal to at least cutoff. Must be between 0 and 1.
#' @param min_num_clusts Integer or numeric; the minimum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns fewer than
#' min_num_clusts clusters, the cutoff will be increased until at least
#' min_num_clusts clusters are selected.)
#' @param max_num_clusts Integer or numeric; the maximum number of clusters to
#' use regardless of cutoff. (That is, if the chosen cutoff returns more than
#' max_num_clusts clusters, the cutoff will be decreased until at most
#' max_num_clusts clusters are selected.)
#' @param newx A numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' the data that will be used to generate the design matrix of cluster
#' representatives. Must contain the same features (in the same
#' number of columns) as the X matrix provided to css, and if the columns of
#' newx are labeled, the names must match the variable names provided to css.
#' newx may be omitted if train_inds were provided to css to set aside
#' observations for model estimation. If this is the case, then when newx is
#' omitted formCssDesign will return a design matrix of cluster representatives
#' formed from the train_inds observations from the matrix X provided to css.
#' (If no train_inds were provided to css, newX must be provided to
#' formCssDesign.)
#' @return A named list with the following elements: \item{newx}{If newx was
#' provided, the provided newx matrix, coerced from a data.frame to a matrix if
#' needed. If newx was not provided, a matrix formed by the train_inds set
#' aside in the original function call to css.} \item{max_num_clusts}{The
#' provided max_num_clusts, coerced to an integer if needed, and coerced to be
#' less than or equal to the total number of clusters.}
#' @author Gregory Faletto, Jacob Bien
checkFormCssDesignInputs <- function(css_results, weighting, cutoff,
    min_num_clusts, max_num_clusts, newx){    
    stopifnot(class(css_results) == "cssr")

    if(length(newx) == 1){
        if(is.na(newx)){
            if(length(css_results$train_inds) == 0){
                stop("If css was not provided with indices to set aside for model training, then newx must be provided to formCssDesign")
            }
            newx <- css_results$X[css_results$train_inds, ]
            # feat_names <- colnames(newx)
        } else{
            results <- checkXInputResults(newx, css_results$X)

            newx <- results$newx
            # feat_names <- results$feat_names

            rm(results)
        }
    } else{
        results <- checkXInputResults(newx, css_results$X)

        newx <- results$newx
        # feat_names <- results$feat_names

        rm(results)
    }

    p <- ncol(newx)

    checkCutoff(cutoff)
    checkWeighting(weighting)
    checkMinNumClusts(min_num_clusts, p, length(css_results$clusters))
    max_num_clusts <- checkMaxNumClusts(max_num_clusts, min_num_clusts, p,
        length(css_results$clusters))

    return(list(newx=newx, max_num_clusts=max_num_clusts))
}

#' Helper function to check that output of getSelectedClusters is as expected
#'
#' @param selected_clusts A named numeric vector containing the selection
#' proportions for the selected clusters. The name of each entry is the name of
#' the corresponding cluster.
#' @param selected_feats A named integer vector; the indices of the features
#' with nonzero weights from all of the selected clusters.
#' @param weights A named list of the same length as the number of selected
#' clusters. Each list element weights[[j]] is a numeric vector of the weights
#' to use for the jth selected cluster, and it has the same name as the cluster
#' it corresponds to.
#' @param n_clusters Integer; the number of clusters in the data (upper bound
#' for the length of selected_clusts)
#' @param p Integer; number of features in the data (all selected_feats should
#' be in 1:p)
#' @author Gregory Faletto, Jacob Bien
checkGetSelectedClustersOutput <- function(selected_clusts, selected_feats,
    weights, n_clusters, p){
    stopifnot(is.numeric(selected_clusts))
    stopifnot(all(selected_clusts >= 0))
    stopifnot(all(selected_clusts <= 1))
    stopifnot(length(selected_clusts) >= 1)
    stopifnot(length(selected_clusts) <= n_clusters)
    stopifnot(length(names(selected_clusts)) ==
        length(unique(names(selected_clusts))))
    stopifnot(!is.null(names(selected_clusts)))
    stopifnot(all(!is.na(names(selected_clusts)) &
        names(selected_clusts) != ""))
    stopifnot(length(names(selected_clusts)) == length(selected_clusts))
    stopifnot(is.integer(selected_feats))
    stopifnot(length(selected_feats) == length(unique(selected_feats)))
    stopifnot(all(selected_feats %in% 1:p))
    stopifnot(length(selected_clusts) <= length(selected_feats))
    stopifnot(identical(names(weights), names(selected_clusts)))
    stopifnot(length(weights) == length(selected_clusts))
}

### BELOW IS DONE AND IN RMD FILE

#' Helper function to confirm that clusters input to css is as expected
#'
#' @param clusters A list of integer vectors; each vector should contain the 
#' indices of a cluster of features (a subset of 1:p). (If there is only one
#' cluster, clusters can either be a list of length 1 or an integer vector.)
#' All of the provided clusters must be non-overlapping. Every feature not
#' appearing in any cluster will be assumed to be unclustered (that is, they
#' will be treated as if they are in a "cluster" containing only themselves). If
#' clusters is a list of length 0 (or a list only containing clusters of length
#' 1), then css() returns the same results as stability selection (so the
#' returned feat_sel_mat will be identical to clus_sel_mat). Names for the
#' clusters will be needed later; any clusters that are not given names in the
#' provided list will be given names automatically by css. Default is list() (so
#' no clusters are specified).
#' @return Same as the input, but all of the clusters will be coerced to
#' integers.
#' @author Gregory Faletto, Jacob Bien
checkCssClustersInput <- function(clusters){
    stopifnot(!is.na(clusters))
    if(is.list(clusters)){
        stopifnot(all(!is.na(clusters)))
        stopifnot(length(clusters) == length(unique(clusters)))

        if(length(clusters) > 0){
            for(i in 1:length(clusters)){
                stopifnot(length(clusters[[i]]) == length(unique(clusters[[i]])))
                stopifnot(all(!is.na(clusters[[i]])))
                stopifnot(is.integer(clusters[[i]]) | is.numeric(clusters[[i]]))
                stopifnot(all(clusters[[i]] == round(clusters[[i]])))
                stopifnot(all(clusters[[i]] >= 1))
                clusters[[i]] <- as.integer(clusters[[i]])
            }

            if(length(clusters) >= 2){
                # Check that clusters are non-overlapping
                for(i in 1:(length(clusters) - 1)){
                    for(j in (i+1):length(clusters)){
                        if(length(intersect(clusters[[i]], clusters[[j]])) != 0){
                            error_mes <- paste("Overlapping clusters detected; clusters must be non-overlapping. Overlapping clusters: ",
                                i, ", ", j, ".", sep="")
                            stop(error_mes)
                        }
                    }
                }
            }
        }
    } else{
        # If clusters is not a list, it should be a vector of indices of
        # features that are in the (one) cluster
        stopifnot(is.numeric(clusters) | is.integer(clusters))
        stopifnot(length(clusters) == length(unique(clusters)))
        stopifnot(all(!is.na(clusters)))
        stopifnot(is.integer(clusters) | is.numeric(clusters))
        stopifnot(all(clusters == round(clusters)))
        stopifnot(all(clusters >= 1))
        clusters <- as.integer(clusters)
    }
    return(clusters)
}

#' Generate randomly sampled data including noisy observations of latent
#' variables, where proxies differ in their relevance (noise level)
#'
#' Generate a data set including latent features Z, observed features X (which
#' may include noisy or noiseless observations of the latent features in Z),
#' an observed response y which is a linear model of features from Z and X as
#' well as independent mean zero noise, and mu (the responses from y without
#' the added noise). Data is generated in the same way as in the simulations
#' from Section 5.3 of Faletto and Bien (2022).
#' @param n Integer or numeric; the number of observations to generate. (The
#' generated X and Z will have n rows, and the generated y and mu will have
#' length n.)
#' @param p Integer or numeric; the number of features to generate. The
#' generated X will have p columns.
#' @param k_unclustered Integer or numeric; the number of features in X that
#' will have nonzero coefficients in the true model for y among those features 
#' not generated from the n_clusters latent variables (called "weak signal" 
#' features in the simulations from Faletto and Bien 2022). The coefficients on
#' these features will be determined by beta_unclustered.
#' @param cluster_size Integer or numeric; for each of the n_clusters latent
#' variables, X will contain cluster_size noisy proxies that are correlated with
#' the latent variable.
#' @param n_strong_cluster_vars Integer or numeric; among the cluster_size
#' proxies in each cluster, n_strong_cluster_vars will have a high covariance
#' (rho_high) with the latent variable and cluster_size - n_strong_cluster_vars
#' will have a low covariance (rho_low) with the latent variable.
#' @param n_clusters Integer or numeric; the number of latent variables to
#' generate, each of which will be associated with an observed cluster in X.
#' Must be at least 1. Default is 1.
#' @param sig_clusters Integer or numeric; the number of generated latent
#' features that will have nonzero coefficients in the true model for y (all of
#' them will have coefficient beta_latent). Must be less than or equal to
#' n_clusters. Default is 1.
#' @param rho_high Integer or numeric; the covariance of the "strong proxies" in
#' each cluster with the latent variable (and each other). Note that the
#' correlation between the "strong proxy" features in the cluster will be
#' rho_high/var. rho_high cannot equal 0 and must be at least as large as
#' rho_low. Default is 0.9.
#' @param rho_low Integer or numeric; the covariance of the "weak proxies" in
#' each cluster with the latent variable (and each other). Note that the
#' correlation between the "weak proxy" features in the cluster will be
#' rho_low/var. rho_low cannot equal 0 and must be no larger than rho_high.
#' Default is 0.5.
#' @param beta_latent Integer or numeric; the coefficient used for all
#' sig_clusters latent variables that have nonzero coefficients in the true
#' model for y. Can't equal 0. Default is 1.5.
#' @param beta_unclustered Integer or numeric; the maximum coefficient in the
#' model for y among the k_unclustered features in X not generated from the
#' latent variables. The coefficients of the features will be
#' beta_unclustered/sqrt(1:k_unclustered). Can't equal 0. Default is 1.
#' @param snr Integer or numeric; the signal-to-noise ratio of the response
#' y. If sigma_eps_sq is not specified, the variance of the noise in y will be
#' calculated using the formula sigma_eps_sq = sum(mu^2)/(n * snr). Only one of
#' snr and sigma_eps_sq must be specified. Default is NA.
#' @param sigma_eps_sq Integer or numeric; the variance on the noise added
#' to y. Only one of snr and sigma_eps_sq must be specified. Default is NA.
#' @return A list of the following elements. \item{X}{An n x p numeric matrix of
#' n observations from a p-dimensional multivariate normal distribution
#' generated using the specified parameters. The first n_clusters times
#' cluster_size features will be the clusters of features correlated with the
#' n_clusters latent variables. The next k_unclustered features will be the
#' "weak signal" features, and the remaining p - n_clusters*cluster_size -
#' k_unclustered features will be the unclustered noise features.} \item{y}{A
#' length n numeric vector; the response generated from X, the latent features
#' from Z, and the coefficient vector, along with additive noise.} \item{Z}{The
#' latent features; either a numeric vector (if n_clusters > 1) or a numeric
#' matrix (if n_clusters > 1). Note that (X, Z) is multivariate Gaussian.}
#' \item{mu}{A length `n` numeric vector; the expected response given X, Z, and
#' the true coefficient vector (equal to y minus the added noise).}
#' @author Gregory Faletto, Jacob Bien
#' @references
<<faletto2022>>
#' @export
genClusteredDataWeighted <- function(n, p, k_unclustered, cluster_size,
    n_strong_cluster_vars, n_clusters=1, sig_clusters=1, rho_high=0.9,
    rho_low=0.5, beta_latent=1.5, beta_unclustered=1, snr=as.numeric(NA),
    sigma_eps_sq=as.numeric(NA)){

    # Check inputs
    checkGenClusteredDataWeightedInputs(p, k_unclustered, cluster_size,
        n_strong_cluster_vars, n_clusters,  sig_clusters, rho_high, rho_low,
        beta_latent, beta_unclustered, snr, sigma_eps_sq)

    ret <- genZmuY(n=n, p=p, k_unclustered=k_unclustered,
        cluster_size=cluster_size, n_clusters=n_clusters,
        sig_clusters=sig_clusters, beta_latent=beta_latent,
        beta_unclustered=beta_unclustered, snr=snr, sigma_eps_sq=sigma_eps_sq)

    Z <- ret$Z
    y <- ret$y
    mu <- ret$mu
    other_X <- ret$other_X

    # Finally, generate clusters of proxies to complete X.
    noise_var_high <- getNoiseVar(rho_high)
    noise_var_low <- getNoiseVar(rho_low)

    # Create matrix of proxies
    Z <- as.matrix(Z)
    proxy_mat <- matrix(as.numeric(NA), n, n_clusters*cluster_size)
    for(i in 1:n_clusters){
        for(j in 1:n_strong_cluster_vars){
            proxy_mat[, (i - 1)*cluster_size + j] <- Z[, i] + rnorm(n,
                mean=0, sd=sqrt(noise_var_high))
        }
        for(j in (n_strong_cluster_vars + 1):cluster_size){
            proxy_mat[, (i - 1)*cluster_size + j] <- Z[, i] + rnorm(n,
                mean=0, sd=sqrt(noise_var_low))
        }
    }

    X <- cbind(proxy_mat, other_X)
    
    # Check output
    stopifnot(length(mu) == n)

    stopifnot(nrow(X) == n)
    stopifnot(ncol(X) == p)

    if(any(!is.na(Z))){
        stopifnot(nrow(Z) == n)
        stopifnot(ncol(Z) == n_clusters)
    }

    return(list(X=X, y=y, Z=Z, mu=mu))
}

#' Generate randomly sampled data including noisy observations of latent
#' variables, where proxies differ in their relevance (noise level)
#'
#' Generate a data set including latent features Z, observed features X (which
#' may include noisy or noiseless observations of the latent features in Z),
#' an observed response y which is a linear model of features from Z and X as
#' well as independent mean zero noise, and mu (the responses from y without
#' the added noise).
#' @param n Integer or numeric; the number of observations to generate. (The
#' generated X and Z will have n rows, and the generated y and mu will have
#' length n.)
#' @param p Integer or numeric; the number of features to generate. The
#' generated X will have p columns.
#' @param k_unclustered Integer or numeric; the number of features in X that
#' will have nonzero coefficients in the true model for y among those features 
#' not generated from the n_clusters latent variables (called "weak signal" 
#' features in the simulations from Faletto and Bien 2022). The coefficients on
#' these features will be determined by beta_unclustered.
#' @param cluster_size Integer or numeric; for each of the n_clusters latent
#' variables, X will contain cluster_size noisy proxies that are correlated with
#' the latent variable.
#' @param n_clusters Integer or numeric; the number of latent variables to
#' generate, each of which will be associated with an observed cluster in X.
#' Must be at least 1. Default is 1.
#' @param sig_clusters Integer or numeric; the number of generated latent
#' features that will have nonzero coefficients in the true model for y (all of
#' them will have coefficient beta_latent). Must be less than or equal to
#' n_clusters. Default is 1.
#' @param rho_high Integer or numeric; the maximum correlation of the proxies
#' each cluster with each other. Default is 1.
#' @param rho_low Integer or numeric; the minimum correlation of the proxies in
#' each cluster with each other. rho_low cannot equal 0 and must be no larger
#' than rho_high. Default is 0.5.
#' @param beta_latent Integer or numeric; the coefficient used for all
#' sig_clusters latent variables that have nonzero coefficients in the true
#' model for y. Can't equal 0. Default is 1.5.
#' @param beta_unclustered Integer or numeric; the maximum coefficient in the
#' model for y among the k_unclustered features in X not generated from the
#' latent variables. The coefficients of the features will be
#' beta_unclustered/sqrt(1:k_unclustered). Can't equal 0. Default is 1.
#' @param snr Integer or numeric; the signal-to-noise ratio of the response
#' y. If sigma_eps_sq is not specified, the variance of the noise in y will be
#' calculated using the formula sigma_eps_sq = sum(mu^2)/(n * snr). Only one of
#' snr and sigma_eps_sq must be specified. Default is NA.
#' @param sigma_eps_sq Integer or numeric; the variance on the noise added
#' to y. Only one of snr and sigma_eps_sq must be specified. Default is NA.
#' @return A list of the following elements. \item{X}{An n x p numeric matrix of
#' n observations from a p-dimensional multivariate normal distribution
#' generated using the specified parameters. The first n_clusters times
#' cluster_size features will be the clusters of features correlated with the
#' n_clusters latent variables. The next k_unclustered features will be the
#' "weak signal" features, and the remaining p - n_clusters*cluster_size -
#' k_unclustered features will be the unclustered noise features.} \item{y}{A
#' length n numeric vector; the response generated from X, the latent features
#' from Z, and the coefficient vector, along with additive noise.} \item{Z}{The
#' latent features; either a numeric vector (if n_clusters > 1) or a numeric
#' matrix (if n_clusters > 1). Note that (X, Z) is multivariate Gaussian.}
#' \item{mu}{A length `n` numeric vector; the expected response given X, Z, and
#' the true coefficient vector (equal to y minus the added noise).}
#' @author Gregory Faletto, Jacob Bien
#' @references
<<faletto2022>>
#' @export
genClusteredDataWeightedRandom <- function(n, p, k_unclustered, cluster_size,
    n_clusters=1, sig_clusters=1, rho_high=1, rho_low=0.5, beta_latent=1.5,
    beta_unclustered=1, snr=as.numeric(NA), sigma_eps_sq=as.numeric(NA)){

    # Check inputs
    checkGenClusteredDataWeightedRandomInputs(p, k_unclustered, cluster_size,
        n_clusters, sig_clusters, rho_high, rho_low,
        beta_latent, beta_unclustered, snr, sigma_eps_sq)

    ret <- genZmuY(n=n, p=p, k_unclustered=k_unclustered,
        cluster_size=cluster_size, n_clusters=n_clusters,
        sig_clusters=sig_clusters, beta_latent=beta_latent,
        beta_unclustered=beta_unclustered, snr=snr, sigma_eps_sq=sigma_eps_sq)

    Z <- ret$Z
    y <- ret$y
    mu <- ret$mu
    other_X <- ret$other_X

    # Finally, generate clusters of proxies to complete X.
    # Create matrix of proxies
    Z <- as.matrix(Z)
    proxy_mat <- matrix(as.numeric(NA), n, n_clusters*cluster_size)
    for(i in 1:n_clusters){
        for(j in 1:cluster_size){
            # Choose correlation at random
            rho_ij <- runif(n=1, min=rho_low, max=rho_high)
            # Get noise variance
            noise_var_ij <- getNoiseVar(rho_ij)
            proxy_mat[, (i - 1)*cluster_size + j] <- Z[, i] + rnorm(n,
                mean=0, sd=sqrt(noise_var_ij))
        }
    }

    X <- cbind(proxy_mat, other_X)
    
    # Check output
    stopifnot(length(mu) == n)

    stopifnot(nrow(X) == n)
    stopifnot(ncol(X) == p)

    if(any(!is.na(Z))){
        stopifnot(nrow(Z) == n)
        stopifnot(ncol(Z) == n_clusters)
    }

    return(list(X=X, y=y, Z=Z, mu=mu))
}

#' Generates Z, weak signal features in X, noise features in X, mu, and y
#' from provided parameters
#'
#' @param n Integer or numeric; the number of observations to generate. (The
#' generated X and Z will have n rows, and the generated y and mu will have
#' length n.)
#' @param p Integer or numeric; the number of features to generate. The
#' generated X will have p columns.
#' @param k_unclustered Integer or numeric; the number of features in X that
#' will have nonzero coefficients in the true model for y among those features 
#' not generated from the n_clusters latent variables (called "weak signal" 
#' features in the simulations from Faletto and Bien 2022). The coefficients on
#' these features will be determined by beta_unclustered. Must be at least 1.
#' @param cluster_size Integer or numeric; for each of the n_clusters latent
#' variables, X will contain cluster_size noisy proxies that are correlated with
#' the latent variable. Must be at least 2.
#' @param n_clusters Integer or numeric; the number of latent variables to
#' generate, each of which will be associated with an observed cluster in X.
#' Must be at least 1. Default is 1.
#' @param sig_clusters Integer or numeric; the number of generated latent
#' features that will have nonzero coefficients in the true model for y (all of
#' them will have coefficient beta_latent). Must be less than or equal to
#' n_clusters. Default is 1.
#' @param beta_latent Integer or numeric; the coefficient used for all
#' sig_clusters latent variables that have nonzero coefficients in the true
#' model for y. Can't equal 0. Default is 1.5.
#' @param beta_unclustered Integer or numeric; the maximum coefficient in the
#' model for y among the k_unclustered features in X not generated from the
#' latent variables. The coefficients of the features will be
#' beta_unclustered/sqrt(1:k_unclustered). Can't equal 0. Default is 1.
#' @param snr Integer or numeric; the signal-to-noise ratio of the response
#' y. If sigma_eps_sq is not specified, the variance of the noise in y will be
#' calculated using the formula sigma_eps_sq = sum(mu^2)/(n * snr). Only one of
#' snr and sigma_eps_sq must be specified. Default is NA.
#' @param sigma_eps_sq Integer or numeric; the variance on the noise added
#' to y. Only one of snr and sigma_eps_sq must be specified. Default is NA.
#' @return A list of the following elements. \item{Z}{The
#' latent features; either a numeric vector (if n_clusters > 1) or a numeric
#' matrix (if n_clusters > 1). Note that (X, Z) is multivariate Gaussian.}
#' item{mu}{A length `n` numeric vector; the expected response given X, Z, and
#' the true coefficient vector (equal to y minus the added noise).} \item{y}{A
#' length n numeric vector; the response generated from X, the latent features
#' from Z, and the coefficient vector, along with additive noise.} 
#' \item{other_X}{A numeric matrix of n observations from a multivariate normal
#' distribution generated using the specified parameters, containing the weak
#' signal features and the noise features that will eventually be in X. (The
#' only missing features are the proxies for the latent features Z.)
#' @author Gregory Faletto, Jacob Bien
#' @references
<<faletto2022>>
genZmuY <- function(n, p, k_unclustered, cluster_size, n_clusters, sig_clusters,
    beta_latent, beta_unclustered, snr, sigma_eps_sq){
    # Generate Z, weak signal features, and noise features (total of
    # p - n_clusters*(cluster_size - 1)) features)
    p_orig_feat_mat <- p - n_clusters*(cluster_size - 1)
    stopifnot(p_orig_feat_mat >= k_unclustered + n_clusters)
    orig_feat_mat <- matrix(stats::rnorm(n*p_orig_feat_mat), n, p_orig_feat_mat)
    # First n_clusters features are Z. Next k_unclustered features are weak
    # signal features. Any remaining features are noise features.
    Z <- orig_feat_mat[, 1:n_clusters]

    other_X <- orig_feat_mat[, (n_clusters + 1):p_orig_feat_mat]

    # Ready to create mu and y
    if(n_clusters > 1){
        if(sig_clusters > 1){
            mu <- Z[, 1:sig_clusters] %*% rep(beta_latent, sig_clusters)
            
        } else{
            mu <- Z[, 1:sig_clusters] * beta_latent
        }
    } else{
        mu <- Z*beta_latent
    }
    for(j in 1:k_unclustered){
        mu <- mu + beta_unclustered/sqrt(j)*other_X[, j]
    }
    mu <- as.numeric(mu)

    # If SNR is null, use sigma_eps_sq
    if(!is.na(sigma_eps_sq)){
        sd <- sqrt(sigma_eps_sq)
    }else{
        sd <- sqrt(sum(mu^2) / (n * snr)) # taking snr = ||mu||^2 /(n * sigma^2)
    }

    stopifnot(is.numeric(sd) | is.integer(sd))
    stopifnot(length(sd) == 1)
    stopifnot(!is.na(sd))
    stopifnot(sd >= 0)

    y <- as.numeric(mu + rnorm(n, mean=0, sd=sd))

    return(list(Z=Z, mu=mu, y=y, other_X=other_X))
}

#' Check inputs to genClusteredDataWeighted
#'
#' @param n Integer or numeric; the number of observations to generate. (The
#' generated X and Z will have n rows, and the generated y and mu will have
#' length n.)
#' @param p Integer or numeric; the number of features to generate. The
#' generated X will have p columns.
#' @param k_unclustered Integer or numeric; the number of features in X that
#' will have nonzero coefficients in the true model for y among those features 
#' not generated from the n_clusters latent variables (called "weak signal" 
#' features in the simulations from Faletto and Bien 2022). The coefficients on
#' these features will be determined by beta_unclustered.
#' @param cluster_size Integer or numeric; for each of the n_clusters latent
#' variables, X will contain cluster_size noisy proxies that are correlated with
#' the latent variable.
#' @param n_strong_cluster_vars Integer or numeric; among the cluster_size
#' proxies in each cluster, n_strong_cluster_vars will have a high covariance
#' (rho_high) with the latent variable and cluster_size - n_strong_cluster_vars
#' will have a low covariance (rho_low) with the latent variable.
#' @param n_clusters Integer or numeric; the number of latent variables to
#' generate, each of which will be associated with an observed cluster in X.
#' Must be at least 1. Default is 1.
#' @param sig_clusters Integer or numeric; the number of generated latent
#' features that will have nonzero coefficients in the true model for y (all of
#' them will have coefficient beta_latent). Must be less than or equal to
#' n_clusters. Default is 1.
#' @param rho_high Integer or numeric; the covariance of the "strong proxies" in
#' each cluster with the latent variable (and each other). Note that the
#' correlation between the "strong proxy" features in the cluster will be
#' rho_high/var. rho_high cannot equal 0 and must be at least as large as
#' rho_low. Default is 0.9.
#' @param rho_low Integer or numeric; the covariance of the "weak proxies" in
#' each cluster with the latent variable (and each other). Note that the
#' correlation between the "weak proxy" features in the cluster will be
#' rho_low/var. rho_low cannot equal 0 and must be no larger than rho_high.
#' Default is 0.5.
#' @param beta_latent Integer or numeric; the coefficient used for all
#' sig_clusters latent variables that have nonzero coefficients in the true
#' model for y. Can't equal 0. Default is 1.5.
#' @param beta_unclustered Integer or numeric; the maximum coefficient in the
#' model for y among the k_unclustered features in X not generated from the
#' latent variables. The coefficients of the features will be
#' beta_unclustered/sqrt(1:k_unclustered). Can't equal 0. Default is 1.
#' @param snr Integer or numeric; the signal-to-noise ratio of the response
#' y. If sigma_eps_sq is not specified, the variance of the noise in y will be
#' calculated using the formula sigma_eps_sq = sum(mu^2)/(n * snr). Only one of
#' snr and sigma_eps_sq must be specified. Default is NA.
#' @param sigma_eps_sq Integer or numeric; the variance on the noise added
#' to y. Only one of snr and sigma_eps_sq must be specified. Default is NA.
#' @return A list of the following elements. \item{X}{An n x p numeric matrix of
#' n observations from a p-dimensional multivariate normal distribution
#' generated using the specified parameters. The first n_clusters times
#' cluster_size features will be the clusters of features correlated with the
#' n_clusters latent variables. The next k_unclustered features will be the
#' "weak signal" features, and the remaining p - n_clusters*cluster_size -
#' k_unclustered features will be the unclustered noise features.} \item{y}{A
#' length n numeric vector; the response generated from X, the latent features
#' from Z, and the coefficient vector, along with additive noise.} \item{Z}{The
#' latent features; either a numeric vector (if n_clusters > 1) or a numeric
#' matrix (if n_clusters > 1). Note that (X, Z) is multivariate Gaussian.}
#' \item{mu}{A length `n` numeric vector; the expected response given X, Z, and
#' the true coefficient vector (equal to y minus the added noise).}
#' @author Gregory Faletto, Jacob Bien
checkGenClusteredDataWeightedInputs <- function(p, k_unclustered, cluster_size,
        n_strong_cluster_vars, n_clusters, sig_clusters, rho_high, rho_low,
        beta_latent, beta_unclustered, snr, sigma_eps_sq){

    stopifnot(is.numeric(sig_clusters) | is.integer(sig_clusters))
    stopifnot(sig_clusters <= n_clusters)
    stopifnot(sig_clusters >= 0)
    stopifnot(sig_clusters == round(sig_clusters))
    
    stopifnot(is.numeric(n_clusters) | is.integer(n_clusters))
    stopifnot(n_clusters == round(n_clusters))
    # TODO(gregfaletto): is it easy to remove the requirement that n_clusters is
    # at least 1 (so that it's possible to generate data with no latent 
    # features)? If so, should only check that cluster_size >= 1 if n_clusters
    # >= 1, and in makeCovarianceMatrix function only need block_size >= 1
    # rather than 2.
    stopifnot(n_clusters >= 1)

    stopifnot(cluster_size >= 2)

    stopifnot(is.integer(n_strong_cluster_vars) |
        is.numeric(n_strong_cluster_vars))
    stopifnot(!is.na(n_strong_cluster_vars))
    stopifnot(length(n_strong_cluster_vars) == 1)
    stopifnot(n_strong_cluster_vars == round(n_strong_cluster_vars))
    stopifnot(n_strong_cluster_vars >= 1)
    stopifnot(n_strong_cluster_vars < cluster_size)

    stopifnot(rho_high <= 1)
    stopifnot(rho_high > 0)
    stopifnot(rho_low > 0)
    stopifnot(rho_high >= rho_low)

    stopifnot(beta_latent != 0)
    stopifnot(beta_unclustered != 0)

    stopifnot(is.numeric(k_unclustered) | is.integer(k_unclustered))
    stopifnot(k_unclustered >= 2)
    stopifnot(k_unclustered == round(k_unclustered))
    
    stopifnot(p >= n_clusters*cluster_size + k_unclustered)

    # Same as make_sparse_blocked_linear_model_random, but ith coefficient
    # of weak signal features is beta_unclustered/sqrt(i) in order to have
    # a definitive ranking of weak signal features.
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }
    
    if(!is.na(snr)){
        stopifnot(snr > 0)
    }
    if(!is.na(sigma_eps_sq)){
        stopifnot(sigma_eps_sq > 0)
    }
}

#' Check inputs to genClusteredDataWeightedRandom
#'
#' @param n Integer or numeric; the number of observations to generate. (The
#' generated X and Z will have n rows, and the generated y and mu will have
#' length n.)
#' @param p Integer or numeric; the number of features to generate. The
#' generated X will have p columns.
#' @param k_unclustered Integer or numeric; the number of features in X that
#' will have nonzero coefficients in the true model for y among those features 
#' not generated from the n_clusters latent variables (called "weak signal" 
#' features in the simulations from Faletto and Bien 2022). The coefficients on
#' these features will be determined by beta_unclustered.
#' @param cluster_size Integer or numeric; for each of the n_clusters latent
#' variables, X will contain cluster_size noisy proxies that are correlated with
#' the latent variable.
#' @param n_clusters Integer or numeric; the number of latent variables to
#' generate, each of which will be associated with an observed cluster in X.
#' Must be at least 1. Default is 1.
#' @param sig_clusters Integer or numeric; the number of generated latent
#' features that will have nonzero coefficients in the true model for y (all of
#' them will have coefficient beta_latent). Must be less than or equal to
#' n_clusters. Default is 1.
#' @param rho_high Integer or numeric; the maximum correlation of the proxies
#' each cluster with each other. Default is 1.
#' @param rho_low Integer or numeric; the minimum correlation of the proxies in
#' each cluster with each other. rho_low cannot equal 0 and must be no larger
#' than rho_high. Default is 0.5.
#' @param beta_latent Integer or numeric; the coefficient used for all
#' sig_clusters latent variables that have nonzero coefficients in the true
#' model for y. Can't equal 0. Default is 1.5.
#' @param beta_unclustered Integer or numeric; the maximum coefficient in the
#' model for y among the k_unclustered features in X not generated from the
#' latent variables. The coefficients of the features will be
#' beta_unclustered/sqrt(1:k_unclustered). Can't equal 0. Default is 1.
#' @param snr Integer or numeric; the signal-to-noise ratio of the response
#' y. If sigma_eps_sq is not specified, the variance of the noise in y will be
#' calculated using the formula sigma_eps_sq = sum(mu^2)/(n * snr). Only one of
#' snr and sigma_eps_sq must be specified. Default is NA.
#' @param sigma_eps_sq Integer or numeric; the variance on the noise added
#' to y. Only one of snr and sigma_eps_sq must be specified. Default is NA.
#' @return A list of the following elements. \item{X}{An n x p numeric matrix of
#' n observations from a p-dimensional multivariate normal distribution
#' generated using the specified parameters. The first n_clusters times
#' cluster_size features will be the clusters of features correlated with the
#' n_clusters latent variables. The next k_unclustered features will be the
#' "weak signal" features, and the remaining p - n_clusters*cluster_size -
#' k_unclustered features will be the unclustered noise features.} \item{y}{A
#' length n numeric vector; the response generated from X, the latent features
#' from Z, and the coefficient vector, along with additive noise.} \item{Z}{The
#' latent features; either a numeric vector (if n_clusters > 1) or a numeric
#' matrix (if n_clusters > 1). Note that (X, Z) is multivariate Gaussian.}
#' \item{mu}{A length `n` numeric vector; the expected response given X, Z, and
#' the true coefficient vector (equal to y minus the added noise).}
#' @author Gregory Faletto, Jacob Bien
checkGenClusteredDataWeightedRandomInputs <- function(p, k_unclustered,
    cluster_size, n_clusters, sig_clusters, rho_high, rho_low, beta_latent,
    beta_unclustered, snr, sigma_eps_sq){

    stopifnot(is.numeric(sig_clusters) | is.integer(sig_clusters))
    stopifnot(sig_clusters <= n_clusters)
    stopifnot(sig_clusters >= 0)
    stopifnot(sig_clusters == round(sig_clusters))
    
    stopifnot(is.numeric(n_clusters) | is.integer(n_clusters))
    stopifnot(n_clusters == round(n_clusters))
    # TODO(gregfaletto): is it easy to remove the requirement that n_clusters is
    # at least 1 (so that it's possible to generate data with no latent 
    # features)? If so, should only check that cluster_size >= 1 if n_clusters
    # >= 1, and in makeCovarianceMatrix function only need block_size >= 1
    # rather than 2.
    stopifnot(n_clusters >= 1)

    stopifnot(cluster_size >= 2)

    stopifnot(rho_high <= 1)
    stopifnot(rho_high > 0)
    stopifnot(rho_low > 0)
    stopifnot(rho_high >= rho_low)

    stopifnot(beta_latent != 0)
    stopifnot(beta_unclustered != 0)

    stopifnot(is.numeric(k_unclustered) | is.integer(k_unclustered))
    stopifnot(k_unclustered >= 2)
    stopifnot(k_unclustered == round(k_unclustered))
    
    stopifnot(p >= n_clusters*cluster_size + k_unclustered)

    # Same as make_sparse_blocked_linear_model_random, but ith coefficient
    # of weak signal features is beta_unclustered/sqrt(i) in order to have
    # a definitive ranking of weak signal features.
    if(is.na(snr) & is.na(sigma_eps_sq)){
        stop("Must specify one of snr or sigma_eps_sq")
    }
    
    if(!is.na(snr)){
        stopifnot(snr > 0)
    }
    if(!is.na(sigma_eps_sq)){
        stopifnot(sigma_eps_sq > 0)
    }
}

#' Generate covariance matrix for simulated clustered data
#'
#' @param p Integer or numeric; the total number of features in the covariance
#' matrix to be created, including latent features, the associated noisy proxies
#' with each latent feature, and other (weak signal and noise) features.
#' @param n_blocks Integer or numeric; the number of latent variables in the
#' data, each of which is associated with an observed cluster in X. Must be at
#' least 1.
#' @param block_size Integer or numeric; for each of the n_blocks latent
#' variables, the covariance matrix will include the original latent feature
#' plus block_size - 1 noisy proxies that are correlated with the latent
#' variable.
#' @param n_strong_block_vars Integer or numeric; the number of proxies that
#' will have a high covariance (rho_high) with the corresponding latent 
#' features. (The remaining block_size - n_strong_block_vars cluster members
#' will have covariance rho_low with the latent feature.)
#' @param rho_high Integer or numeric; the covariance of the "strong proxies" in
#' each cluster with the latent variable (and each other). Note that the
#' correlation between the "strong proxy" features in the cluster will be
#' rho_high/var. rho_high cannot equal 0 and must be at least as large as
#' rho_low.
#' @param rho_low Integer or numeric; the covariance of the "weak proxies" in
#' each cluster with the latent variable (and each other). Note that the
#' correlation between the "weak proxy" features in the cluster will be
#' rho_low/var. rho_low cannot equal 0 and must be no larger than rho_high.
#' @param var Integer or numeric; the variance of all of the observed features
#' in X (both the proxies for the latent variables and the k_unclustered other
#' features). var cannot equal 0.
#' @return A `p` x `p` numeric matrix representing the covariance matrix for
#' the latent features, the associated proxies, and the remaining features. All
#' features not in a block will be independent from each other and the blocks
#' and have variance var.
#' @author Gregory Faletto, Jacob Bien
makeCovarianceMatrixWeighted <- function(p, nblocks, block_size,
    n_strong_block_vars, rho_high, rho_low, var) {
    # Check inputs

    stopifnot(nblocks >= 1)
    stopifnot(rho_high != 0)
    stopifnot(rho_low != 0)
    stopifnot(abs(rho_high) >= abs(rho_low))
    stopifnot(var != 0)
    stopifnot(abs(rho_high) <= abs(var))
    stopifnot(block_size >= 2)
    stopifnot(p >= nblocks*block_size)
    stopifnot(n_strong_block_vars <= block_size)

    # start with p x p identity matrix
    Sigma <- var*diag(p)

    # create matrix with nblocks rows, each containing a vector of
    # indices of highly correlated features
    block_feats <- matrix(seq(nblocks*block_size), nrow=nblocks, byrow=TRUE)

    stopifnot(length(unique(block_feats)) == length(block_feats))

    # add covariances of highly correlated features to sigma
    for(i in 1:nblocks){
        for(j in 1:(block_size - 1)){
            for(k in (j+1):block_size){
                feat_1 <- block_feats[i, j]
                feat_2 <- block_feats[i, k]
                Sigma[feat_1, feat_2] <- rho_low
                Sigma[feat_2, feat_1] <- rho_low
            }
        }

        for(j in 1:(n_strong_block_vars - 1)){
            for(k in (j+1):n_strong_block_vars){
                feat_1 <- block_feats[i, j]
                feat_2 <- block_feats[i, k]
                Sigma[feat_1, feat_2] <- rho_high
                Sigma[feat_2, feat_1] <- rho_high
            }
        }
    }

    stopifnot(is.numeric(Sigma))
    stopifnot(is.matrix(Sigma))
    stopifnot(nrow(Sigma) == p & ncol(Sigma) == p)
    stopifnot(all(Sigma == t(Sigma)))

    return(Sigma)
}

