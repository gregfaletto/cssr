#' Converts the provided design matrix to an appropriate format for either the
#' protolasso or the cluster representative lasso.
#'
#' Creates design matrix for glmnet by dealing with clusters (for
#' type="protolasso", discards all cluster members except prototype; for
#' type="clusterRepLasso", replaces all cluster members with a simple
#' average of all the cluster members).
#' @param x A numeric matrix; the provided matrix with n observations and p
#' features.
#' @param clusters A named list where each entry is an integer vector of indices
#' of features that are in a common cluster. (The length of list clusters should
#' be equal to the number of clusters.) All identified clusters should be
#' non-overlapping. All features should appear in exactly one cluster (any
#' unclustered features should be put in their own "cluster" of size 1).
#' @param type Character; "protolasso" for the protolasso or "clusterRepLasso"
#' for the cluster representative lasso.
#' @param prototypes Only required for type "protolasso". An integer vector
#' whose length is equal to the number of clusters. Entry i should be the
#' prototype for cluster i (the feature belonging to cluster i that is most
#' highly correlated with y; see Reid and Tibshirani 2016).
#' @return A numeric matrix; the design matrix as required for the protolasso or
#' cluster representative lasso, prepared for input to glmnet. The number of 
#' columns will be equal to length(clusters), and column j will be a feature
#' corresponding to clusters[[j]] (either the prototype from clusters[[j]] or
#' a simple average of the features in clusters[[j]]).
#' @author Gregory Faletto, Jacob Bien
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}.
getXglmnet <- function(x, clusters, type, prototypes=NA){
    
    # Check inputs
    checkGetXglmnetInputs(x, clusters, type, prototypes)

    n <- nrow(x)
    p <- ncol(x)

    # if(n_clusters > 0){
    for(i in 1:length(clusters)){
        cluster_i <- clusters[[i]]

        if(length(cluster_i) == 1){
            X_glmnet_i <- x[, cluster_i]
        } else{
            stopifnot(length(cluster_i) > 1)

            if(type == "protolasso"){
                prototype_ind_i <- which(prototypes %in% cluster_i)
                stopifnot(length(prototype_ind_i) == 1)
                prototype_i <- prototypes[prototype_ind_i]
                X_glmnet_i <- x[, prototype_i]
            } else {
                stopifnot(type == "clusterRepLasso")
                X_glmnet_i <- rowMeans(x[, cluster_i])
            }
        }
        
        stopifnot(length(X_glmnet_i) == n)
        
        if(i == 1){
            X_glmnet <- as.matrix(X_glmnet_i)
        } else{
            X_glmnet <- cbind(X_glmnet, X_glmnet_i)
        }
    }
    stopifnot(ncol(X_glmnet) == length(clusters))
    
    # if(length(non_cluster_feats) > 0){
    #     X_glmnet <- cbind(X_glmnet, x[, non_cluster_feats])
    # } else{
    #     X_glmnet <- X_glmnet
    # }
    # stopifnot(ncol(X_glmnet) == length(non_cluster_feats) + n_clusters)
    # } else{
    #     X_glmnet <- x
    #     if(type == "protolasso"){
    #         stopifnot(length(prototypes) == 0)
    #     }
    #     cluster_members <- integer()
    #     non_cluster_feats <- 1:p
    # }

    colnames(X_glmnet) <- character()

    # Check output
    stopifnot(is.matrix(X_glmnet))
    stopifnot(nrow(X_glmnet) == n)
    stopifnot(ncol(X_glmnet) <= p)
    stopifnot(ncol(X_glmnet) >= 1)
    
    return(X_glmnet)

}


#' Verifies the inputs for getXglmnet.
#'
#' @param x A numeric matrix; the provided matrix with n observations and p
#' features.
#' @param clusters A named list where each entry is an integer vector of indices
#' of features that are in a common cluster. (The length of list clusters should
#' be equal to the number of clusters.) All identified clusters should be
#' non-overlapping. All features should appear in exactly one cluster (any
#' unclustered features should be put in their own "cluster" of size 1).
#' @param type Character; "protolasso" for the protolasso or "clusterRepLasso"
#' for the cluster representative lasso.
#' @param prototypes Only required for type "protolasso". An integer vector
#' whose length is equal to the number of clusters. Entry i should be the
#' prototype for cluster i (the feature belonging to cluster i that is most
#' highly correlated with y; see Reid and Tibshirani 2016).
#' @author Gregory Faletto, Jacob Bien
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}.
checkGetXglmnetInputs <- function(x, clusters, type, prototypes){
    stopifnot(is.matrix(x))

    stopifnot(is.list(clusters))
    stopifnot(all(lengths(clusters) >= 1))

    # stopifnot(length(n_clusters) == 1)
    # stopifnot(is.integer(n_clusters) | is.numeric(n_clusters))
    # stopifnot(n_clusters >= 0)
    # if(type=="protolasso"){
    #     stopifnot(n_clusters == length(prototypes))
    # }
    # stopifnot(n_clusters == length(clusters))

    stopifnot(length(type) == 1)
    stopifnot(is.character(type))
    stopifnot(!is.na(type))
    stopifnot(type %in% c("protolasso", "clusterRepLasso"))

    if(type=="protolasso"){
        stopifnot(!is.na(prototypes))
        stopifnot(is.integer(prototypes))
        # if(is.list(clusters)){
            # stopifnot(length(prototypes) == length(clusters))
        # } else{
        #     stopifnot(length(prototypes) == 1)
        # }
        stopifnot(all(!is.na(prototypes)))
        stopifnot(length(prototypes) == length(unique(prototypes)))
        stopifnot(all(prototypes %in% 1:ncol(x)))
    }
    for(i in 1:length(clusters)){
        cluster_i <- clusters[[i]]
        if(length(cluster_i) > 1){
            stopifnot(sum(prototypes %in% cluster_i) == 1)
        }
    }
}

#' Extracts selected clusters and cluster prototypes from the glmnet lasso
#' output
#'
#' @param lasso_sets A list of integer vectors. Each vector represents a set of
#' features selected by the lasso for a given value of the penalty parameter
#' lambda.
#' @param clusters A named list where each entry is an integer vector of indices
#' of features that are in a common cluster. (The length of list clusters is
#' equal to the number of clusters.) All identified clusters must be
#' non-overlapping. All features appear in exactly one cluster (any unclustered
#' features must be in their own "cluster" of size 1).
#' @param prototypes An integer vector whose length must be equal to the number
#' of clusters. Entry i should be the index of the feature belonging to cluster
#' i that is most highly correlated with y (that is, the prototype for the
#' cluster, as in the protolasso; see Reid and Tibshirani 2016).
#' @param feat_names Character vector; the names of the features in X. (If the
#' X provided to protolasso or clusterRepLasso did not have feature names,
#' feat_names will be NA.)
#' @return A list containing the following items: \item{selected_sets}{A list of
#' integer vectors. Entry k of this list contains a selected set of size k
#' yielded by glmnet--each member of the set is the index of a single feature
#' from a cluster selected by either the protolasso or the cluster
#' representative lasso (the prototype from that cluster--the cluster member
#' most highly correlated with y). (If no set of size k was selected, entry k
#' will be NULL.)} \item{selected_clusts_list}{A list of lists; entry k of this
#' list is a list of length k of clusters (the clusters that were selected by
#' the cluster representative lasso). Again, if no set of size k was selected,
#' entry k will be NULL.}
#' @author Gregory Faletto, Jacob Bien
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}. \cr Bühlmann, P.,
#' Rütimann, P., van de Geer, S., & Zhang, C. H. (2013). Correlated variables in
#' regression: Clustering and sparse estimation.
#' \emph{Journal of Statistical Planning and Inference}, 143(11), 1835–1858.
#' \url{https://doi.org/10.1016/j.jspi.2013.05.019}.
getClusterSelsFromGlmnet <- function(lasso_sets, clusters, prototypes,
    feat_names){
    # checkGetClusterSelsFromGlmnetInput(clusters, prototypes, averaging)

    if(any(!is.na(feat_names))){
        stopifnot(all(!is.na(feat_names)))
    }

    # Largest selected set among all those in lasso_sets
    max_length <- max(vapply(lasso_sets, length, integer(1)))

    # Preparing lists to store 
    selected_sets <- list()
    # to_avg_list <- list()
    selected_clusts_list <- list()
    # weights_list <- list()
    
    for(j in 1:max_length){
        # Lasso selected set of size j
        lasso_sets_j <- lasso_sets[lapply(lasso_sets, length) == j]
        # Are there any lasso selected sets of size j? (If not, we will skip to
        # the next j, and slot j in the list will be empty.)
        if(length(lasso_sets_j) > 0){

            # Select the first set of size j
            lasso_set_j <- lasso_sets_j[[1]]
            stopifnot(length(lasso_set_j) == j)
            
            ret <- getSelectedSets(lasso_set=lasso_set_j, clusters=clusters,
                prototypes=prototypes, feat_names=feat_names)

            selected_sets[[j]] <- ret$selected_set
            selected_clusts_list[[j]] <- ret$selected_clusts_list

            rm(ret)
        }
    }

    stopifnot(length(selected_sets) <= max_length)
    stopifnot(length(selected_clusts_list) <= max_length)

    # if(averaging){
    return(list(selected_sets=selected_sets,
        selected_clusts_list=selected_clusts_list))
    # }
    # else{
    #     return(selected_sets)
    # }
}

#' Converts a selected set from X_glmnet to selected sets and selected clusters
#' from the original feature space of X.
#'
#' @param lasso_set A vector containing the indices of selected cluster
#' representatives or prototypes.
#' @param clusters A named list where each entry is an integer vector of indices
#' of features that are in a common cluster. (The length of list clusters is
#' equal to the number of clusters.) All identified clusters must be
#' non-overlapping. All features appear in exactly one cluster (any unclustered
#' features must be in their own "cluster" of size 1).
#' @param prototypes An integer vector whose length must be equal to the number
#' of clusters. Entry i should be the index of the feature belonging to cluster
#' i that is most highly correlated with y (that is, the prototype for the
#' cluster, as in the protolasso).
#' @param feat_names Character vector; the names of the features in X.
#' @return A list containing two items: \item{selected_set}{An integer vector
#' with length equal to lasso_set containing a set of selected features in the
#' original X matrix. (Selections in lasso_set corresponding to a cluster will
#' be replaced by the cluster's prototype from X.)}
#' \item{selected_clusts_list}{A named list of integer vectors with length equal
#' to selected_set. selected_clusts_list[[k]] will be an integer vector
#' containing the indices of the features in X that are in the cluster
#' containing prototype selected_set[k].}
#' @author Gregory Faletto, Jacob Bien
getSelectedSets <- function(lasso_set, clusters, prototypes, feat_names){
    
    model_size <- length(lasso_set)
    stopifnot(model_size > 0)

    stopifnot(length(unique(lasso_set)) == model_size)
    stopifnot(all(lasso_set <= length(clusters)))

    selected_set <- integer()
    selected_clusts_list <- list()
    # Recover features from original feature space
    for(k in 1:model_size){
        selected_cluster_k <- clusters[[lasso_set[k]]]
        stopifnot(is.integer(selected_cluster_k))
        selected_clusts_list[[k]] <- selected_cluster_k

        if(length(selected_cluster_k) == 1){
            stopifnot(!(selected_cluster_k %in% selected_set))
            selected_set <- c(selected_set, selected_cluster_k)
        } else{
            sel_prototype <- which(prototypes %in% selected_cluster_k)
            stopifnot(length(sel_prototype) == 1)
            stopifnot(!(prototypes[sel_prototype] %in% selected_set))
            selected_set <- c(selected_set, prototypes[sel_prototype])
        }
    }


    # # : deal with
    # # prototype and non-prototype features separately
    # proto_inds <- lasso_set %in% prototypes
    # lasso_set_protos <- lasso_set[proto_inds]
    # lasso_set_non_protos <- lasso_set[!proto_inds]

    # # Recover indices of non-prototype features

    # if(length(lasso_set_non_protos) > 0){
    #     stopifnot(length(lasso_set_non_protos) <=
    #         length(non_cluster_feats))
    #     stopifnot(all(lasso_set_non_protos - n_clusters %in%
    #         1:length(non_cluster_feats)))

    #     recovered_feats <- non_cluster_feats[lasso_set_non_protos -
    #         n_clusters]

    #     stopifnot(length(lasso_set_non_protos) ==
    #         length(recovered_feats))
    #     stopifnot(all.equal(!proto_inds, lasso_set %in%
    #         lasso_set_non_protos))
    #     stopifnot(length(recovered_feats) ==
    #         length(unique(recovered_feats)))

    #     lasso_set_non_protos <- recovered_feats

    #     stopifnot(length(lasso_set_non_protos) == sum(!proto_inds))
    # }

    # stopifnot(length(unique(lasso_set_protos)) == length(lasso_set_protos))
    # stopifnot(length(lasso_set_protos) == sum(proto_inds))

    # proto_selected <- length(lasso_set_protos) > 0

    # stopifnot(all(lasso_set_protos %in% 1:n_clusters))

    # if(proto_selected){
    #     lasso_set_protos <- prototypes[lasso_set_protos]
    #     stopifnot(length(lasso_set_protos) == sum(proto_inds))
    # }

    # stopifnot(length(lasso_set) == length(lasso_set_protos) +
    #     length(lasso_set_non_protos))

    # lasso_set[proto_inds] <- lasso_set_protos
    # lasso_set[!proto_inds] <- lasso_set_non_protos

    # stopifnot(length(unique(lasso_set)) == length(lasso_set))
    
    # if(var_names_provided){
    #     stopifnot(max(lasso_set) <= length(var_names))
    #     selected_sets <- var_names[lasso_set]
    # } else{
    #     selected_sets <- lasso_set
    # }

    stopifnot(length(selected_set) == model_size)
    stopifnot(length(unique(selected_set)) == model_size)
    
    if(any(!is.na(feat_names))){
        names(selected_set) <- feat_names[selected_set]
    }

    # if(averaging){

    # to_avg_list[[model_size]] <- logical(model_size)

    # Finally, identify selected clusters
    
    # weights_list[[model_size]] <- list()

    stopifnot(length(selected_clusts_list) == model_size)
    all_feats <- unlist(selected_clusts_list)
    stopifnot(length(all_feats) == length(unique(all_feats)))

    return(list(selected_set=selected_set,
        selected_clusts_list=selected_clusts_list))


}

# getCrlOutputs <- function(to_avg_list, selected_clusts_list, weights_list,
#         j, lasso_set_j, n_clusters, clusters){

#     ind_j <- matrix(FALSE, nrow=length(lasso_set_j),
#         ncol=n_clusters)
#     stopifnot(is.list(clusters))
#     # if(is.list(clusters)){
        
#     for(k in 1:n_clusters){
#         n_cluster_members_k <- length(clusters[[k]])

#         stopifnot(n_cluster_members_k >= 1)

#         ind_j[, k] <- lasso_set_j == prototypes[k]

#         stopifnot(sum(ind_j[, k]) %in% c(0, 1))

#         if(sum(ind_j[, k]) == 1){
#             to_avg_list[[j]][which(ind_j[, k])] <- TRUE

#             # selected_sets[[j]][lasso_set_j == 1] <- cluster_prototype
#             if(var_names_provided){
#                 selected_clusts_list[[j]][[which(ind_j[, k])]] <-
#                     var_names[clusters[[k]]]
#             } else{
#                 selected_clusts_list[[j]][[which(ind_j[, k])]] <-
#                     clusters[[k]]
#             }

#             weights_list[[j]][[which(ind_j[, k])]] <-
#                 rep(1/n_cluster_members_k, n_cluster_members_k)
#         }
#     }
#     # } else{
#     #     stopifnot(n_clusters <= 1)
#     #     n_cluster_members_k <- length(clusters)
#     #     stopifnot(n_cluster_members_k != 1)

#     #     ind_j[, 1] <- lasso_set_j == prototypes

#     #     stopifnot(sum(ind_j[, 1]) %in% c(0, 1))

#     #     if(sum(ind_j[, 1]) == 1){
#     #         to_avg_list[[j]][which(ind_j[, 1])] <- TRUE

#     #         # selected_sets[[j]][lasso_set_j == 1] <- cluster_prototype
#     #         if(var_names_provided){
#     #             selected_clusts_list[[j]][[which(ind_j[, 1])]] <-
#     #                 var_names[clusters]
#     #         } else{
#     #             selected_clusts_list[[j]][[which(ind_j[, 1])]] <-
#     #                 clusters
#     #         }


            
#     #         weights_list[[j]][[which(ind_j[, 1])]] <-
#     #             rep(1/n_cluster_members_k, n_cluster_members_k)

#     #         stopifnot(length(weights_list[[j]][[which(ind_j[, 1])]]) ==
#     #             length(selected_clusts_list[[j]][[which(ind_j[, 1])]]))
#     #     }
#     # }
        
#     stopifnot(sum(ind_j) != 0)
#     stopifnot(all(rowSums(ind_j) <= 1))

#     return(list(to_avg_list=to_avg_list,
#         selected_clusts_list=selected_clusts_list, weights_list=weights_list))
# }

# checkGetClusterSelsFromGlmnetInput <- function(clusters, prototypes, averaging){

#     stopifnot(is.list(clusters))
#     stopifnot(all(lengths(clusters) >= 1))
    
#     stopifnot(is.integer(prototypes))
#     stopifnot(all(!is.na(prototypes)))
#     stopifnot(length(prototypes) <= length(clusters))
#     stopifnot(length(prototypes) == length(unique(prototypes)))

#     # stopifnot(length(var_names_provided) == 1)
#     # stopifnot(is.logical(var_names_provided))
#     # if(var_names_provided){
#     #     stopifnot(is.character(var_names))
#     #     if(any(is.na(var_names))){
#     #         stop("must provide var_names (with no NAs) if var_names_provided=TRUE")
#     #     }
#     # }

#     stopifnot(length(averaging) == 1)
#     stopifnot(is.logical(averaging))

#     # n_clusters <- sum(lengths(clusters) > 1)
#     # stopifnot(length(prototypes) == n_clusters)

#     # return(n_clusters)

# }



#' Select features via the protolasso (Reid and Tibshirani 2016)
#'
#' @param X An n x p numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' p >= 2 features/predictors
#' @param y The response; A length n numeric (or integer) real-valued vector.
#' @param clusters A list of integer vectors; each vector should contain the 
#' indices of a cluster of features (a subset of 1:p). (If there is only one
#' cluster, clusters can either be a list of length 1 or an integer vector.)
#' All of the provided clusters must be non-overlapping. Every feature not
#' appearing in any cluster will be assumed to be unclustered (that is, they
#' will be treated as if they are in a "cluster" containing only themselves).
#' Default is list() (so no clusters are specified).
#' @param nlambda Integer; the number of lambda values to use in the lasso fit
#' for the protolasso. Default is 100 (following the default for glmnet). For
#' now, nlambda must be at least 2 (using a single lambda is not supported).
#' @return A list with three elements. \item{selected_sets}{A list of integer
#' vectors. Entry k of this list contains a selected set (an integer vector) of
#' size k yielded by the protolasso (If no set of size k was selected, entry k
#' will be empty.)} \item{selected_clusts_list}{A list; each element of the list
#' is a named list of selected clusters. (That is, if a selected set of size k
#' was yielded by the protolasso, then selected_clusts_list[[k]] is a named
#' list of length k, where each member of the list is an integer vector
#' of cluster members. In particular, selected_clusts_lists[[k]][[j]] will be
#' the cluster that contains feature selected_sets[[k]][j].)} \item{beta}{The
#' beta output from glmnet when the lasso was estimated on a matrix of
#' prototypes. (See documentation for the function glmnet from the glmnet
#' package for details.)}
#' @author Gregory Faletto, Jacob Bien
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}.
#' @export
protolasso <- function(X, y, clusters, nlambda=100){

    # Handle and format inputs; get cluster prototypes
    ret <- processClusterLassoInputs(X, y, clusters, nlambda)

    x <- ret$x
    clusters <- ret$clusters
    prototypes <- ret$prototypes
    feat_names <- ret$var_names

    rm(ret)

    # Format the design matrix for glmnet according to the protolasso procedure
    X_glmnet <- getXglmnet(x, clusters, type="protolasso",
        prototypes=prototypes)

    #  <- getXglmnet_results$X_glmnet
    # cluster_members <- getXglmnet_results$cluster_members
    # non_cluster_feats <- getXglmnet_results$non_cluster_feats

    # rm(getXglmnet_results)

    # Estimate the lasso on the cluster prototypes
    fit <- glmnet::glmnet(x=X_glmnet, y=y, family="gaussian", nlambda=nlambda)
    lasso_sets <- unique(glmnet::predict.glmnet(fit, type="nonzero"))

    # Finally, obtain a tidy list of selected sets--one for each model size
    cluster_sel_results <- getClusterSelsFromGlmnet(lasso_sets, clusters,
        prototypes, feat_names)

    return(list(selected_sets=cluster_sel_results$selected_sets,
        selected_clusts_list=cluster_sel_results$selected_clusts_list,
        beta=fit$beta))
}

#' Check the inputs to protolasso and clusterRepLasso, format clusters, and
#' identify prototypes for each cluster
#'
#' @param X An n x p numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' p >= 2 features/predictors
#' @param y The response; A length n numeric (or integer) real-valued vector.
#' @param clusters A list of integer vectors; each vector should contain the 
#' indices of a cluster of features (a subset of 1:p). (If there is only one
#' cluster, clusters can either be a list of length 1 or an integer vector.)
#' All of the provided clusters must be non-overlapping. Every feature not
#' appearing in any cluster will be assumed to be unclustered (that is, they
#' will be treated as if they are in a "cluster" containing only themselves).
#' Default is list() (so no clusters are specified).
#' @param nlambda Integer; the number of lambda values to use in the lasso fit
#' for the protolasso. Default is 100 (following the default for glmnet). For
#' now, nlambda must be at least 2 (using a single lambda is not supported).
#' @return A list with four elements. \item{x}{The provided X, converted to a
#' matrix if it was provided as a data.frame, and with column names removed.}
#' \item{clusters}{A named list where each entry is an integer vector of indices
#' of features that are in a common cluster. (The length of list clusters is
#' equal to the number of clusters.) All identified clusters are
#' non-overlapping. All features appear in exactly one cluster (any unclustered
#' features will be put in their own "cluster" of size 1).}
#' \item{prototypes}{An integer vector whose length is equal to the number of
#' clusters. Entry i is the index of the feature belonging to cluster i that is
#' most highly correlated with y (that is, the prototype for the cluster, as in
#' the protolasso; see Reid and Tibshirani 2016).} \item{var_names}{If the
#' provided X matrix had column names, the names of the featurrs in the provided
#' X matrix. If no names were provided, feat_names will be NA.}
#' @author Gregory Faletto, Jacob Bien
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}.
processClusterLassoInputs <- function(X, y, clusters, nlambda){

    stopifnot(is.matrix(X) | is.data.frame(X))

    # Check if x is a matrix; if it's a data.frame, convert to matrix.
    if(is.data.frame(X)){
        X <- stats::model.matrix(~ ., X)
        X <- X[, colnames(X) != "(Intercept)"]
    }

    stopifnot(is.matrix(X))
    stopifnot(all(!is.na(X)))

    feat_names <- as.character(NA)
    if(!is.null(colnames(X))){
        feat_names <- colnames(X)
        if(any(is.na(feat_names))){
            stop("Some features in provided X matrix had valid names and some had NA names; please neither name all features in X or remove the names altogether.")
        }
    }

    n <- nrow(X)

    colnames(X) <- character()

    stopifnot(is.numeric(y) | is.integer(y))
    stopifnot(n == length(y))
    stopifnot(all(!is.na(y)))

    # Check clusters argument
    clusters <- checkCssClustersInput(clusters)

    # Format clusters into a list where all features are in exactly one
    # cluster (any unclustered features are put in their own "cluster" of size
    # 1).
    clust_names <- as.character(NA)
    if(!is.null(names(clusters)) & is.list(clusters)){
        clust_names <- names(clusters)
    }

    cluster_results <- formatClusters(clusters, p=ncol(X),
        clust_names=clust_names, get_prototypes=TRUE, x=X, y=y)

    clusters <- cluster_results$clusters
    prototypes <- cluster_results$prototypes

    rm(cluster_results)

    stopifnot(length(clusters) == length(prototypes))

    stopifnot(is.numeric(nlambda) | is.integer(nlambda))
    stopifnot(length(nlambda) == 1)
    stopifnot(!is.na(nlambda))
    stopifnot(nlambda >= 2)
    stopifnot(nlambda == round(nlambda))

    return(list(x=X, clusters=clusters, prototypes=prototypes,
        var_names=feat_names))
}

#' Select features via the cluster representative lasso (Bühlmann et. al. 2013)
#'
#' @param X An n x p numeric matrix (preferably) or a data.frame (which will
#' be coerced internally to a matrix by the function model.matrix) containing
#' p >= 2 features/predictors
#' @param y The response; A length n numeric (or integer) real-valued vector.
#' @param clusters A list of integer vectors; each vector should contain the 
#' indices of a cluster of features (a subset of 1:p). (If there is only one
#' cluster, clusters can either be a list of length 1 or an integer vector.)
#' All of the provided clusters must be non-overlapping. Every feature not
#' appearing in any cluster will be assumed to be unclustered (that is, they
#' will be treated as if they are in a "cluster" containing only themselves).
#' Default is list() (so no clusters are specified).
#' @param nlambda Integer; the number of lambda values to use in the lasso fit
#' for the cluster representative lasso. Default is 100 (following the default
#' for glmnet). For now, nlambda must be at least 2 (using a single lambda is
#' not supported).
#' @return A list with three elements. \item{selected_sets}{A list of integer
#' vectors. Entry k of this list contains a selected set (an integer vector) of
#' size k yielded by the lasso--each member of the set is the index of a single
#' feature from a cluster selected by the cluster representative lasso (the
#' prototype from that cluster--the cluster member most highly correlated with
#' y). (If no set of size k was selected, entry k will be empty.)}
#' \item{selected_clusts_list}{A list; each element of the list is a named list
#' of selected clusters. (That is, if a selected set of size k was yielded by
#' the cluster representative lasso, then selected_clusts_list[[k]] is a named
#' list of length k, where each member of the list is an integer vector
#' of cluster members. Note that selected_clusts_lists[[k]][[j]] will be the
#' cluster that contains feature selected_sets[[k]][j].)} \item{beta}{The beta
#' output from glmnet when the lasso was estimated on a matrix of prototypes.
#' (See documentation for the function glmnet from the glmnet package for
#' details.)}
#' @references Bühlmann, P., Rütimann, P., van de Geer, S., & Zhang, C. H.
#' (2013). Correlated variables in regression: Clustering and sparse estimation.
#' \emph{Journal of Statistical Planning and Inference}, 143(11), 1835–1858.
#' \url{https://doi.org/10.1016/j.jspi.2013.05.019}. \cr Jerome Friedman, Trevor
#' Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear
#' Models via Coordinate Descent. \emph{Journal of Statistical Software}, 33(1)
#' ' 1-22. URL \url{https://www.jstatsoft.org/v33/i01/}.
clusterRepLasso <- function(X, y, clusters=list(), nlambda=100){

    # Handle and format inputs; get cluster prototypes
    ret <- processClusterLassoInputs(X, y, clusters, nlambda)

    x <- ret$x
    clusters <- ret$clusters
    prototypes <- ret$prototypes
    feat_names <- ret$var_names

    rm(ret)

    # Format the design matrix for glmnet according to the cluster
    # representative lasso procedure
    X_glmnet <- getXglmnet(x, clusters, type="clusterRepLasso")

    #  <- getXglmnet_results$X_glmnet
    # cluster_members <- getXglmnet_results$cluster_members
    # non_cluster_feats <- getXglmnet_results$non_cluster_feats

    # rm(getXglmnet_results)

    # Estimate the lasso on the cluster representatives
    fit <- glmnet::glmnet(x=X_glmnet, y=y, family="gaussian", nlambda=nlambda)
    lasso_sets <- unique(glmnet::predict.glmnet(fit, type="nonzero"))

    # Finally, extract the desired information from the lasso fit--all the
    # sets of selected clusters (one for each observed model size), and
    # corresponding sets of selected features
    cluster_sel_results <- getClusterSelsFromGlmnet(lasso_sets, clusters,
        prototypes, feat_names)

    return(list(selected_sets=cluster_sel_results$selected_sets,
        selected_clusts_list=cluster_sel_results$selected_clusts_list,
        beta=fit$beta))
}