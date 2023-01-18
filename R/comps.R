getXglmnet <- function(x, clusters, type, prototypes=NA){
    # Creates design matrix for glmnet by dealing with clusters (for
    # type="protolasso", discards all cluster members except prototype; for
    # type="clusterRepLasso", replaces all cluster members with a simple
    # average of all the cluster members).

    # Check inputs
    checkGetXglmnetInputs(x, clusters, type, prototypes)

    n <- nrow(x)
    p <- ncol(x)

    # if(n_clusters > 0){
    for(i in 1:length(clusters)){
        cluster_i <- clusters[[i]]
        
        if(type == "protolasso"){
            if(length(cluster_i) > 1){
                prototype_ind_i <- which(prototypes %in% cluster_i)
                stopifnot(length(prototype_ind_i) == 1)
                prototype_i <- prototypes[prototype_ind_i]
            } else{
                prototype_i <- cluster_i
            }
            
            X_glmnet_i <- x[, prototype_i]
        } else {
            stopifnot(type == "clusterRepLasso")
            if(length(cluster_i) > 1){
                X_glmnet_i <- rowMeans(x[, cluster_i])
            } else{
                X_glmnet_i <- x[, cluster_i]
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
    checkGetXglmnetOutput(n, p, X_glmnet)
    
    return(X_glmnet)

}

checkGetXglmnetOutput <- function(n, p, X_glmnet){
    stopifnot(is.matrix(X_glmnet))
    stopifnot(nrow(X_glmnet) == n)
    stopifnot(ncol(X_glmnet) <= p)
    stopifnot(ncol(X_glmnet) >= 1)

    # stopifnot(is.numeric(cluster_members) | is.integer(cluster_members))
    # stopifnot(length(cluster_members) <= p)
    # stopifnot(length(cluster_members) >= 0)
    # if(length(cluster_members) >= 1){
    #     stopifnot(length(cluster_members) == length(unique(cluster_members)))
    #     stopifnot(all(round(cluster_members) == cluster_members))
    #     stopifnot(all(cluster_members) %in% 1:p)
    # }

    # stopifnot(is.numeric(non_cluster_feats) | is.integer(non_cluster_feats))
    # stopifnot(length(non_cluster_feats) <= p)
    # stopifnot(length(non_cluster_feats) >= 0)
    # if(length(non_cluster_feats) >= 1){
    #     stopifnot(length(non_cluster_feats) == length(unique(non_cluster_feats)))
    #     stopifnot(all(round(non_cluster_feats) == non_cluster_feats))
    #     stopifnot(all(non_cluster_feats) %in% 1:p)
    # }
    
    # stopifnot(length(intersect(cluster_members, non_cluster_feats)) == 0)
    # stopifnot(length(cluster_members) + length(non_cluster_feats) == p)

}

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

# TODO(gregfaletto): figure out how to eliminate need for prototypes argument in
# getClusterSelsFromGlmnet when called by clusterRepLasso (shouldn't be
# necessary for anything)

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
#' @param averaging Logical; if TRUE, then the features within each cluster will
#' be averaged (as in the cluster representative lasso); if FALSE, then only
#' the prototype from each cluster will be selected (as in the protolasso).
#' @return If averaging is FALSE, a list of integer vectors. Entry k of this
#' list contains a selected set of size k yielded by glmnet--each member of the
#' set is the index of a single feature from a cluster selected by either the
#' protolasso or the cluster representative lasso (the prototype from that
#' cluster--the cluster member most highly correlated with y). (If no set of
#' size k was selected, entry k will be empty.) If averaging is TRUE, a list
#' containing the following itemsis returned: \item{selected_sets}{The same list
#' described above that is returned if averaging is FALSE.}
#' \item{selected_clusts_list}{A list of lists; entry k of this list is a list
#' of length k of clusters (the clusters that were selected by the cluster
#' representative lasso).}
#' @author Gregory Faletto, Jacob Bien
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}. \cr Bühlmann, P.,
#' Rütimann, P., van de Geer, S., & Zhang, C. H. (2013). Correlated variables in
#' regression: Clustering and sparse estimation.
#' \emph{Journal of Statistical Planning and Inference}, 143(11), 1835–1858.
#' \url{https://doi.org/10.1016/j.jspi.2013.05.019}.
getClusterSelsFromGlmnet <- function(lasso_sets, clusters, prototypes,
    averaging){

    # Check inputs
    checkGetClusterSelsFromGlmnetInput(clusters, prototypes, averaging)

    # Largest selected set among all those in lasso_sets
    max_length <- max(vapply(lasso_sets, length, integer(1)))

    # Preparing lists to store 
    selected_sets <- list()
    if(averaging){
        # to_avg_list <- list()
        selected_clusts_list <- list()
        # weights_list <- list()
    }
    
    for(j in 1:max_length){
        # Lasso selected set of size j
        lasso_sets_j <- lasso_sets[lapply(lasso_sets, length) == j]
        # Are there any lasso selected sets of size j? (If not, we will skip to
        # the next j, and slot j in the list will be empty.)
        if(length(lasso_sets_j) > 0){
            # Select the first set of size j
            lasso_set_j <- lasso_sets_j[[1]]
            stopifnot(length(unique(lasso_set_j)) == j)
            stopifnot(length(lasso_set_j) == j)
            stopifnot(all(lasso_set_j <= length(clusters)))


            selected_set_j <- integer()
            # Recover features from original feature space
            for(k in 1:j){
                selected_cluster_k <- clusters[[lasso_set_j[k]]]
                if(length(selected_cluster_k) == 1){
                    stopifnot(!(selected_cluster_k %in% selected_set_j))
                    selected_set_j <- c(selected_set_j, selected_cluster_k)
                } else{
                    sel_prototype <- which(prototypes %in% selected_cluster_k)
                    stopifnot(length(sel_prototype) == 1)
                    stopifnot(!(prototypes[sel_prototype] %in% selected_set_j))
                    selected_set_j <- c(selected_set_j,
                        prototypes[sel_prototype])
                }
                
            }


            # # : deal with
            # # prototype and non-prototype features separately
            # proto_inds <- lasso_set_j %in% prototypes
            # lasso_set_j_protos <- lasso_set_j[proto_inds]
            # lasso_set_j_non_protos <- lasso_set_j[!proto_inds]

            # # Recover indices of non-prototype features

            # if(length(lasso_set_j_non_protos) > 0){
            #     stopifnot(length(lasso_set_j_non_protos) <=
            #         length(non_cluster_feats))
            #     stopifnot(all(lasso_set_j_non_protos - n_clusters %in%
            #         1:length(non_cluster_feats)))

            #     recovered_feats <- non_cluster_feats[lasso_set_j_non_protos -
            #         n_clusters]

            #     stopifnot(length(lasso_set_j_non_protos) ==
            #         length(recovered_feats))
            #     stopifnot(all.equal(!proto_inds, lasso_set_j %in%
            #         lasso_set_j_non_protos))
            #     stopifnot(length(recovered_feats) ==
            #         length(unique(recovered_feats)))

            #     lasso_set_j_non_protos <- recovered_feats

            #     stopifnot(length(lasso_set_j_non_protos) == sum(!proto_inds))
            # }

            # stopifnot(length(unique(lasso_set_j_protos)) == length(lasso_set_j_protos))
            # stopifnot(length(lasso_set_j_protos) == sum(proto_inds))

            # proto_selected_j <- length(lasso_set_j_protos) > 0

            # stopifnot(all(lasso_set_j_protos %in% 1:n_clusters))

            # if(proto_selected_j){
            #     lasso_set_j_protos <- prototypes[lasso_set_j_protos]
            #     stopifnot(length(lasso_set_j_protos) == sum(proto_inds))
            # }

            # stopifnot(length(lasso_set_j) == length(lasso_set_j_protos) +
            #     length(lasso_set_j_non_protos))

            # lasso_set_j[proto_inds] <- lasso_set_j_protos
            # lasso_set_j[!proto_inds] <- lasso_set_j_non_protos

            # stopifnot(length(unique(lasso_set_j)) == length(lasso_set_j))
            
            # if(var_names_provided){
            #     stopifnot(max(lasso_set_j) <= length(var_names))
            #     selected_sets[[j]] <- var_names[lasso_set_j]
            # } else{
            #     selected_sets[[j]] <- lasso_set_j
            # }

            stopifnot(length(selected_set_j) == j)
            stopifnot(length(unique(selected_set_j)) == j)
            selected_sets[[j]] <- selected_set_j

            if(averaging){

                # to_avg_list[[j]] <- logical(j)
                selected_clusts_list[[j]] <- list()
                # weights_list[[j]] <- list()

                for(k in 1:j){
                    selected_cluster_k <- lasso_set_j[k]
                    cluster_k <- clusters[[selected_cluster_k]]
                    stopifnot(is.integer(cluster_k))
                    selected_clusts_list[[j]][[k]] <- cluster_k
                }

                stopifnot(length(selected_clusts_list[[j]]) == j)
                all_feats_j <- unlist(selected_clusts_list[[j]])
                stopifnot(length(all_feats_j) == length(unique(all_feats_j)))

                # if(proto_selected_j){
                #     crl_results <- getCrlOutputs(to_avg_list,
                #         selected_clusts_list, weights_list, j, lasso_set_j,
                #         n_clusters, clusters)

                #     # to_avg_list <- crl_results$to_avg_list
                #     selected_clusts_list <- crl_results$selected_clusts_list
                #     # weights_list <- crl_results$weights_list

                #     rm(crl_results)
                # } 
            }
        }
    }

    stopifnot(length(selected_sets) == max_length)
    stopifnot(length(selected_clusts_list) == max_length)

    if(averaging){
        return(list(selected_sets=selected_sets,
            selected_clusts_list=selected_clusts_list))
    }
    else{
        return(selected_sets)
    }
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

checkGetClusterSelsFromGlmnetInput <- function(clusters, prototypes, averaging){

    stopifnot(is.list(clusters))
    stopifnot(all(lengths(clusters) >= 1))
    
    stopifnot(is.integer(prototypes))
    stopifnot(all(!is.na(prototypes)))
    stopifnot(length(prototypes) <= length(clusters))
    stopifnot(length(prototypes) == length(unique(prototypes)))

    # stopifnot(length(var_names_provided) == 1)
    # stopifnot(is.logical(var_names_provided))
    # if(var_names_provided){
    #     stopifnot(is.character(var_names))
    #     if(any(is.na(var_names))){
    #         stop("must provide var_names (with no NAs) if var_names_provided=TRUE")
    #     }
    # }

    stopifnot(length(averaging) == 1)
    stopifnot(is.logical(averaging))

    # n_clusters <- sum(lengths(clusters) > 1)
    # stopifnot(length(prototypes) == n_clusters)

    # return(n_clusters)

}




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
#' @return A list with two elements. \item{selected_sets}{A list of integer
#' vectors. Entry k of this list contains a selected
#' set of size k yielded by the lasso--each member of the
#' set is the index of a single feature from a cluster selected by either the
#' protolasso or the cluster representative lasso (the prototype from that
#' cluster--the cluster member most highly correlated with y). (If no set of
#' size k was selected, entry k will be empty.)} \item{beta}{The beta output
#' from glmnet when the lasso was estimated on a matrix of prototypes. (See
#' documentation for the function glmnet from the glmnet package for details.)}
#' @references Reid, S., & Tibshirani, R. (2016). Sparse regression and marginal
#' testing using cluster prototypes. \emph{Biostatistics}, 17(2), 364–376.
#' \url{https://doi.org/10.1093/biostatistics/kxv049}.
protolasso <- function(x, y, clusters, var_names=NA, nlambda=4000){

    # Handle and format inputs
    ret <- processClusterLassoInputs(x, y, var_names, clusters)

    x <- ret$x
    var_names_provided <- ret$var_names_provided
    clusters <- ret$clusters
    prototypes <- ret$prototypes
    n_clusters <- ret$n_clusters

    rm(ret)

    n <- nrow(x)

    # Format the design matrix for glmnet according to the protolasso procedure
    X_glmnet <- getXglmnet(x, clusters, type="protolasso",
        prototypes=prototypes)

    #  <- getXglmnet_results$X_glmnet
    # cluster_members <- getXglmnet_results$cluster_members
    # non_cluster_feats <- getXglmnet_results$non_cluster_feats

    # rm(getXglmnet_results)

    stopifnot(nrow(X_glmnet) == n)

    # Estimate the lasso on the cluster prototypes
    fit <- glmnet::glmnet(x=X_glmnet, y=y, family="gaussian", nlambda=nlambda)
    lasso_sets <- unique(glmnet::predict.glmnet(fit, type="nonzero"))

    # Finally, obtain a tidy list of selected sets--one for each model size
    selected_sets <- getClusterSelsFromGlmnet(lasso_sets, clusters,
        prototypes, averaging=FALSE)

    return(list(selected_sets=selected_sets, beta=fit$beta))

}

processClusterLassoInputs <- function(x, y, var_names, clusters){
    if(is.data.frame(x)){
        x <- stats::model.matrix(~ ., x)
        x <- x[, colnames(x) != "(Intercept)"]
    }

    colnames(x) <- character()

    stopifnot(is.matrix(x))
    stopifnot(is.numeric(y) | is.integer(y))
    stopifnot(nrow(x) == length(y))

    var_names_provided <- FALSE
    # If var_names is provided, convert avg_feats entries to character vectors
    if(any(!is.na(var_names))){
        stopifnot(all(!is.na(var_names)))
        var_names_provided <- TRUE
        stopifnot(length(var_names) == ncol(x))
        stopifnot(is.character(var_names))
    }

    # Check clusters argument
    clusters <- checkCssClustersInput(clusters)

    # Format clusters into a list where all features are in exactly one
    # cluster (any unclustered features are put in their own "cluster" of size
    # 1).
    clust_names <- as.character(NA)
    if(!is.null(names(clusters)) & is.list(clusters)){
        clust_names <- names(clusters)
    }

    cluster_results <- formatClusters(clusters, p=ncol(x),
        clust_names=clust_names, get_prototypes=TRUE, x=x, y=y)

    clusters <- cluster_results$clusters
    # multiple <- cluster_results$multiple
    prototypes <- cluster_results$prototypes

    rm(cluster_results)

    stopifnot(is.list(clusters))
    n_clusters <- sum(lengths(clusters) > 1)
    stopifnot(n_clusters == length(prototypes))

    return(list(x=x, var_names_provided=var_names_provided, clusters=clusters,
        prototypes=prototypes, n_clusters=n_clusters))
}

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
#' @return A list with three elements. \item{selected_sets}{A list of integer
#' vectors. Entry k of this list contains a selected
#' set of size k yielded by the lasso--each member of the
#' set is the index of a single feature from a cluster selected by either the
#' protolasso or the cluster representative lasso (the prototype from that
#' cluster--the cluster member most highly correlated with y). (If no set of
#' size k was selected, entry k will be empty.)}
#' @references Bühlmann, P., Rütimann, P., van de Geer, S., & Zhang, C. H.
#' (2013). Correlated variables in regression: Clustering and sparse estimation.
#' \emph{Journal of Statistical Planning and Inference}, 143(11), 1835–1858.
#' \url{https://doi.org/10.1016/j.jspi.2013.05.019}. \cr Jerome Friedman, Trevor
#' Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear
#' Models via Coordinate Descent. \emph{Journal of Statistical Software}, 33(1)
#' ' 1-22. URL \url{https://www.jstatsoft.org/v33/i01/}.
clusterRepLasso <- function(x, y, clusters, var_names=NA, nlambda=4000){

    # Handle and format inputs
    ret <- processClusterLassoInputs(x, y, var_names, clusters)

    x <- ret$x
    var_names_provided <- ret$var_names_provided
    clusters <- ret$clusters
    prototypes <- ret$prototypes
    n_clusters <- ret$n_clusters

    rm(ret)

    n <- nrow(x)

    # Format the design matrix for glmnet according to the cluster
    # representative lasso procedure
    X_glmnet <- getXglmnet(x, clusters, type="clusterRepLasso")

    #  <- getXglmnet_results$X_glmnet
    # cluster_members <- getXglmnet_results$cluster_members
    # non_cluster_feats <- getXglmnet_results$non_cluster_feats

    # rm(getXglmnet_results)

    stopifnot(nrow(X_glmnet) == n)

    # Estimate the lasso on the cluster representatives
    fit <- glmnet::glmnet(x=X_glmnet, y=y, family="gaussian", alpha=1, nlambda=nlambda)
    lasso_sets <- unique(glmnet::predict.glmnet(fit, type="nonzero"))

    # Finally, extract the desired information from the lasso fit--all the
    # sets of selected clusters (one for each observed model size), and
    # corresponding sets of selected features
    cluster_sel_results <- getClusterSelsFromGlmnet(lasso_sets, clusters,
        prototypes, averaging=TRUE)

    return(list(selected_sets=cluster_sel_results$selected_sets,
        selected_clusts_list=cluster_sel_results$selected_clusts_list,
        beta=fit$beta))
}