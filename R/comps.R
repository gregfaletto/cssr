getXglmnet <- function(x, clusters, n_clusters, type, prototypes=NA){
    # Creates design matrix for glmnet by dealing with clusters (for
    # type="protolasso", discards all cluster members except prototype; for
    # type="clusterRepLasso", replaces all cluster members with a simple
    # average of all the cluster members).

    # Check inputs
    checkGetXglmnetInputs(x, clusters, n_clusters, type, prototypes)

    n <- nrow(x)
    p <- ncol(x)

    if(n_clusters > 0){
        for(i in 1:n_clusters){

            if(type == "protolasso"){
                X_cluster_i <- x[, prototypes[i]]
            } else {
                stopifnot(type == "clusterRepLasso")
                if(length(clusters[[i]]) > 1){
                    X_cluster_i <- rowMeans(x[, clusters[[i]]])
                } else{
                    X_cluster_i <- x[, clusters[[i]]]
                }    
            }

            stopifnot(length(X_cluster_i) == n)
            
            if(i == 1){
                X_cluster <- as.matrix(X_cluster_i)
                cluster_members <- clusters[[i]]
            } else{
                X_cluster <- cbind(X_cluster, X_cluster_i)
                cluster_members <- c(cluster_members, clusters[[i]])
            }
        }
        non_cluster_feats <- setdiff(1:p, cluster_members)

        if(length(non_cluster_feats) > 0){
            X_glmnet <- cbind(X_cluster, x[, non_cluster_feats])
        } else{
            X_glmnet <- X_cluster
        }
        stopifnot(ncol(X_glmnet) == length(non_cluster_feats) + n_clusters)
    } else{
        X_glmnet <- x
        if(type == "protolasso"){
            stopifnot(length(prototypes) == 0)
        }
        cluster_members <- integer()
        non_cluster_feats <- 1:p
    }

    colnames(X_glmnet) <- character()

    # Check output
    checkGetXglmnetOutput(n, p, X_glmnet, cluster_members, non_cluster_feats)
    
    return(list(X_glmnet=X_glmnet, cluster_members=cluster_members,
        non_cluster_feats=non_cluster_feats))

}

checkGetXglmnetOutput <- function(n, p, X_glmnet, cluster_members,
    non_cluster_feats){
    stopifnot(is.matrix(X_glmnet))
    stopifnot(nrow(X_glmnet) == n)
    stopifnot(ncol(X_glmnet) <= p)
    stopifnot(ncol(X_glmnet) >= 1)

    stopifnot(is.numeric(cluster_members) | is.integer(cluster_members))
    stopifnot(length(cluster_members) <= p)
    stopifnot(length(cluster_members) >= 0)
    if(length(cluster_members) >= 1){
        stopifnot(length(cluster_members) == length(unique(cluster_members)))
        stopifnot(all(round(cluster_members) == cluster_members))
        stopifnot(all(cluster_members) %in% 1:p)
    }

    stopifnot(is.numeric(non_cluster_feats) | is.integer(non_cluster_feats))
    stopifnot(length(non_cluster_feats) <= p)
    stopifnot(length(non_cluster_feats) >= 0)
    if(length(non_cluster_feats) >= 1){
        stopifnot(length(non_cluster_feats) == length(unique(non_cluster_feats)))
        stopifnot(all(round(non_cluster_feats) == non_cluster_feats))
        stopifnot(all(non_cluster_feats) %in% 1:p)
    }
    
    stopifnot(length(intersect(cluster_members, non_cluster_feats)) == 0)
    stopifnot(length(cluster_members) + length(non_cluster_feats) == p)

}

checkGetXglmnetInputs <- function(x, clusters, n_clusters, type, prototypes){
    stopifnot(is.matrix(x))

    stopifnot(is.list(clusters))
    stopifnot(all(lengths(clusters) >= 1))

    stopifnot(length(n_clusters) == 1)
    stopifnot(is.integer(n_clusters) | is.numeric(n_clusters))
    stopifnot(n_clusters >= 0)
    if(type=="protolasso"){
        stopifnot(n_clusters == length(prototypes))
    }
    stopifnot(n_clusters == length(clusters))

    stopifnot(length(type) == 1)
    stopifnot(is.character(type))
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
    }
}

# TODO(gregfaletto): figure out how to eliminate need for prototypes argument in
# getClusterSelsFromGlmnet when called by clusterRepLasso (shouldn't be
# necessary for anything)
getClusterSelsFromGlmnet <- function(lasso_sets, clusters, prototypes,
    n_cluster_members, non_cluster_feats, var_names_provided=FALSE,
    var_names=NA, averaging=FALSE){

    # Check inputs

    stopifnot(!is.list(clusters) | all(lengths(clusters) >= 1))
    stopifnot(is.list(clusters) | length(clusters) >= 1)

    stopifnot(is.integer(prototypes))
    if(is.list(clusters)){
        stopifnot(length(prototypes) == length(clusters))
    } else{
        stopifnot(length(prototypes) == 1)
    }
    stopifnot(all(!is.na(prototypes)))
    stopifnot(length(prototypes) == length(unique(prototypes)))

    stopifnot(is.numeric(n_cluster_members) | is.integer(n_cluster_members))
    stopifnot(n_cluster_members == round(n_cluster_members))
    stopifnot(n_cluster_members >= 0)
    # stopifnot(n_cluster_members <= p)

    stopifnot(is.numeric(non_cluster_feats) | is.integer(non_cluster_feats))
    # stopifnot(length(non_cluster_feats) <= p)
    stopifnot(length(non_cluster_feats) >= 0)
    if(length(non_cluster_feats) >= 1){
        stopifnot(length(non_cluster_feats) == length(unique(non_cluster_feats)))
        stopifnot(all(round(non_cluster_feats) == non_cluster_feats))
        # stopifnot(all(non_cluster_feats) %in% 1:p)
    }

    stopifnot(length(var_names_provided) == 1)
    stopifnot(is.logical(var_names_provided))
    if(var_names_provided){
        stopifnot(is.character(var_names))
        if(any(is.na(var_names))){
            stop("must provide var_names (with no NAs) if var_names_provided=TRUE")
        }
    }

    stopifnot(length(averaging) == 1)
    stopifnot(is.logical(averaging))

    if(is.list(clusters)){
        n_clusters <- length(clusters)
    } else{
        n_clusters <- 1
    }
    
    stopifnot(length(prototypes) == n_clusters)




    max_length <- max(vapply(lasso_sets, length, integer(1)))

    selected_sets <- list()
    if(averaging){
        to_avg_list <- list()
        selected_clusts_list <- list()
        weights_list <- list()
    }
    

    for(j in 1:max_length){
        # Lasso selected set of size j
        lasso_sets_j <- lasso_sets[lapply(lasso_sets, length) == j]
        if(length(lasso_sets_j) > 0){
            lasso_set_j <- lasso_sets_j[[1]]
            stopifnot(length(unique(lasso_set_j)) == length(lasso_set_j))

            # Recover features from original feature space: deal with
            # prototype and non-prototype features separately
            proto_inds <- lasso_set_j %in% 1:n_clusters
            lasso_set_j_protos <- lasso_set_j[proto_inds]
            lasso_set_j_non_protos <- lasso_set_j[!proto_inds]

            # Recover indices of non-prototype features

            if(length(lasso_set_j_non_protos) > 0){
                stopifnot(length(lasso_set_j_non_protos) <=
                    length(non_cluster_feats))
                stopifnot(all(lasso_set_j_non_protos - n_clusters %in%
                    1:length(non_cluster_feats)))

                recovered_feats <- non_cluster_feats[lasso_set_j_non_protos -
                    n_clusters]

                stopifnot(length(lasso_set_j_non_protos) ==
                    length(recovered_feats))
                stopifnot(all.equal(!proto_inds, lasso_set_j %in%
                    lasso_set_j_non_protos))
                stopifnot(length(recovered_feats) ==
                    length(unique(recovered_feats)))

                lasso_set_j_non_protos <- recovered_feats

                stopifnot(length(lasso_set_j_non_protos) == sum(!proto_inds))
            }

            stopifnot(length(unique(lasso_set_j_protos)) == length(lasso_set_j_protos))
            stopifnot(length(lasso_set_j_protos) == sum(proto_inds))

            proto_selected_j <- length(lasso_set_j_protos) > 0

            stopifnot(all(lasso_set_j_protos %in% 1:n_clusters))

            if(proto_selected_j){
                lasso_set_j_protos <- prototypes[lasso_set_j_protos]
                stopifnot(length(lasso_set_j_protos) == sum(proto_inds))
            }

            stopifnot(length(lasso_set_j) == length(lasso_set_j_protos) +
                length(lasso_set_j_non_protos))

            lasso_set_j[proto_inds] <- lasso_set_j_protos
            lasso_set_j[!proto_inds] <- lasso_set_j_non_protos

            stopifnot(length(unique(lasso_set_j)) == length(lasso_set_j))
            
            if(var_names_provided){
                stopifnot(max(lasso_set_j) <= length(var_names))
                selected_sets[[j]] <- var_names[lasso_set_j]
            } else{
                selected_sets[[j]] <- lasso_set_j
            }

            if(averaging){
                to_avg_list[[j]] <- logical(j)
                selected_clusts_list[[j]] <- list()
                weights_list[[j]] <- list()
            
                if(proto_selected_j){
                    ind_j <- matrix(FALSE, nrow=length(lasso_set_j),
                        ncol=n_clusters)
                    if(is.list(clusters)){
                        for(k in 1:n_clusters){
                            n_cluster_members_k <- length(clusters[[k]])

                            stopifnot(n_cluster_members_k >= 1)

                            ind_j[, k] <- lasso_set_j == prototypes[k]

                            stopifnot(sum(ind_j[, k]) %in% c(0, 1))

                            if(sum(ind_j[, k]) == 1){
                                to_avg_list[[j]][which(ind_j[, k])] <- TRUE

                                # selected_sets[[j]][lasso_set_j == 1] <- cluster_prototype
                                if(var_names_provided){
                                    selected_clusts_list[[j]][[which(ind_j[, k])]] <-
                                        var_names[clusters[[k]]]
                                } else{
                                    selected_clusts_list[[j]][[which(ind_j[, k])]] <-
                                        clusters[[k]]
                                }

                                weights_list[[j]][[which(ind_j[, k])]] <-
                                    rep(1/n_cluster_members_k, n_cluster_members_k)
                            }
                        }
                    } else{
                        stopifnot(n_clusters <= 1)
                        n_cluster_members_k <- length(clusters)
                        stopifnot(n_cluster_members_k != 1)

                        ind_j[, 1] <- lasso_set_j == prototypes

                        stopifnot(sum(ind_j[, 1]) %in% c(0, 1))

                        if(sum(ind_j[, 1]) == 1){
                            to_avg_list[[j]][which(ind_j[, 1])] <- TRUE

                            # selected_sets[[j]][lasso_set_j == 1] <- cluster_prototype
                            if(var_names_provided){
                                selected_clusts_list[[j]][[which(ind_j[, 1])]] <-
                                    var_names[clusters]
                            } else{
                                selected_clusts_list[[j]][[which(ind_j[, 1])]] <-
                                    clusters
                            }


                            
                            weights_list[[j]][[which(ind_j[, 1])]] <-
                                rep(1/n_cluster_members_k, n_cluster_members_k)

                            stopifnot(length(weights_list[[j]][[which(ind_j[, 1])]]) ==
                                length(selected_clusts_list[[j]][[which(ind_j[, 1])]]))
                        }
                    }
                    
                    stopifnot(sum(ind_j) != 0)
                    stopifnot(all(rowSums(ind_j) <= 1))
                } 
            }
        }
    }

    if(averaging){
        return(list(selected_sets=selected_sets,
            to_avg_list=to_avg_list, selected_clusts_list=selected_clusts_list,
            weights_list=weights_list))
    }
    else{
        return(selected_sets)
    }
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
protolasso <- function(x, y, clusters, var_names=NA, nlambda=4000){

    if(is.data.frame(x)){
        x <- stats::model.matrix(~ ., x)
        x <- x[, colnames(x) != "(Intercept)"]
    }

    stopifnot(is.matrix(x))
    n <- nrow(x)

    stopifnot(is.numeric(y) | is.integer(y))
    stopifnot(n == length(y))

    p <- ncol(x)

    stopifnot(is.numeric(y))

    colnames(x) <- character()

    var_names_provided <- FALSE
    # If var_names is provided, convert avg_feats entries to character vectors
    if(all(!is.na(var_names))){
        stopifnot(all(!is.na(var_names)))
        var_names_provided <- TRUE
        stopifnot(length(var_names) == ncol(x))
    }

    # Check clusters argument
    clusters <- checkCssClustersInput(clusters)

    ### Format clusters into a list where all features are in exactly one
    # cluster (any unclustered features are put in their own "cluster" of size
    # 1).
    clust_names <- as.character(NA)
    # if(!is.null(names(clusters))){
    if(!is.null(names(clusters)) & is.list(clusters)){
        clust_names <- names(clusters)
    }

    cluster_results <- formatClusters(clusters, p=p, clust_names=clust_names,
    	get_prototypes=TRUE, x=x, y=y)

    clusters <- cluster_results$clusters
    # multiple <- cluster_results$multiple
    prototypes <- cluster_results$prototypes

    rm(cluster_results)

    if(is.list(clusters)){
        n_clusters <- length(clusters)
    } else{
        n_clusters <- 1
    }
    
    stopifnot(n_clusters == length(prototypes))

    getXglmnet_results <- getXglmnet(x, clusters, n_clusters,
        type="protolasso", prototypes=prototypes)

    X_glmnet <- getXglmnet_results$X_glmnet
    cluster_members <- getXglmnet_results$cluster_members
    non_cluster_feats <- getXglmnet_results$non_cluster_feats

    rm(getXglmnet_results)

    stopifnot(nrow(X_glmnet) == n)

    fit <- glmnet(x=X_glmnet, y=y, family="gaussian", alpha=1, nlambda=nlambda)
    lasso_sets <- unique(predict(fit, type="nonzero"))

    selected_sets <- getClusterSelsFromGlmnet(lasso_sets, clusters,
        prototypes, n_cluster_members=length(cluster_members),
        non_cluster_feats=non_cluster_feats, 
        var_names_provided=var_names_provided, var_names=var_names,
        averaging=FALSE)

    return(list(selected_sets=selected_sets,
        non_cluster_feats=non_cluster_feats, beta=fit$beta))

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
clusterRepLasso <- function(x, y, clusters, var_names=NA, nlambda=4000){

    if(is.data.frame(x)){
        x <- stats::model.matrix(~ ., x)
        x <- x[, colnames(x) != "(Intercept)"]
    }

    colnames(x) <- character()

    stopifnot(is.matrix(x))
    n <- nrow(x)

    stopifnot(is.numeric(y) | is.integer(y))
    stopifnot(n == length(y))

    p <- ncol(x)

    stopifnot(is.numeric(y))

    var_names_provided <- FALSE
    # If var_names is provided, convert avg_feats entries to character vectors
    if(all(!is.na(var_names))){
        stopifnot(all(!is.na(var_names)))
        var_names_provided <- TRUE
        stopifnot(length(var_names) == ncol(x))
    }

    # Check clusters argument
    clusters <- checkCssClustersInput(clusters)

    ### Format clusters into a list where all features are in exactly one
    # cluster (any unclustered features are put in their own "cluster" of size
    # 1).
    clust_names <- as.character(NA)
    # if(!is.null(names(clusters))){
    if(!is.null(names(clusters)) & is.list(clusters)){
        clust_names <- names(clusters)
    }

    cluster_results <- formatClusters(clusters, p=p, clust_names=clust_names,
    	get_prototypes=TRUE, x=x, y=y)

    clusters <- cluster_results$clusters
    # multiple <- cluster_results$multiple
    prototypes <- cluster_results$prototypes

    rm(cluster_results)

    if(is.list(clusters)){
        n_clusters <- length(clusters)
    } else{
        n_clusters <- 1
    }

    stopifnot(n_clusters == length(prototypes))

    getXglmnet_results <- getXglmnet(x, clusters, n_clusters,
        type="clusterRepLasso")

    X_glmnet <- getXglmnet_results$X_glmnet
    cluster_members <- getXglmnet_results$cluster_members
    non_cluster_feats <- getXglmnet_results$non_cluster_feats

    rm(getXglmnet_results)

    stopifnot(nrow(X_glmnet) == n)

    fit <- glmnet(x=X_glmnet, y=y, family="gaussian", alpha=1, nlambda=nlambda)
    lasso_sets <- unique(predict(fit, type="nonzero"))

    cluster_sel_results <- getClusterSelsFromGlmnet(lasso_sets, clusters,
        prototypes, n_cluster_members=length(cluster_members),
        non_cluster_feats=non_cluster_feats,
        var_names_provided=var_names_provided, var_names=var_names,
        averaging=TRUE)

    return(list(selected_sets=cluster_sel_results$selected_sets,
        to_avg_list=cluster_sel_results$to_avg_list,
        selected_clusts_list=cluster_sel_results$selected_clusts_list,
        weights_list=cluster_sel_results$weights_list, 
        non_cluster_feats=non_cluster_feats, beta=fit$beta))
}