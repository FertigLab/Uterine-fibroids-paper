network2df <- function(dom, sending_cluster = NULL, rec_clusters = NULL) {
  if(is.null(sending_cluster)) {
    sending_cluster <- levels(dom@clusters)
  }
  if(is.null(rec_clusters)) {
    rec_clusters <- levels(dom@clusters)
  }
  #Create a matrix of all ligands expressed by the sending clusters 
  cl_ligands <- all_cluster_ligands(dom, sending_cluster)
  cl_ligands.melt <- cl_ligands %>%
    pivot_longer(
      cols = !ligand, 
      names_to = "sending_cluster", 
      values_to = "ligand_mean_z_score"
    )
  cl_ligands.melt <- cl_ligands.melt[cl_ligands.melt$ligand_mean_z_score > 0, ] #Can replace with filter
  network.list <- map(seq_len(length(rec_clusters)), function(x) {
    cl <- rec_clusters[x]
    tfs <- dom@linkages$clust_tf[[cl]]
    if(length(tfs)) {
      nt <- fetch_rec_lig(dom, cl, tfs)
    }
    else {
      nt <- NULL
    }
    nt
  })
  names(network.list) <- rec_clusters
  network <- list_rbind(network.list, names_to = "receiving_cluster")
  if(nrow(network)) {
    merged.lig_exp <- network %>%
      left_join(cl_ligands.melt, by = join_by("ligand"), relationship = "many-to-many") %>%
      filter(! is.na(ligand_mean_z_score))
    
    #Add receptor expression
    clusts <- unique(merged.lig_exp$receiving_cluster)
    receptors <- unique(merged.lig_exp$receptor)
    cl_receptors <- cluster_receptors_exp(dom, receptors, clusts)
    cl_recs.melt <- cl_receptors %>%
      pivot_longer(
        cols = !receptor, 
        names_to = "receiving_cluster", 
        values_to = "receptor_mean_z_score"
      )
    #cl_recs.melt <- cl_recs.melt %>%
    #  filter(receptor_mean_z_score > 0)
    
    merged.rec_exp <- merged.lig_exp %>%
      left_join(cl_recs.melt, by = join_by("receptor", "receiving_cluster"), relationship = "many-to-one")
    
    #Add median TF AUC
    tfs <- unique(merged.rec_exp$transcription_factor)
    clusts <- unique(merged.rec_exp$receiving_cluster)
    auc.df <- get_median_auc(dom, tfs, clusts)
    auc.df.long <- auc.df %>%
      rownames_to_column(var = "transcription_factor") %>%
      pivot_longer(cols = !transcription_factor,
                   names_to = "receiving_cluster",
                   values_to = "median_tf_auc")
    
    merged.tf_auc <- merged.rec_exp %>%
      left_join(auc.df.long, by = join_by("transcription_factor", "receiving_cluster"), relationship = "many-to-one")
    
    #Add receptor TF correlation
    tf_rec_df <- merged.tf_auc %>%
      select(receiving_cluster, transcription_factor, receptor) %>%
      unique()
    
    cor.df <- get_tf_rec_cor(dom, tf_rec_df)
    
    merged <- merged.tf_auc %>%
      left_join(cor.df, by = join_by("receiving_cluster", "transcription_factor", "receptor"), relationship = "many-to-one")
    
  } else {
    merged <- data.frame()
  }
  
  return(merged)
}

#Fetch receptors and ligands for a set of tfs in a cluster
fetch_rec_lig <- function(dom, cl, tfs) {
  rec_lig.list <- map(seq_len(length(tfs)), function(x) {
    t <- tfs[x]
    r <- dom@linkages$clust_tf_rec[[cl]][[t]]
    if(length(r)) {
      rl.df <- fetch_ligands(dom, r)
    }
    else {
      rl.df <- NULL
    }
    rl.df
  })
  names(rec_lig.list) <- tfs
  rec_lig <- list_rbind(rec_lig.list, names_to = "transcription_factor")
  return(rec_lig)
}

#Fetch ligands from domino object for a set of receptors
fetch_ligands <- function(dom, rec) {
  ligs.list <- map(seq_len(length(rec)), function(x) {
    r <- rec[x]
    l <- dom@linkages$rec_lig[[r]]
    l <- resolve_names(dom, l)
    names(l) <- NULL
    if(length(l)) {
      rl <- data.frame(ligand = l)
    } else {
      rl <- NULL
    }
    rl
  })
  names(ligs.list) <- rec
  ligs <- list_rbind(ligs.list, names_to = "receptor")
  return(ligs)
}


#' Title
#'
#' @param dom 
#' @param sending_cluster 
#'
#' @returns A dataframe of average z-score ligand expression by sending cluster. 
#' The rows correspond to ligand genes and the columns to sending clusters. 
#' The function also outputs additional column labeled "ligand".
#' @export
#'
#' @examples
all_cluster_ligands <- function(dom, sending_cluster) {
  all_lig <- list_c(dom@linkages$rec_lig)
  all_lig <- unique(all_lig)
  all_lig <- all_lig[!all_lig == ""]
  all_lig_names_resolved <- resolve_names(dom, all_lig)
  all_lig_names_resolved <- unique(all_lig_names_resolved)
  if(length(dom@linkages$complexes) > 0){
    all_lig_complexes_resolved_list <- resolve_complexes(dom, all_lig_names_resolved)
    all_lig_names_resolved <- list_c(all_lig_complexes_resolved_list)
  }
  
  lig_genes <- intersect(all_lig_names_resolved, rownames(dom@z_scores))
  cl_ligands <- mean_exp_by_cluster(dom, sending_cluster, lig_genes)
  if(length(dom@linkages$complexes) > 0){
    cl_ligands_coll_list <- avg_exp_for_complexes(cl_ligands, all_lig_complexes_resolved_list)
    #if(length(cl_ligands_coll_list)>1){mat <- do.call(rbind, cl_ligands_coll_list)}
    df <- list_rbind(cl_ligands_coll_list, names_to = "ligand")
    cl_ligands <- df
  } else {
    cl_ligands$ligand <- rownames(cl_ligands)
    cl_ligands <- cl_ligands |> relocate(ligand)
  }
  return(cl_ligands) #filter out the rows with all 0s before returning
}

#' Title
#'
#' @param dom 
#' @param receptors 
#' @param clusts 
#'
#' @returns A dataframe of average z-score receptor expression by receiving cluster. 
#' The rows correspond to receptor genes and the columns to receiving clusters. 
#' The function also outputs additional column labeled "receptor".
#' @export
#'
#' @examples
cluster_receptors_exp <- function(dom, receptors, clusts) {
  all_names_resolved <- receptors
  if(length(dom@linkages$complexes) > 0){
    all_complexes_resolved_list <- resolve_complexes(dom, all_names_resolved)
    all_names_resolved <- list_c(all_complexes_resolved_list)
  }
  rec_genes <- intersect(all_names_resolved, rownames(dom@z_scores))
  cl_recs <- mean_exp_by_cluster(dom, clusts, rec_genes)
  if(length(dom@linkages$complexes) > 0){
    cl_coll_list <- avg_exp_for_complexes(cl_recs, all_complexes_resolved_list)
    #if(length(cl_ligands_coll_list)>1){mat <- do.call(rbind, cl_ligands_coll_list)}
    df <- list_rbind(cl_coll_list, names_to = "receptor")
    cl_recs <- df
  } else {
    cl_recs$receptor <- rownames(cl_recs)
    cl_recs <- cl_recs |> relocate(receptor)
  }
  return(cl_recs)
}

#Average expression for a set of genes over cluster(s)
mean_exp_by_cluster <- function(dom, clusts, genes) {
  #Replace with a dplyr command with cluster as a grouping variable. Will be faster,cleaner. this code may be slowing things down
  #df %>% group_by(clusters) %>% summarise_all(mean) 
  gene_exp_list <- map(seq_len(length(clusts)), function(x) {
    cl <- clusts[x]
    n_cell = length(which(dom@clusters == cl))
    if(n_cell > 1){
      sig = rowMeans(dom@z_scores[genes, which(dom@clusters == cl)])
    } else if(n_cell == 1){
      sig = dom@z_scores[genes, which(dom@clusters == cl)]
    } else {
      sig = rep(0, length(genes))
      names(sig) = genes
    }
    sig[which(sig < 0)] <- 0
    sig <- as.data.frame(sig)
    colnames(sig) <- cl
    return(sig)
  })
  gene_exp <- list_cbind(gene_exp_list)
  return(gene_exp)
}

resolve_names <- function(dom, genes) {
  rl_map = dom@misc[["rl_map"]]
  genes_resolved <- sapply(genes, function(l){
    int <- rl_map[rl_map$L.name == l, ][1,] 
    if((int$L.name != int$L.gene) & !grepl("\\,", int$L.gene)){
      int$L.gene
    } else { 
      int$L.name
    }
  })
  return(genes_resolved)
}

resolve_complexes <- function(dom, genes) {
  genes_list <- lapply(genes, function(l){
    if(l %in% names(dom@linkages$complexes)){
      return(dom@linkages$complexes[[l]])
    } else {
      return(l)
    }
  })
  names(genes_list) <- genes
  #genes_complex_resolved <- unlist(genes_list)
  return(genes_list)
}

#exp_mat mat of genesXclusters, values are z-scores averaged over the clusters
#complexes_list list similar to dom@linkages$complexes
#genes a vector of gene names present in the data
avg_exp_for_complexes <- function(exp_mat, complexes_list) {
  trim_list <- complexes_list %>% keep(~{all(.x %in% rownames(exp_mat))})
  gene_exp_list <- lapply(seq_along(trim_list), function(x) {
    if (length(trim_list[[x]]) > 1) {
      mean_exp <- exp_mat %>%
        filter(rownames(exp_mat) %in% trim_list[[x]]) %>%
        summarise(across(all_of(colnames(exp_mat)), mean))
      return(mean_exp)
      #return(colMeans(exp_mat[trim_list[[x]], ]))
    } else {
      return(exp_mat[trim_list[[x]], ,drop=FALSE])
    }
  })
  names(gene_exp_list) <- names(trim_list)
  #if(length(gene_exp_list)>1){mat <- do.call(rbind, gene_exp_list)}
  return(gene_exp_list)
}

cor_by_cluster <- function(dom, rec, tf, cluster) {
  rec_resolved <- resolve_complexes(dom, rec)
  cor.list <- map(list_c(rec_resolved), function(r) {
    df <- data.frame(rec = dom@z_scores[r, ], tf = dom@features[tf, ], cluster = dom@clusters)
    df <- df[df$cluster == cluster, ]
    cor <- stats::cor.test(df$rec, df$tf, method = "spearman", alternative = "greater", exact = FALSE)
    cor$estimate
  })
  return(median(list_c(cor.list)))
}

get_tf_rec_cor <- function(dom, df) {
  dat.cor.list <- lapply(seq_len(nrow(df)), function(x) {
    tf <- df[["transcription_factor"]][x]
    rec <- df[["receptor"]][x]
    dom@cor[rec, tf]
  })
  dat.cor <- list_c(dat.cor.list)
  df$data_wide_cor <- dat.cor
  
  cor.list <- lapply(seq_len(nrow(df)), function(x) {
    tf <- df[["transcription_factor"]][x]
    rec <- df[["receptor"]][x]
    cluster <- df[["receiving_cluster"]][x]
    cor_by_cluster(dom, rec = rec, tf = tf, cluster = cluster)
  })
  cor <- list_c(cor.list)
  df$cluster_cor <- cor
  return(df)
}

get_tf_auc <- function(dom, tfs) {
  mytfs <- dom@features[tfs, , drop=FALSE]
  mytfs <- t(mytfs)
  mytfs <- as.data.frame(mytfs)
  mytfs$cluster <- dom@clusters
  return(mytfs)
}

get_median_auc <- function(dom, tfs, groups) {
  df <- get_tf_auc(dom, tfs)
  df <- df %>%
    filter(cluster %in% groups)
  median_auc_by_group <- df %>% 
    group_by(cluster) %>% 
    summarise_all("median") %>%
    column_to_rownames(var = "cluster") %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column() %>% #Need this if the tibble has only 1 column
    column_to_rownames(var = "rowname")
  return(median_auc_by_group)
}