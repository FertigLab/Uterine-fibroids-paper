
tf_vln <- function(dom, regulon, cols = NULL) {
  auc_mat <- get_tf_auc(dom, regulon)
  plot <- ggplot(auc_mat, aes(x = cluster, y = get(regulon), col = cluster)) + 
    geom_violin(trim = TRUE, scale = "width", adjust = 1, fill = NA) + 
    geom_beeswarm(size = 0.1, method = "compactswarm", corral = "random", corral.width = 0.9) +
    #geom_point(aes(x = cluster, y = get(regulon)), position = "jitter", size = 0.5) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + 
    labs(x = NULL, y = regulon) 
  if(! is.null(cols)) {
    myplot <- plot + scale_color_manual(values = cols)
    return(myplot)
  }
  else {
    return(plot)
  }
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
  #df <- as.data.frame(t(df))
  #df$group <- groups
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

get_median_auc2 <- function(df, groups) {
  df <- as.data.frame(t(df))
  df$cluster <- groups
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
