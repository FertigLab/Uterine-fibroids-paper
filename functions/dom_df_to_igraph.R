dom_df_to_igraph_lig_rec <- function(sub_df, send_cols, rec_cols, label_additional = NULL,
                                     ligand_col = "#C94D4D", receptor_col = "#2783B2",
                                     min_node = 1, max_node = 5, min_edge = 1, max_edge = 5) {
  nodes = sub_df[ , c("ligand", "receptor", "sending_cluster", "receiving_cluster")]
  nodes = pivot_longer(nodes, names_to = "type", values_to = "values", cols = everything())
  nodes = dplyr::distinct(nodes)
  
  keep_send_cols = send_cols[which(send_cols %in% nodes$values)]
  keep_rec_cols = rec_cols[which(rec_cols %in% nodes$values)]
  
  # Add colors for node properties:
  nodes$values = factor(nodes$values)
  nodes$color = nodes$values
  
  ligands = as.character(nodes$values[which(nodes$type == "ligand")])
  names(ligands) = rep(ligand_col, length(ligands))
  
  receptors = as.character(nodes$values[which(nodes$type == "receptor")])
  names(receptors) = rep(receptor_col, length(receptors))
  
  node_cols <- c(keep_send_cols, keep_rec_cols, ligands, receptors)
  nodes$color = fct_recode(nodes$color, !!!node_cols)
  
  if(!is.null(label_additional)) {
    keep_additional <- label_additional[which(label_additional %in% nodes$values[nodes$type == "ligand"])]
    nodes$color <- as.character(nodes$color)
    nodes$color[nodes$values %in% keep_additional] <- as.character(nodes$values[nodes$values %in% keep_additional])
    nodes$color <- as.factor(nodes$color)
    nodes$color <- fct_recode(nodes$color, !!!keep_additional)
  }
  
  nodes = nodes[ ,c("values", "color")]
  
  # And edges:
  sub_df %>%
    select(sending_cluster, ligand, ligand_mean_z_score) %>%
    distinct() %>%
    rename(from = sending_cluster, 
           to = ligand,
           weight = ligand_mean_z_score) %>%
    mutate(color = from) -> weighted_to_lig
  weighted_to_lig$color <- factor(weighted_to_lig$color)
  weighted_to_lig$color <- fct_recode(weighted_to_lig$color, !!!keep_send_cols)
  
  sub_df %>%
    select(ligand, receptor) %>%
    distinct() %>%
    rename(from = ligand, to = receptor) %>%
    mutate(weight = 2, color = from) -> weighted_lig_to_rec
  weighted_lig_to_rec$color <- factor(weighted_lig_to_rec$color)
  weighted_lig_to_rec$color <- fct_recode(weighted_lig_to_rec$color, !!!ligands)
  
  sub_df %>%
    select(receptor, receiving_cluster, receptor_mean_z_score) %>%
    distinct() %>%
    rename(from = receptor, 
           to = receiving_cluster, 
           weight = receptor_mean_z_score) %>%
    mutate(color = to) -> weighted_rec_to_cl
  weighted_rec_to_cl$color <- factor(weighted_rec_to_cl$color)
  weighted_rec_to_cl$color <- fct_recode(weighted_rec_to_cl$color, !!!keep_rec_cols)
  
  edges = rbind(weighted_to_lig, weighted_lig_to_rec, weighted_rec_to_cl)
  edges$color <- as.character(edges$color)
  #edges$color = factor(edges$from)
  #edges$color = fct_recode(edges$color,
  #                         !!!keep_send_cols,
  #                         !!!ligands,
  #                         !!!receptors)
  
  edges$weight = scales::rescale(edges$weight, to = c(min_edge, max_edge))
  
  # Add size scaled nodes?
  sub_df %>%
    select(sending_cluster) %>%
    distinct() %>%
    rename(node = sending_cluster) %>%
    mutate(weight = 2) -> send_weight
  
  sub_df %>%
    select(receiving_cluster) %>%
    distinct() %>%
    rename(node = receiving_cluster) %>%
    mutate(weight = 2) -> rec_cl_weight
  
  sub_df %>%
    select(ligand) %>%
    distinct() %>%
    rename(node = ligand) %>% 
    mutate(weight = 2) -> ligand_weight
  
  sub_df %>%
    select(receptor) %>%
    distinct() %>%
    rename(node = receptor) %>%
    mutate(weight = 2) -> receptor_weight
  
  all_weights = rbind(send_weight, rec_cl_weight, ligand_weight, receptor_weight)
  nodes$weight = all_weights[match(nodes$values, all_weights$node), "weight"]
  # Need to scale weights so things aren't tooo ridiculous; min to max as function arguments
  #nodes$weight = scales::rescale(nodes$weight, to = c(min_node, max_node))
  nodes$color <- as.character(nodes$color)
  edges$color <- as.character(edges$color)
  
  # Make a graph:
  net = igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  E(net)$weight = edges$weight
  V(net)$size = nodes$weight
  return(net)
}


dom_df_to_igraph_1tomany <- function(sub_df, send_cols, rec_cols, label_additional = NULL,
                                     ligand_col = "#C94D4D", receptor_col = "#2783B2", tf_col = "#44A043",
                                     min_node = 1, max_node = 5, min_edge = 1, max_edge = 5) {
  nodes = sub_df[ , c("ligand", "receptor", "transcription_factor", "sending_cluster", "receiving_cluster")]
  nodes = pivot_longer(nodes, names_to = "type", values_to = "values", cols = everything())
  nodes = dplyr::distinct(nodes)
  
  keep_send_cols = send_cols[which(send_cols %in% nodes$values)]
  keep_rec_cols = rec_cols[which(rec_cols %in% nodes$values)]
  
  # Add colors for node properties:
  nodes$values = factor(nodes$values)
  nodes$color = nodes$values
  
  tfs = as.character(nodes$values[which(nodes$type == "transcription_factor")])
  names(tfs) = rep(tf_col, length(tfs))
  
  ligands = as.character(nodes$values[which(nodes$type == "ligand")])
  names(ligands) = rep(ligand_col, length(ligands))
  
  receptors = as.character(nodes$values[which(nodes$type == "receptor")])
  names(receptors) = rep(receptor_col, length(receptors))
  
  node_cols <- c(keep_send_cols, keep_rec_cols, tfs, ligands, receptors)
  nodes$color = fct_recode(nodes$color, !!!node_cols)
  
  if(!is.null(label_additional)) {
    keep_additional <- label_additional[which(label_additional %in% nodes$values)]
    nodes$color <- as.character(nodes$color)
    nodes$color[nodes$values %in% keep_additional] <- as.character(nodes$values[nodes$values %in% keep_additional])
    nodes$color <- as.factor(nodes$color)
    nodes$color <- fct_recode(nodes$color, !!!keep_additional)
  }
  
  nodes = nodes[ ,c("values", "color")]
  
  # And edges:
  sub_df %>%
    select(sending_cluster, ligand) %>%
    distinct() %>%
    rename(from = sending_cluster, to = ligand) %>%
    mutate(weight = 2) -> weighted_to_lig
  
  sub_df %>%
    select(ligand, receptor) %>%
    distinct() %>%
    rename(from = ligand, to = receptor) %>%
    mutate(weight = 2) -> weighted_lig_to_rec
  
  sub_df %>%
    select(receptor, transcription_factor, data_wide_cor) %>%
    distinct() %>%
    rename(from = receptor, 
           to = transcription_factor, 
           weight = data_wide_cor) -> weighted_rec_to_tf
  weighted_rec_to_tf$weight <- scales::rescale(weighted_rec_to_tf$weight, to = c(min_edge, max_edge))
  
  sub_df %>%
    select(transcription_factor, receiving_cluster, median_tf_auc) %>%
    distinct() %>%
    rename(from = transcription_factor, 
           to = receiving_cluster, 
           weight = median_tf_auc) -> weighted_tf_to_rec
  weighted_tf_to_rec$weight <- scales::rescale(weighted_tf_to_rec$weight, to = c(0.5, 3))
  
  edges = rbind(weighted_to_lig, weighted_lig_to_rec, weighted_rec_to_tf, weighted_tf_to_rec)
  edges$color = factor(edges$from)
  edges$color = fct_recode(edges$color,
                           !!!keep_send_cols,
                           !!!tfs,
                           !!!ligands,
                           !!!receptors)
  #edges = ungroup(edges)
  #edges$weight = scales::rescale(edges$weight, to = c(min_edge, max_edge))
  
  # Add size scaled nodes?
  sub_df %>%
    select(sending_cluster) %>%
    distinct() %>%
    rename(node = sending_cluster) %>%
    mutate(weight = 1) -> send_weight
  
  sub_df %>%
    select(receiving_cluster) %>%
    distinct() %>%
    rename(node = receiving_cluster) %>%
    mutate(weight = 1) -> rec_cl_weight
  
  sub_df %>%
    select(ligand, ligand_mean_z_score) %>%
    distinct() %>%
    rename(node = ligand, 
           weight = ligand_mean_z_score) -> ligand_weight
  
  sub_df %>%
    select(receptor) %>%
    distinct() %>%
    rename(node = receptor) %>%
    mutate(weight = 1) -> receptor_weight
  
  sub_df %>%
    select(transcription_factor) %>%
    distinct() %>%
    rename(node = transcription_factor) %>%
    mutate(weight = 1) -> tf_weight
  
  all_weights = rbind(send_weight, rec_cl_weight, ligand_weight, receptor_weight, tf_weight)
  nodes$weight = all_weights[match(nodes$values, all_weights$node), "weight"]
  # Need to scale weights so things aren't tooo ridiculous; min to max as function arguments
  nodes$weight = scales::rescale(nodes$weight, to = c(min_node, max_node))
  nodes$color <- as.character(nodes$color)
  edges$color <- as.character(edges$color)
  
  # Make a graph:
  net = igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  E(net)$weight = edges$weight
  V(net)$size = nodes$weight
  return(net)
}

dom_df_to_igraph <- function(sub_df, send_cols, rec_cols, label_additional = NULL,
                                     ligand_col = "#C94D4D", receptor_col = "#2783B2", tf_col = "#44A043",
                                     min_node = 1, max_node = 5, min_edge = 1, max_edge = 5) {
  nodes <- sub_df[ , c("ligand", "receptor", "transcription_factor", "sending_cluster", "receiving_cluster")]
  nodes <- pivot_longer(nodes, names_to = "type", values_to = "values", cols = everything())
  nodes = dplyr::distinct(nodes)
  
  keep_send_cols = send_cols[which(send_cols %in% nodes$values)]
  keep_rec_cols = rec_cols[which(rec_cols %in% nodes$values)]
  
  # Add colors for node properties:
  nodes$values = factor(nodes$values)
  nodes$color = nodes$values
  
  tfs = as.character(nodes$values[which(nodes$type == "transcription_factor")])
  names(tfs) = rep(tf_col, length(tfs))
  
  ligands = as.character(nodes$values[which(nodes$type == "ligand")])
  names(ligands) = rep(ligand_col, length(ligands))
  
  receptors = as.character(nodes$values[which(nodes$type == "receptor")])
  names(receptors) = rep(receptor_col, length(receptors))
  
  node_cols <- c(keep_send_cols, keep_rec_cols, tfs, ligands, receptors)
  nodes$color = fct_recode(nodes$color, !!!node_cols)
  
  if(!is.null(label_additional)) {
    keep_additional <- label_additional[which(label_additional %in% nodes$values)]
    nodes$color <- as.character(nodes$color)
    nodes$color[nodes$values %in% keep_additional] <- as.character(nodes$values[nodes$values %in% keep_additional])
    nodes$color <- as.factor(nodes$color)
    nodes$color <- fct_recode(nodes$color, !!!keep_additional)
  }
  
  nodes = nodes[ ,c("values", "color")]
  
  # And edges:
  sub_df %>%
    select(sending_cluster, ligand, ligand_mean_z_score) %>%
    distinct() %>%
    rename(from = sending_cluster, 
           to = ligand,
           weight = ligand_mean_z_score) %>%
    mutate(color = from) -> weighted_to_lig
  weighted_to_lig$color <- factor(weighted_to_lig$color)
  weighted_to_lig$color <- fct_recode(weighted_to_lig$color, !!!keep_send_cols)
  weighted_to_lig$weight <- scales::rescale(weighted_to_lig$weight, to = c(min_edge, max_edge))
  
  sub_df %>%
    select(ligand, receptor) %>%
    distinct() %>%
    rename(from = ligand, to = receptor) %>%
    mutate(weight = 2, color = from) -> weighted_lig_to_rec
  weighted_lig_to_rec$color <- factor(weighted_lig_to_rec$color)
  weighted_lig_to_rec$color <- fct_recode(weighted_lig_to_rec$color, !!!ligands)
  
  sub_df %>%
    select(receptor, transcription_factor, data_wide_cor) %>%
    distinct() %>%
    rename(from = receptor, 
           to = transcription_factor, 
           weight = data_wide_cor) %>%
    mutate(color = from) -> weighted_rec_to_tf
  weighted_rec_to_tf$color <- factor(weighted_rec_to_tf$color)
  weighted_rec_to_tf$color <- fct_recode(weighted_rec_to_tf$color, !!!receptors)
  weighted_rec_to_tf$weight <- scales::rescale(weighted_rec_to_tf$weight, to = c(min_edge, max_edge))
   
  
  sub_df %>%
    select(transcription_factor, receiving_cluster, median_tf_auc) %>%
    distinct() %>%
    rename(from = transcription_factor, 
           to = receiving_cluster, 
           weight = median_tf_auc) %>%
    mutate(color = to) -> weighted_tf_to_rec
  weighted_tf_to_rec$color <- factor(weighted_tf_to_rec$color)
  weighted_tf_to_rec$color <- fct_recode(weighted_tf_to_rec$color, !!!keep_rec_cols)
  weighted_tf_to_rec$weight <- scales::rescale(weighted_tf_to_rec$weight, to = c(0.5, 3))
  
  
  edges = rbind(weighted_to_lig, weighted_lig_to_rec, weighted_rec_to_tf, weighted_tf_to_rec)
  #edges$color = factor(edges$from)
  #edges$color = fct_recode(edges$color,
  #                         !!!keep_send_cols,
  #                         !!!tfs,
  #                         !!!ligands,
  #                         !!!receptors)
  #edges = ungroup(edges)
  edges$weight = scales::rescale(edges$weight, to = c(min_edge, max_edge))
  
  # Add size scaled nodes?
  sub_df %>%
    select(sending_cluster) %>%
    distinct() %>%
    rename(node = sending_cluster) %>%
    mutate(weight = 1) -> send_weight
  
  sub_df %>%
    select(receiving_cluster) %>%
    distinct() %>%
    rename(node = receiving_cluster) %>%
    mutate(weight = 1) -> rec_cl_weight
  
  sub_df %>%
    select(ligand) %>%
    distinct() %>%
    rename(node = ligand) %>%
    mutate(weight = 1) -> ligand_weight
  
  sub_df %>%
    select(receptor) %>%
    distinct() %>%
    rename(node = receptor) %>%
    mutate(weight = 1) -> receptor_weight
  
  sub_df %>%
    select(transcription_factor) %>%
    distinct() %>%
    rename(node = transcription_factor) %>%
    mutate(weight = 1) -> tf_weight
  
  all_weights = rbind(send_weight, rec_cl_weight, ligand_weight, receptor_weight, tf_weight)
  nodes$weight = all_weights[match(nodes$values, all_weights$node), "weight"]
  # Need to scale weights so things aren't tooo ridiculous; min to max as function arguments
  nodes$weight = scales::rescale(nodes$weight, to = c(min_node, max_node))
  nodes$color <- as.character(nodes$color)
  edges$color <- as.character(edges$color)
  
  # Make a graph:
  net = igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  E(net)$weight = edges$weight
  V(net)$size = nodes$weight
  return(net)
}

# Layout in grid form:
net_layout <- function(graph, df) {
  l <- matrix(0, ncol = 2, nrow = length(igraph::V(graph)))
  rownames(l) <- names(igraph::V(graph))
  # get unique items:
  all_ligs = unique(df$ligand)
  all_recs = unique(df$receptor)
  all_tfs = unique(df$transcription_factor)
  all_sends = unique(as.character(df$sending_cluster))
  all_rec_cl = unique(df$receiving_cluster)
  # Vertical gridding
  l[all_ligs, 1] <- -0.75
  l[all_recs, 1] <- 0
  l[all_tfs, 1] <- 0.75
  l[all_rec_cl, 1] <- 1.5
  l[all_sends, 1] <- -1.5
  # Horizontal gridding
  l[all_ligs, 2] <- (seq_along(all_ligs)/mean(seq_along(all_ligs)) - 1) * 5
  l[all_recs, 2] <- (seq_along(all_recs)/mean(seq_along(all_recs)) - 1) * 5
  l[all_tfs, 2] <- (seq_along(all_tfs)/mean(seq_along(all_tfs)) - 1) * 5
  l[all_rec_cl, 2] <- (seq_along(all_rec_cl)/mean(seq_along(all_rec_cl)) - 1) * 5
  l[all_sends, 2] <- (seq_along(all_sends)/mean(seq_along(all_sends)) - 1) * 5
  
  return(l)
}

# Layout in grid form:
net_layout_lig_rec <- function(graph, df) {
  l <- matrix(0, ncol = 2, nrow = length(igraph::V(graph)))
  rownames(l) <- names(igraph::V(graph))
  # get unique items:
  all_ligs = unique(df$ligand)
  all_recs = unique(df$receptor)
  all_sends = unique(as.character(df$sending_cluster))
  all_rec_cl = unique(df$receiving_cluster)
  # Vertical gridding
  l[all_ligs, 1] <- -0.75
  l[all_recs, 1] <- 0
  l[all_rec_cl, 1] <- 0.75
  l[all_sends, 1] <- -1.5
  # Horizontal gridding
  l[all_ligs, 2] <- (seq_along(all_ligs)/mean(seq_along(all_ligs)) - 1) * 5
  l[all_recs, 2] <- (seq_along(all_recs)/mean(seq_along(all_recs)) - 1) * 5
  l[all_rec_cl, 2] <- (seq_along(all_rec_cl)/mean(seq_along(all_rec_cl)) - 1) * 5
  l[all_sends, 2] <- (seq_along(all_sends)/mean(seq_along(all_sends)) - 1) * 5
  
  return(l)
}

plot_network <- function(net, coords, aspect) {
  plot(net, 
       layout = coords, 
       vertex.size = V(net)$weight*2, 
       vertex.label.size = 5,
       vertex.label.dist = 0, 
       #edge.arrow.size = 0.1 * E(net)$weight,
       edge.arrow.size = 0.3,
       #edge.arrow.width = 2*E(net)$weight, 
       edge.arrow.width = 1,
       vertex.label.family = "Arial",
       vertex.label.color = "black", 
       mode = "all", 
       asp = aspect, 
       edge.width = E(net)$weight * 0.5)
}

