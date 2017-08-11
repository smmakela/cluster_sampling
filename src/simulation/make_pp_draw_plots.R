make_pp_draw_plots <- function(Mj_new_df, ybar_new, ybar_true,
                               num_clusters, num_units, use_sizes, outcome_type,
                               size_model, rootdir, stanmod_name,
                               pop_cluster_data, sampled_cluster_data) {
  # num_clusters -- number of clusters to sample
  # num_units -- number of units to sample
  # use_sizes -- 0/1 for whether cluster sizes used in pop data
  # outcome_type -- whether outcome is continuous or binary
  # rootdir -- root directory where src, output, etc folders are
  # simno -- current iteration; used so that multiple instances aren't trying to write to the same file
  # stanmod -- compiled stan model
  # stanmod_name -- string for name of stan model so we know which parts of the code to run
  # sim_data -- data to use for simulation
  # num_iter -- number of iterations stan should run for
  # num_chains -- number of chains to run in stan

  if (num_units <= 1) {
    nunits <- paste(100*num_units, "pct", sep = "")
  } else {
    nunits <- num_units
  }

  J <- length(unique(Mj_new_df$cluster_id))
  K <- num_clusters

  pop_cluster_data <- dplyr::arrange(pop_cluster_data, cluster_id)
  sampled_cluster_data <- dplyr::arrange(sampled_cluster_data, cluster_id)
print("str(Mj_new_df)")
print(str(Mj_new_df))
  Mj_pop <- pop_cluster_data$Mj
  Mj_sample <- sampled_cluster_data$Mj

  ##########################################
  ### Plot draws of cluster sizes vs truth
  ##########################################
  if ((grepl("bb", stanmod_name) || grepl("lognormal", stanmod_name) ||
       grepl("negbin", stanmod_name))) {
    Mj_new_df$in_sample <- Mj_new_df$cluster_id <= K
    tmpdf <- data.frame(cluster_id = c(1:J), draw_num = 9999,
                        Mj_new = Mj_pop, in_sample = FALSE, yj_new = NA)
    Mj_new_df <- rbind(Mj_new_df, tmpdf)
    Mj_new_df$draw_num <- as.integer(Mj_new_df$draw_num)
    Mj_new_df$is_truth <- ifelse(Mj_new_df$draw_num == 9999, "truth", "draws")
    Mj_sam_df <- data.frame(Mj_sample)
    Mj_pop_df <- data.frame(Mj_pop)
    xmax <- max(Mj_new_df$Mj_new)
    plt <- ggplot(Mj_new_df, aes(x = Mj_new)) +
      geom_line(aes(group = draw_num, colour = is_truth), stat = "density") +
      geom_line(data = Mj_sam_df, stat = "density",
                aes(x = Mj_sample, colour = "sample")) +
      scale_colour_manual("", values = c("truth" = "black", "draws" = "grey50",
                                         "sample" = "grey80")) +
      #scale_x_continuous(limits = c(0, max(Mj_new))) +
      scale_x_continuous(limits = c(0, xmax)) +
      xlab("Cluster size") +
      ylab("Density") +
      ggtitle(paste0("Model: ", stanmod_name, ", ", size_model)) +
      theme_bw()
    ggsave(plt, file = paste0(rootdir, "/output/figures/Mj_draws_usesizes_",
                              use_sizes, "_", outcome_type, "_", size_model, "_",
                              model_name, "_nclusters_", num_clusters,
                              "_nunits_", nunits, "_sim_1.png"),
           width = 10, height = 8)
  }
  print("done plotting Mj_new")
  print(Sys.time())      

  ##########################################
  ### Plot draws of ybar_new
  ##########################################
  tmpdf <- data.frame(ybar_new, truth = ybar_true)
  plt <- ggplot(tmpdf, aes(x = ybar_new)) +
    geom_vline(aes(xintercept = truth), colour = "grey80") +
    geom_line(stat = "density") +
    xlab("Draws of ybar_new") +
    ylab("Density") +
    ggtitle(paste0("Model: ", stanmod_name, ", ", size_model)) +
    theme_bw()
  ggsave(plt, file = paste0(rootdir, "/output/figures/ybar_new_draws_usesizes_",
                            use_sizes, "_", outcome_type, "_", size_model, "_",
                            model_name, "_nclusters_", num_clusters,
                            "_nunits_", nunits, "_sim_1.png"),
         width = 10, height = 8)
  print("done making plots")
  print(Sys.time())

  return(NULL)

}


