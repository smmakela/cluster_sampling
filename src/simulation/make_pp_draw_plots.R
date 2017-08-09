
  ##########################################
  ### Plot draws of cluster sizes vs truth
  ##########################################
  if ((grepl("bb", stanmod_name) || grepl("lognormal", stanmod_name) ||
       grepl("negbin", stanmod_name)) && simno == 1) {
    Nj_new_df$in_sample <- Nj_new_df$cluster_id <= K
    tmpdf <- data.frame(cluster_id = c(1:J), draw_num = 9999,
                        Nj_new = Nj_pop, in_sample = FALSE)
    Nj_new_df <- rbind(Nj_new_df, tmpdf)
    Nj_new_df$draw_num <- as.integer(Nj_new_df$draw_num)
    Nj_new_df$is_truth <- ifelse(Nj_new_df$draw_num == 9999, "truth", "draws")
    Nj_sam_df <- data.frame(Nj_sample)
    Nj_pop_df <- data.frame(Nj_pop)
    plt <- ggplot(Nj_new_df, aes(x = Nj_new)) +
      geom_line(aes(group = draw_num, colour = is_truth), stat = "density") +
      geom_line(data = Nj_sam_df, stat = "density",
                aes(x = Nj_sample, colour = "sample")) +
      scale_colour_manual("", values = c("truth" = "black", "draws" = "grey50",
                                         "sample" = "grey80")) +
      scale_x_continuous(limits = c(0, max(Nj_new))) +
      xlab("Cluster size") +
      ylab("Density") +
      ggtitle(paste0("Model: ", stanmod_name, ", ", size_model)) +
      theme_bw()
    ggsave(plt, file = paste0(rootdir, "/output/figures/Nj_draws_usesizes_",
                              use_sizes, "_", outcome_type, "_", size_model, "_",
                              model_name, "_nclusters_", num_clusters,
                              "_nunits_", nunits, "_sim_", simno, ".png"),
           width = 10, height = 8)
  }
  print("done plotting Nj_new")
  print(Sys.time())      

  ##########################################
  ### Plot draws of ybar_new
  ##########################################
  if (simno == 1) {
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
                              "_nunits_", nunits, "_sim_", simno, ".png"),
           width = 10, height = 8)
  }
  print("done making plots")
  print(Sys.time())      

