add_contri_pie_to_tree_branchnames <- function (plot, tree, signature, branch_names, cex, remove_min) 
{
  fraction = pull(tree@data, signature)
  perc_df = tibble(node = tree@data$node, fraction = fraction)
  perc_df = perc_df[!is.na(perc_df$fraction), ]
  perc_df$other = 1 - perc_df$fraction
  #Add branch name
  if (branch_names){
    branch_names = pull(tree@data, "branch_id")
    perc_df$branch_names <- branch_names
  }
  #Remove branches that have less then remove_min mutations
  perc_df <- perc_df[colSums(branch_mm) >= remove_min,]
  type <- value <- NULL
  cols = 2:3
  ldf <- gather(perc_df, type, value, !!cols) %>% split(., 
                                                        .$node)
  #Plot
  pies = lapply(ldf, function(df) {
    ggplot(df, aes(x = 1, y = value)) +
      geom_segment(y = 0, yend = 1, x = 1.45, xend = 1.45, size = 1) +
      geom_bar(stat = "identity", aes(fill = type)) +
      scale_fill_manual(values = c(fraction = "black", other = "white")) +
      scale_y_continuous(limits = c(0, 1)) +
      coord_polar(theta = "y", clip = "off") +
      theme_inset() +
      labs(tag = df$branch_names[1], size = 2) +
      theme(plot.tag.location = "plot",
            plot.tag.position = c(0.5,0),
            plot.tag = element_text(size = 6))
  })
  ggtree::inset(plot, pies, width = 0.08 * cex, height = 0.08 * 
                  cex, x = "branch")
}

plot_tree_contribution_branchnames <- function (tree, signature, type = "pie", pie_size = 1, branch_names = F, remove_min = 0, ...) 
{
  if (type == "pie") {
    p = plot_gg_tree(tree, ...) + new_scale_color()
    add_contri_pie_to_tree_branchnames(plot = p, tree = tree, signature = signature, 
                                       cex = pie_size, branch_names = branch_names,
                                       remove_min = remove_min)
  }
  else if (type == "color_dot") {
    plot_gg_tree(tree, ...) + geom_point(aes(x = branch, 
                                             color = !!sym(signature)), size = pie_size * 2) + 
      scale_color_gradientn(colors = grad_cols)
  }
  else if (type == "color") {
    plot_gg_tree(tree, branch_color_param = signature, ...)
  }
  else {
    stop("type must be 'pie', 'color', or 'color_dot'")
  }
}
