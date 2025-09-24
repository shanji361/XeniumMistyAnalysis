### (Optional) Heatmap Customization Function 

# Because MISTy does not support customizable heatmap colors, the plot_misty_heatmap() function was developed to allow independent visualization of interaction results. It directly accesses the raw data from the results object, bypassing MISTy’s built-in plotting functions. When collect_results() is executed, MISTy stores importance values in structured data frames, including $importances.aggregated, which contains predictor-target interaction strengths for each view.

# plot_misty_heatmap() filters this data for a specified view and generates a heatmap using ggplot2, with predictors on the X-axis, targets on the Y-axis, and fill colors representing importance values. The function supports customization of color gradients, midpoint values, and applies a slight tilt to the X-axis labels for readability. Because the plot is built entirely from scratch, users can further modify visual elements using additional ggplot2 commands. Example usage: a basic heatmap for the view "para.50" can be generated with plot_misty_heatmap(misty_results_complete_linear, view_name = "para.50"), while a customized version with green-to-red coloring and a specific midpoint can be created using plot_misty_heatmap(misty_results_complete_linear, view_name = "para.50", low_color = "blue", mid_color = "green", high_color = "orange", midpoint = 0.1). 

# Function to create a customized colored MISTy interaction heatmap:

#  plot_misty_heatmap <- function(misty_results, view_name, low_color = "blue", mid_color = "white", high_color = "red", midpoint = 0) {
#    
#    # Filter interaction data for the specified view
#    interaction_data <- misty_results$importances.aggregated %>%
#      filter(view == view_name)
#    
#    # Create heatmap
#    ggplot(interaction_data, aes(x = Predictor, y = Target, fill = Importance)) +
#      geom_tile() +
#      scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = midpoint) +
#      theme_minimal() +
#      theme(axis.text.x = element_text(angle = 45, hjust = 1))
#  }
