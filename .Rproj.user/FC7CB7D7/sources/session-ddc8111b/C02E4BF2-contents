### Visualize fused matrix in chord diagram


## Preparation for plotting

# List to keep colors of the grid based on the family of traits
col_chord <- list()

for (i in rownames(consensus_total_case)) {
  list_index <- which(sapply(families, function(x) i %in% x))
  col_chord <- append(col_chord,family_colors[[list_index]])
}

# Color of grids
state_col2 = unlist(col_chord)
names(state_col2) = c(rownames(consensus_total_case))

# Trait Family structure
group <- structure(
  rep(names(families), sapply(families, length)),
  names = unlist(families))



## Function to plot chord diagram

cor_plot2 <- function(df, cutoff=0.05,highlight_chord = NULL) {
  
  #cutoff <- quantile(df,0.75)
  
  df_cut <- ifelse(abs(df) > cutoff, df, 0)
  diag(df_cut) <- 0
  
  col_fun <- colorRamp2(c(cutoff, max(df_cut)), c("#C9B2E4", "#3A0F73"))
  
  
  rows <- rownames(df_cut)
  cols <- colnames(df_cut)
  combinations <- expand.grid(Row = rows, Col = cols)
  combinations$Value <- as.vector(df_cut)
  combinations$Color <- col_fun(combinations$Value)
  result <- combinations[, c("Row", "Col", "Color")]
  
  if (!is.null(highlight_chord)) {
    result$Color[result$Row == highlight_chord[1] & result$Col == highlight_chord[2]] <- "#A6DBA0"
    result$Color[result$Row == highlight_chord[2] & result$Col == highlight_chord[1]] <- "#A6DBA0"
  }
  
  circos.clear()
  
  chordDiagram(df_cut, big.gap=2,  symmetric = TRUE, order = rownames(df_cut),
               directional = 1, transparency=0.5,
               group = group, annotationTrack = c("grid", "axis"),
               col = result,
               grid.col = state_col2, preAllocateTracks = list(track.height = mm_h(4), track.margin = c(mm_h(4), 0)),
               link.sort = TRUE, link.decreasing = TRUE)
  
  # Inner circles indicating trait name and sclaes
  circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, niceFacing = TRUE)
  }, bg.border = NA)
  
  traits_in_diagram <- unique(rownames(df_cut)[rowSums(df_cut) != 0])
  
  # Outer circles indicating trait family
  for (n in names(families)) {
    if(length(traits_in_diagram[traits_in_diagram %in% families[[n]]])>0) {
      highlight.sector(traits_in_diagram[traits_in_diagram %in% families[[n]]], track.index = 1, col = family_colors[[n]], 
                       text = n, cex = 0.8, text.col = "white", niceFacing = TRUE)
    }
  }
  
}




## Plot the chord diagram

png("Results/SNF/chord_diagrams_total-euc.png", width = 5000, height = 3000, res = 500)

par(mfrow = c(1, 2), mar = c(1, 1, 1, 1)) 

cor_plot2(consensus_total_case)
title("Unified Consensus Matrix - Case")

cor_plot2(consensus_total_control)
title("Unified Consensus Matrix - Control")

dev.off()


png("Results/SNF/chord_diagrams_diseasetypes-euc.png", width = 8000, height = 4000, res = 500)

par(mfrow = c(2, 4), mar = c(1, 1, 1, 1)) 

cor_plot2(consensus_cancer_case, highlight_chord=c('CA1','MM'))
title("Unified Consensus Matrix - Cancer - Case")

cor_plot2(consensus_ai_case, highlight_chord=c('CA2','A3S'))
title("Unified Consensus Matrix - AI - Case")

cor_plot2(consensus_inf_case, highlight_chord=c('CA1','A3S'))
title("Unified Consensus Matrix - INF - Case")

cor_plot2(consensus_meta_case, highlight_chord=c('THy','A3S'))
title("Unified Consensus Matrix - M - Case")

cor_plot2(consensus_cancer_control)
title("Unified Consensus Matrix - Cancer - Control")

cor_plot2(consensus_ai_control)
title("Unified Consensus Matrix - AI - Control")

cor_plot2(consensus_inf_control)
title("Unified Consensus Matrix - INF - Control")

cor_plot2(consensus_meta_control)
title("Unified Consensus Matrix - M - Control")

dev.off()

