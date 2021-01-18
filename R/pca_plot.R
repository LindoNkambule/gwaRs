#' PCA Plot
#'
#' Creates a PCA plot
#'
#' Creates a PCA plot from PLINK pca output, EIGENSTRAT smartpca,
#'   (or any tab-delimited file or data.frame with the same format as PLINK pca or EIGENSTRAT smartpca output).
#'
#' @param data PLINK pca or EIGENSTRAT smartpca output,
#'   (or any tab-delimited file or data.frame with the same format as PLINK pca or EIGENSTRAT smartpca output)
#' @param xComponent A character vector indicating the principal component value to use for the x-axis.
#'   Default is "PC1"
#' @param yComponent A character vector indicating the principal component value to use for the x-axis.
#'   Default is "PC2"
#' @param legendPos A character vector indicating the legend position. Default is "right".
#' @param soft A character vector indicating the software output format. Default is "PLINK".
#'   If you have a tab-delimited file or data.frame with the same format as EIGENSTRAT output, use "EIGENSTRAT"
#' @param colPalette A character vector indicating the color palette to use. Default is "Accent".
#' @param title A string denoting the title to use for the plot. Default is 'PCA Plot'
#'
#' @importFrom data.table fread
#' @importFrom grDevices colorRampPalette
#' @import dplyr
#' @import ggplot2
#' @import RColorBrewer
#'
#' @author Lindokuhle Nkambule
#'
#' @return A PCA plot.
#'
#' @examples
#' pca_plot(pcaData)
#'
#' @export

pca_plot <-
function(data, xComponent = "PC1", yComponent = "PC2", legendPos = "right", soft = "PLINK",
         colPalette = "Accent", title = NULL){

  # binding for global variable(s)
  POP <- NULL

  input_type <- typeof(data)

  if(input_type == "list"){
    eigenvec_file <- data
  }else if(input_type == "character"){
    if(!file.exists(data)){
      stop(paste("File", data, "does not exist", sep = " "))
    }else{
      eigenvec_file <- data.table::fread(data, header = FALSE)
    }
  }else{
    stop(paste("unsupported data input"))
  }

  # data header
  switch(soft,
         PLINK = {names(eigenvec_file)[1:12] = c("POP", "Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")},
         EIGENSTRAT = {names(eigenvec_file)[1:12] = c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "POP")}
  )

  # select only the specified eigenvectors
  eigenvec_file <- eigenvec_file %>% dplyr::select("Sample", all_of(xComponent), all_of(yComponent), "POP")

  # point shapes (for files with a large number of populations)
  point_shapes <- c(21, 22, 23, 24, 25, 21, 22, 23, 25, 25,
                    21, 22, 23, 24, 25, 21, 22, 23, 25, 25,
                    21, 22, 23, 24, 25, 21, 22, 23, 25, 25,
                    21, 22, 23, 24, 25, 21, 22, 23, 25, 25,
                    21, 22, 23, 24, 25, 21, 22, 23, 25, 25,
                    21, 22, 23, 24, 25, 21, 22, 23, 25, 25)

  # total number of populations
  number_of_populations <- nrow(data.frame(table(eigenvec_file$POP)))

  # color palette for inputs with many populations
  colourCount = length(unique(eigenvec_file$POP))

  switch(colPalette,
         Paired = {getPalette = grDevices::colorRampPalette(brewer.pal(12, "Paired"))},
         Set1 = {getPalette = grDevices::colorRampPalette(brewer.pal(9, "Set1"))},
         Set2 = {getPalette = grDevices::colorRampPalette(brewer.pal(8, "Set2"))},
         Set3 = {getPalette = grDevices::colorRampPalette(brewer.pal(12, "Set3"))},
         Dark2 = {getPalette = grDevices::colorRampPalette(brewer.pal(8, "Dark2"))},
         Accent = {getPalette = grDevices::colorRampPalette(brewer.pal(8, "Accent"))}
  )

  if(is.null(title)){
    pcatitle <- paste('PCA Plot')
  }else{
    pcatitle<- title
  }

  # ggplot
  ggplot2::ggplot(eigenvec_file, aes_string(x = xComponent, y = yComponent, shape = "POP", fill = "POP")) +
    geom_point(aes(fill = POP)) +
    scale_shape_manual(values = rep(point_shapes[1:number_of_populations], len = number_of_populations)) +
    scale_fill_manual(values = getPalette(colourCount)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(x = xComponent, y = yComponent) +
    theme_bw() +
    theme(legend.position = legendPos) +
    theme(legend.title=element_blank()) +
    theme(legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", colour = "black")) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    ggtitle(pcatitle) +
    theme(plot.title = element_text(hjust=0.5))

}
