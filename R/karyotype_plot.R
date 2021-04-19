#' Karyotype Plot
#'
#' Creates a Karyotype plot
#'
#' Creates a SNP Karyotype or Density plot from an R dataframe with "CHR" and "BP" columns.
#'
#' @param data A data.frame with "CHR" and "BP"columns.
#' @param density.col A character vector with colors to use for gradients.
#' @param window.size A double precision numeric value indicating the window size.
#' @param title A string denoting the title to use for the plot. Default is 'Manhattan Plot'
#'
#' @importFrom stats median setNames
#' @import dplyr
#' @import ggplot2
#' @import scales
#'
#' @author Lindokuhle Nkambule
#'
#' @return A SNP Karyotype plot.
#'
#' @examples
#' karyotype_plot(gwasData)
#'
#' @export

karyotype_plot <-
  function(data, density.col=c("darkgreen","yellow","red"),
           window.size=1e6, title=NULL){

    # binding for global variable(s)
    BP <- CHR <- chromStart <- cont <- window <- . <- NULL

    # error handling
    colheaders <- c("CHR", "BP")
    for(header in colheaders){
      if(!(header %in% colnames(data))){
        stop(paste("The data you provided does not have the", header, "column header. Check the Information tab for required file format", sep = " "))
      }
    }

    plt_data <- data %>% dplyr::select(CHR, BP) %>% mutate(chromStart = BP, chromEnd = BP) %>%
      select(-BP) %>% dplyr::rename(chrom = CHR)

    plt_data$chrom <- as.character(plt_data$chrom)

    chroms_in_data <- unique(plt_data$chrom)

    rename_chromosomes <- c("X", "Y", "M")

    for (chrom in rename_chromosomes){
      if (chrom %in% chroms_in_data){
        if(chrom == "X"){
          plt_data$chrom[plt_data$chrom == "X"] <- 23
        }
        if(chrom == "Y"){
          plt_data$chrom[plt_data$chrom == "Y"] <- 24
        }
        if(chrom == "M"){
          plt_data$chrom[plt_data$chrom == "M"] <- 25
        }
      }
    }

    plt_data$chrom <- as.integer(plt_data$chrom)

    sizes <- as.data.frame(table(plt_data$chrom))
    names(sizes) <- c("chrom", "size")

    # create an ordered factor level to use for the chromosomes in all the datasets
    sizes <- sizes %>% dplyr::mutate(chrom_order = paste0("chr", chrom, sep=""))

    chrom_order <- sizes[['chrom_order']]
    chrom_order_renames <- c("chr23", "chr24", "chr25")
    for (chrom in chrom_order_renames){
      if (chrom %in% chrom_order){
        if(chrom == "chr23"){
          chrom_order[chrom_order=="chr23"] <- "chrX"
        }
        if(chrom == "chr24"){
          chrom_order[chrom_order=="chr24"] <- "chrY"
        }
        if(chrom == "chr25"){
          chrom_order[chrom_order=="chr25"] <- "chrM"
        }
      }
    }

    chrom_key <- stats::setNames(object = as.character(c(1:nrow(sizes))),
                          nm = chrom_order)

    plt_data$chrom <- as.character(plt_data$chrom)
    plt_data <- plt_data %>% base::merge(sizes, by = "chrom")

    bin <- window.size

    formatterbin <- function(x){
      base::paste(base::round(x/1e6, digits = 0), "Mb", sep = "")
      # x/bin
    }

    # create a window column so we can merge df with window_cont_df
    plt_data <- plt_data %>% dplyr::mutate(window = (.$chromStart %/% bin) + 1)

    window_cont_df <- as.data.frame(plt_data %>%
                                      dplyr::mutate(window = (.$chromStart %/% bin) + 1) %>%
                                      dplyr::group_by(window,chrom) %>%
                                      dplyr::summarise(cont = n(), .groups = 'drop'))

    merged_data <- base::merge(plt_data, window_cont_df, by=c("chrom","window"))

    # This is for making sure that chromosomes like X, Y, M are plotted based on the number of chromosome
    # and not on x = 23, 24, 25 respectively. Doing this because some data might be from chrom1-18, then X, Y
    # so we have to make sure X and Y are plotted as chrom19 and 20, so there's no space.
    merged_chroms <- base::unique(merged_data$chrom)
    remove_chrom_vec <- c("23", "24", "25")
    new_chrom_vec <-  merged_chroms[! merged_chroms %in% remove_chrom_vec]
    chrom_n_without_chars <- length(new_chrom_vec)

    for (chrom in remove_chrom_vec){
      if(chrom == "23"){
        merged_data$chrom[merged_data$chrom == "23"] <- (chrom_n_without_chars + 1)
      }
      if(chrom == "24"){
        merged_data$chrom[merged_data$chrom == "24"] <- (chrom_n_without_chars + 2)
      }
      if(chrom == "25"){
        merged_data$chrom[merged_data$chrom == "25"] <- (chrom_n_without_chars + 3)
      }
    }

    # set up scale for gradient
    if((nrow(merged_data) %% 2) == 0) {
      middle_num <- stats::median(scales::rescale(merged_data$cont))
    } else {
      middle_num <- base::mean(scales::rescale(merged_data$cont))
    }
    max_num <- base::max(scales::rescale(merged_data$cont))
    col.gradient <- c(0, middle_num, max_num)

    if ((bin/1e6 < 0.001) || (bin/1e6 > 100000)) {
      title_window <- base::formatC(bin/1e6,format="e", digits=0)
    }else{
      title_window <- bin/1e6
    }

    if(is.null(title)){
      plt_title <- base::paste("The number of SNPs within ", title_window, "Mb window size", sep="")
    }else{
      plt_title<- title
    }

    # This is for making sure we always have 10 y axis/BP labels
    y_axis_labs_breaks <- base::round(seq(0, max(merged_data$chromEnd), by = max(merged_data$chromEnd)/10),0)

    karyotype_plt <- ggplot2::ggplot(data = merged_data) +
      # geom_bar(stat = "identity", x = chrom, y = barmax, fill = "gray8", alpha = 0.2, width=0.4) +
      #geom_rect(aes(xmin = as.numeric(chrom) - 0.2,
                    #xmax = as.numeric(chrom) + 0.2,
                    #ymax = size, ymin = 0),
                #color = NA, fill = NA) +

      geom_segment(aes(x = as.numeric(chrom) - 0.19, y = chromStart,
                       xend = as.numeric(chrom) + 0.19, yend = chromStart, colour = as.integer(cont))) +
      coord_flip() +

      scale_x_reverse(name = "chromosome", expand = c(0.01, 0),
                      breaks = seq(1, length(chroms_in_data), by = 1), labels = names(chrom_key))  +

      # convert the y scale to specified bin value
      scale_y_continuous(labels = formatterbin, expand = c(0.001, 0),
                         breaks = y_axis_labs_breaks, position = "right") +
      scale_colour_gradientn(colours=density.col,
                             values=col.gradient,
                             space = "Lab", na.value = "grey80") +

      labs(title = plt_title) +
      # black & white color theme
      theme(plot.title = element_text(colour = "black", face = "bold", hjust = 0.5),
            axis.text.x = element_text(colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.title=element_blank(),
            axis.line.x = element_line())

    karyotype_plt

}
