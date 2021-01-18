#' QQ Plot
#'
#' Creates a Q-Q plot
#'
#' Creates a Q-Q plot from PLINK assoc output (or any tab-delimited file or data frame with "P" column).
#'
#' @param data PLINK assoc output, tab-delimited, or a data.frame with "P" column.
#' @param point_col A character vector indicating the color to use for the SNP p-values. Default is "black".
#' @param diag_col A character vector indicating the color to use for the diagonal line. Default is "red".
#' @param diag_line A character vector indicating the line type to use for the diagonal line. Default is "solid".
#' @param title A string denoting the title to use for the plot. Default is 'Q-Q Plot'
#'
#' @importFrom  data.table fread
#' @importFrom  stats median
#' @import ggplot2
#'
#' @author Lindokuhle Nkambule
#'
#' @return A Q-Q plot.
#'
#' @examples
#' qq_plot(gwasData)
#'
#' @export

qq_plot <-
function(data, point_col = "black", diag_col = "red", diag_line = "solid", title = NULL){

  # error handling
  #1. input data type
  input_type <- typeof(data)
  if(input_type == "list"){
    data <- data
  }else if(input_type == "character"){
    if(!file.exists(data)){
      stop(paste("File", data, "does not exist", sep = " "))
    }else{
      data <- fread(data, header = TRUE)
    }
  }else if(input_type == "double"){
    data <- data.frame(data)
    names(data)[1] <- "P"
  }else{
    stop(paste("unsupported data input"))
  }

  # check if the file has the correct 'P' header label
  if(!("P" %in% names(data))){
    stop("The file you provided does not have the 'P' column header label. Check the Information tab for required file format.")
  }

  # data
  observed_p <- sort(data$P)
  log_obs <- -(log10(observed_p))
  expected_p <- c(1:length(observed_p))
  log_exp <- -(log10(expected_p / (length(expected_p)+1)))
  xlim <- as.numeric(as.integer(max(log_obs) + 1))
  ylim <- xlim

  if(is.null(title)){
    qqtitle <- 'QQ Plot'
  }else{
    qqtitle<- title
  }

  # plot
  if(!("CHISQ" %in% names(data))){

    qqplot <- ggplot2::ggplot(data, aes(x = log_exp, y = log_obs)) +
      geom_point(size = 0.5, color = point_col) + geom_abline(color = diag_col, linetype = diag_line) +
      theme_classic() +
      ggtitle(qqtitle) +
      theme(plot.title = element_text(hjust=0.5)) +
      labs(x = expression(paste("Expected ", -log[10]~italic((p)))),
           y = expression(paste("Observed ", -log[10]~italic((p))))) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, xlim)) + scale_y_continuous(expand = c(0, 0), limits = c(0, ylim))

  }else{

    gc <- stats::median(data$CHISQ)/0.455

    qqplot <- ggplot2::ggplot(data, aes(x = log_exp, y = log_obs)) +
      geom_point(size = 0.5, color = point_col) + geom_abline(color = diag_col, linetype = diag_line) +
      theme_classic() +
      ggtitle(qqtitle) +
      theme(plot.title = element_text(hjust=0.5)) +
      labs(x = expression(paste("Expected ", -log[10]~italic((p)))),
           y = expression(paste("Observed ", -log[10]~italic((p))))) +
      annotate("text", x = 2 , y = 5, label = paste("lambda[GC] ==", gc, sep = " "), parse = TRUE) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, xlim)) + scale_y_continuous(expand = c(0, 0), limits = c(0, ylim))
  }

  qqplot

}
