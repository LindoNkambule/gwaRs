#' Manhattan Plot
#'
#' Creates a Manhattan plot
#'
#' Creates a Manhattan plot from PLINK assoc output (or any tab-delimited file or data frame with
#' "SNP", "CHR", "BP", and "P" columns).
#'
#' @param data PLINK assoc output, tab-delimited, or a data.frame with "SNP", "CHR", "BP", and "P" columns.
#' @param chromCol A character vector indicating which colors to alternate for the chromosomes.
#' @param genomewideline Where to draw the "genome-wide significant" line. Default
#'   -log10(5e-8). Set to FALSE or F to disable
#' @param suggestiveline Where to draw the "suggestive" line. Default
#'   -log10(1e-5). Set to FALSE or F to disable.
#' @param chromosome An integer indicating which chromosome to plot. Default is "ALL".
#' @param annotatePval If set, SNPs with p-value less than or equal to this p-value will be annotated on the plot.
#' @param annotateSNP A character vector of SNPs in your dataset to annotate.
#'   If some of the SNPs are not in your dataset, gwaRs will throw a warning message.
#' @param annotateCol A string denoting the color to use for the annotations.
#' @param highlight A character vector of SNPs in the dataset to highlight.
#'   If some of the SNPs are not in your dataset, gwaRs will throw a warning message. Default is NULL.
#' @param highlightCol A string denoting the color to use to highlight the SNPs.
#' @param title A string denoting the title to use for the plot. Default is 'Manhattan Plot'
#'
#' @importFrom  data.table fread
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import scales
#' @import tidyr
#'
#' @author Lindokuhle Nkambule
#'
#' @return A Manhattan plot.
#'
#' @examples
#' man_plot(gwasData)
#'
#' @export

man_plot <-
function(data, chromCol=c("gray44", "black"), genomewideline=-log10(5e-8),
                     suggestiveline=-log10(1e-5), chromosome = "ALL",
                     annotatePval=FALSE, annotateSNP=NULL, annotateCol="red",
                     highlight=NULL, highlightCol="green3", title=NULL){

  # binding for global variable(s)
  BP <- BPcum <- CHR <- chr_len <- P <- SNP <- tot <- . <- NULL

  # error handling
  #1. input data type
  input_type <- typeof(data)

  if(input_type == "list"){
    data <- data
  }else if(input_type == "character"){
    if(!file.exists(data)){
      stop(paste("File", data, "does not exist", sep = " "))
    }else{
      data <- data.table::fread(data, header = TRUE)
    }
  }else{
    stop(paste("unsupported data input"))
  }

  #2. column headers
  colheaders <- c("SNP", "CHR", "BP", "P")
  for(header in colheaders){
    if(!(header %in% colnames(data))){
      stop(paste("The file you provided does not have the", header, "column header. Check the Information tab for required file format", sep = " "))
    }
  }

  #3. check chromosomes
  # 23 = X
  # 24 = Y
  # 25 = XY
  # 26 = MT
  chromosomes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                   21, 22, 23, 24, 25, 26)

  for (chrom in unique(data$CHR)){
    if(!(chrom %in% chromosomes)){
      stop(paste("The file you provided has an unrecognized chromosome", chrom, ". The CHR column should be numeric. If you have 'X', 'Y', 'XY', 'MT' chromosome etc. in your file, change them to numbers and try again.", sep = " "))
    }
  }

  #4. if plotting only one chromosome, check if the file has the specified chromosome
  chroms <- unique(data$CHR)
  if (!(chromosome %in% chroms) & !(chromosome == "ALL")){
    stop(paste("Unrecognized chromsome", chromosome, sep = " "))
  }

  #5. Annotations method
  if(!is.logical(annotatePval) & !is.null(annotateSNP)){
    stop(paste("Cannot use both annotatePval and annotateSNP. Choose one method to annotate SNPs"))
  }

  pfile <- data %>%

    dplyr::group_by(CHR) %>%
    dplyr::summarise(chr_len=max(BP), .groups = 'drop') %>%

    dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%

    dplyr::left_join(data, ., by=c("CHR"="CHR")) %>%

    dplyr::arrange(CHR, BP) %>%
    dplyr::mutate(BPcum=BP+tot)

  # sort the chromosomes in ascending order
  sfile <- pfile[order(as.numeric(as.character(pfile$CHR))), ]

  if(chromosome != "ALL"){
    plt_file <- sfile[ which(sfile$CHR == chromosome),]
    nCHR <- length(unique(plt_file$CHR))
    x_lab <- paste("Chromosome", chromosome, "(Mb)", sep = " ")
    ylim <- abs(floor(log10(min(plt_file$P)))) + 3

    # get SNPs to annotate
    if(!is.logical(annotatePval)){
      annotatedSNPS <- plt_file[ which(plt_file$P <= annotatePval), ]
    }else if(!is.null(annotateSNP)){
      annnotations <- data.frame(annotateSNP)
      names(annnotations)[1] <- "SNP"
      annotatedSNPS <- annnotations %>% dplyr::left_join(plt_file, by = "SNP")
      notAnnotated <- annotatedSNPS[is.na(annotatedSNPS$CHR),]
      if(nrow(notAnnotated) > 0){
        warning(paste("Some of the SNPs are not in chromosome",nCHR, ",so they won't be annotated."))
      }
      annotatedSNPS <- annotatedSNPS %>% tidyr::drop_na()
    }

    # highlight selected SNPs
    if(!is.null(highlight)){
      highlighted <- data.frame(highlight)
      names(highlighted)[1] <- "SNP"
      highlightedSNPs <- highlighted %>% left_join(plt_file, by = "SNP")
      notHighlighted <- highlightedSNPs[is.na(highlightedSNPs$CHR),]
      if(nrow(notHighlighted) > 0){
        warning(paste("Some of the SNPs are not in chromosome",nCHR, ",so they won't be annotated."))
      }
      highlightedSNPs <- highlightedSNPs %>% tidyr::drop_na()
    }

    if(is.null(title)){
      mantitle <- paste('Chromosome',nCHR, 'Manhattan Plot')
    }else{
      mantitle<- title
    }

    # default params for plot
    manplot <- ggplot2::ggplot(plt_file, aes(x = BPcum/1e6, y = -log10(P), color = as.factor(CHR), label = SNP)) +
      geom_point(size = 0.5) +

      # optional arguments
      # genomewide line
      {
        if(!is.logical(genomewideline)) geom_hline(yintercept = genomewideline, color = "red", linetype = "dashed")
      } +

      # suggestive line
      {
        if(!is.logical(suggestiveline)) geom_hline(yintercept = suggestiveline, color = "blue", linetype = "dashed")
      } +

      # p-value annotations
      # annotate by p-value
      {
        if(!is.logical(annotatePval)) geom_point(data = annotatedSNPS, color = annotateCol, size = 1)
      } +

      {
        if(!is.logical(annotatePval)) ggrepel::geom_label_repel(data = annotatedSNPS, label.size = 0.1, size = 3, color = "black")
      } +

      # annotate by SNP ID/rsid
      {
        if(!is.null(annotateSNP)) geom_point(data = annotatedSNPS, color = annotateCol, size = 1)
      } +

      {
        if(!is.null(annotateSNP)) ggrepel::geom_label_repel(data = annotatedSNPS, label.size = 0.1, size = 2, color = "black")
      } +

      # highlight SNPs
      {
        if(!is.null(highlight)) geom_point(data = highlightedSNPs, color = highlightCol, size = 1)
      } +
      ##

      scale_y_continuous(expand = c(0,0), limits = c(0, ylim), labels = scales::comma) +
      scale_x_continuous(labels = scales::comma) +
      scale_color_manual(values = rep(chromCol, nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = x_lab, y = bquote(-log[10]~italic((p)))) +
      theme_classic() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5),
        axis.text.y = element_text(size = 7)
      ) +
      ggtitle(mantitle) +
      theme(plot.title = element_text(hjust=0.5))

  }else{
    plt_file <- sfile
    nCHR <- length(unique(plt_file$CHR))
    axis.set <- plt_file %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 , .groups = 'drop')
    axis_label <- axis.set[order(as.numeric(as.character(axis.set$CHR))), ]
    x_lab <- "Genomic Position (chromosome)"
    ylim <- abs(floor(log10(min(plt_file$P)))) + 1

    # get SNPs to annotate
    if(!is.logical(annotatePval)){
      annotatedSNPS <- plt_file[ which(plt_file$P <= annotatePval), ]
    }else if(!is.null(annotateSNP)){
      annnotations <- data.frame(annotateSNP)
      names(annnotations)[1] <- "SNP"
      annotatedSNPS <- annnotations %>% dplyr::left_join(plt_file, by = "SNP")
      notAnnotated <- annotatedSNPS[is.na(annotatedSNPS$CHR),]
      if(nrow(notAnnotated) > 0){
        warning(paste("Some of the SNPs are not in the data, so they won't be annotated."))
      }
      annotatedSNPS <- annnotations %>% left_join(plt_file, by = "SNP") %>% tidyr::drop_na()
    }

    # highlight selected SNPs
    if(!is.null(highlight)){
      highlighted <- data.frame(highlight)
      names(highlighted)[1] <- "SNP"
      highlightedSNPs <- highlighted %>% left_join(plt_file, by = "SNP")
      notHighlighted <- highlightedSNPs[is.na(highlightedSNPs$CHR),]
      if(nrow(notHighlighted) > 0){
        warning(paste("Some of the SNPs are not in the data,so they won't be annotated."))
      }
      highlightedSNPs <- highlightedSNPs %>% tidyr::drop_na()
    }

    if(is.null(title)){
      mantitle <- 'Manhattan Plot'
    }else{
      mantitle<- title
    }

    manplot <- ggplot2::ggplot(plt_file, aes(x = BPcum, y = -log10(P), color = as.factor(CHR), label = SNP)) +
      geom_point(size = 0.5) +

      ## optional arguments
      # genomewide line
      {
        if(!is.logical(genomewideline)) geom_hline(yintercept = genomewideline, color = "red", linetype = "dashed")
      } +

      # suggestive line
      {
        if(!is.logical(suggestiveline)) geom_hline(yintercept = suggestiveline, color = "blue", linetype = "dashed")
      } +

      # p-value annotations
      # annotate by p-value
      {
        if(!is.logical(annotatePval)) geom_point(data = annotatedSNPS, color = annotateCol, size = 1)
      } +

      {
        if(!is.logical(annotatePval)) ggrepel::geom_label_repel(data = annotatedSNPS, label.size = 0.1, size = 3, color = "black")
      } +

      # annotate by SNP ID/rsid
      {
        if(!is.null(annotateSNP)) geom_point(data = annotatedSNPS, color = annotateCol, size = 1)
      } +

      {
        if(!is.null(annotateSNP)) ggrepel::geom_label_repel(data = annotatedSNPS, label.size = 0.1, size = 2, color = "black")
      } +

      # highlight SNPs
      {
        if(!is.null(highlight)) geom_point(data = highlightedSNPs, color = highlightCol, size = 1)
      } +
      ##

      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      scale_x_continuous(labels = axis_label$CHR, breaks = axis_label$center) +
      scale_color_manual(values = rep(chromCol, nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = x_lab, y = bquote(-log[10]~italic((p)))) +
      theme_classic() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5),
        axis.text.y = element_text(size = 7)
      ) +
      ggtitle(mantitle) +
      theme(plot.title = element_text(hjust=0.5))
  }

  manplot

}
