#' Mirrored Manhattan Plot
#'
#' Creates a Mirrored Manhattan Plot for two traits
#'
#' Create a Mirrored Manhattan Plot from a tab-delimited file or data frame with the compulsory
#' columns: "CHR", "SNP", "BP", "P", "Trait" .
#'
#' @param data A tab-delimited or data frame with the compulsory columns:
#'   "CHR", "SNP", "BP", "P", "Trait".
#' @param trait1 A character string of the trait1 as it appears in the input data.
#' @param trait2 A character string of the trait2 as it appears in the input data.
#' @param trait1_chromCols A character vector indicating which colors to alternate
#'   for trait1 chromosomes.
#' @param trait2_chromCols A character vector indicating which colors to alternate
#'   for trait2 chromosomes.
#' @param xlab A character string to be used as the x-axis label.
#' @param title A character string to be used as the plot title
#' @param annotate_trait1_pval If set, trait1 SNPs with p-value less than or equal to this
#'   p-value will be annotated on the plot.
#' @param annotate_trait1_color A character string indicating the color to be used for annotating
#'   trait1 SNPs by p-value
#' @param annotate_trait2_pval If set, trait2 SNPs with p-value less than or equal to this
#'   p-value will be annotated on the plot.
#' @param annotate_trait2_color A character string indicating the color to be used for annotating
#'   trait2 SNPs by p-value
#' @param annotateSNP A character vector of SNPs in your dataset to annotate.
#'   If some of the SNPs are not in your dataset, gwaRs will throw a warning message.
#' @param annotateSNPcolor A character string denoting the color to use for the annotations.
#' @param highlight A character vector of SNPs in the dataset to highlight.
#'   If some of the SNPs are not in your dataset, gwaRs will throw a warning message. Default is NULL.
#' @param highlightcolor A character string denoting the color to use to highlight the SNPs.
#' @param genomewideline_trait1 Where to draw the "genome-wide significant" line for trait1
#' @param genomewideline_trait2 Where to draw the "genome-wide significant" line for trait2
#' @param genomewideline_type A character string denoting the type of line to be used for the
#'   "genome-wide significant" line. This is the same for both traits. Default is dashed.
#' @param genomewideline_color A character string denoting the color to be used for the
#'   "genome-wide significant" line. This is the same for both traits. Default is red.
#' @param suggestiveline_trait1 Where to draw the "suggestive" line for trait1.
#' @param suggestiveline_trait2 Where to draw the "suggestive" line for trait2.
#' @param suggestiveline_type A character string denoting the type of line to be used for the
#'   "suggestive" line. This is the same for both traits. Default is dashed
#' @param suggestiveline_color A character string denoting the color to be used for the
#'   "suggestive" line. This is the same for both traits. Default is blue.
#'
#' @importFrom  data.table fread
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import tidyr
#'
#' @author Lindokuhle Nkambule
#'
#' @return A Mirrored Manhattan plot for two traits.
#'
#' @examples
#' \dontrun{
#' mirrored_man_plot(inputData)
#' }
#'
#' @export

mirrored_man_plot <-
  function(data, trait1 = NULL, trait2 = NULL,
           trait1_chromCols = c("gray66", "grey36"),
           trait2_chromCols = c("steelblue1", "steelblue4"),
           xlab = "Genomic Position (chromosome)", title = "Manhattan Plot",
           annotate_trait1_pval = FALSE, annotate_trait1_color = "red",
           annotate_trait2_pval = FALSE, annotate_trait2_color = "red",
           annotateSNP=NULL, annotateSNPcolor = "red",
           highlight=NULL, highlightcolor="green3",
           genomewideline_trait1=NULL, genomewideline_trait2=NULL,
           genomewideline_type = "dashed", genomewideline_color = "red",
           suggestiveline_trait1=NULL, suggestiveline_trait2=NULL,
           suggestiveline_type = "dashed", suggestiveline_color = "blue"){

  # binding for global variable(s)
  BP <- BPcum <- CHR <- chr_len <- P <- SNP <- tot <- Trait <- x <- y <- label <- . <- NULL

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

  # 2. Number of traits
  if(is.null(trait1) | is.null(trait2)){
    stop("Two traits are required for mirrored Manhattan plots")
  }

  #3. column headers
  colheaders <- c("SNP", "CHR", "BP", "P", "Trait")
  for(header in colheaders){
    if(!(header %in% colnames(data))){
      stop(paste("The file you provided does not have the", header, "column header.
                 Run ?mirrored_man_plot to check the required input format", sep = " "))
    }
  }

  #4. check if specified traits are in the data
  traits <- c(trait1, trait2)
  for(tr in traits){
    if(!(tr %in% unique(data$Trait))){
      stop(paste("The file you provided does not have the", tr, "trait.
                 Only use traits specified in the data", sep = " "))
    }
  }

  #5. Annotation methods
  if((!is.logical(annotate_trait1_pval) & !is.null(annotateSNP)) |
     ((!is.logical(annotate_trait2_pval) & !is.null(annotateSNP)))){
    stop("Cannot annotate SNPs by both p-value and rsid. Choose one method to annotate SNPs")
  }
  data <- data %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(chr_len=max(BP), .groups = 'drop') %>%
    dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    dplyr::left_join(data, ., by=c("CHR"="CHR")) %>%
    dplyr::arrange(CHR, BP) %>%
    dplyr::mutate(BPcum=BP+tot)

  # sort the chromosomes in ascending order
  data <- data[order(as.numeric(as.character(data$CHR))), ]

  t1 <- data %>% filter(Trait == trait1)
  t1 <- t1 %>% mutate(log = -log10(P))
  t1df <-
    with(t1,
         data.frame(CHR = unique(t1$CHR),
                    col = rep(trait1_chromCols, length(unique(t1$CHR)))))
  t1colors <- t1df[!duplicated(t1df$CHR), ] # trait1 chromosome colors
  t1Data <- merge(t1,t1colors, by  = "CHR")

  t2 <- data %>% filter(Trait == trait2)
  t2 <- t2 %>% mutate(log = -(-log10(P)))
  t2df <-
    with(t2,
         data.frame(CHR = unique(t2$CHR),
                    col = rep(trait2_chromCols, length(unique(t2$CHR)))))
  t2colors <- t2df[!duplicated(t2df$CHR), ]
  t2Data <- merge(t2,t2colors, by  = "CHR")

  manData <- rbind(t1Data, t2Data)

  nCHR <- length(unique(manData$CHR))
  axis.set <- manData %>% dplyr::group_by(CHR) %>% dplyr::summarize(
    center=( max(BPcum) + min(BPcum) ) / 2 , .groups = 'drop')
  axis_label <- axis.set[order(as.numeric(as.character(axis.set$CHR))), ]
  x_lab <- xlab
  mantitle <- title
  ylim <- abs(floor(log10(min(manData$P)))) + 1
  annotations <- data.frame(
    x = c(10,10),
    y = c(ylim, -ylim),
    label = c(trait1, trait2)
  )

  if(!is.null(annotateSNP)){
    annnotations <- data.frame(annotateSNP)
    names(annnotations)[1] <- "SNP"
    annotatedSNPS <- annnotations %>% dplyr::left_join(manData, by = "SNP")
    notAnnotated <- annotatedSNPS[is.na(annotatedSNPS$CHR),]
    if(nrow(notAnnotated) > 0){
      warning(paste("Some of the SNPs are not in chromosome",nCHR, ",so they won't be annotated."))
    }
    annotatedSNPS <- annotatedSNPS %>% tidyr::drop_na()
  }

  if(!is.null(highlight)){
    highlighted <- data.frame(highlight)
    names(highlighted)[1] <- "SNP"
    highlightedSNPs <- highlighted %>% left_join(manData, by = "SNP")
    notHighlighted <- highlightedSNPs[is.na(highlightedSNPs$CHR),]
    if(nrow(notHighlighted) > 0){
      warning(paste("Some of the SNPs are not in chromosome",nCHR, ",so they won't be highlighted"))
    }
    highlightedSNPs <- highlightedSNPs %>% tidyr::drop_na()
  }

  manplot <- ggplot2::ggplot(manData, aes(x = BPcum, y = log, color = col, label = SNP)) +
    geom_point(size = 0.5) +

    # optional arguments
    # genomewide line
    {
      if(!is.null(genomewideline_trait1))
        geom_hline(yintercept = genomewideline_trait1, color = genomewideline_color,
                   linetype = genomewideline_type)
    } +
    {
      if(!is.null(genomewideline_trait2))
        geom_hline(yintercept = -(genomewideline_trait2), color = genomewideline_color,
                   linetype = genomewideline_type)
    } +

    # suggestive line
    {
      if(!is.null(suggestiveline_trait1))
        geom_hline(yintercept = suggestiveline_trait1, color = suggestiveline_color,
                   linetype = suggestiveline_type)
    } +
    {
      if(!is.null(suggestiveline_trait2))
        geom_hline(yintercept = -(suggestiveline_trait2), color = suggestiveline_color,
                   linetype = suggestiveline_type)
    } +

    # trait1 annotations
    {
      if(!is.logical(annotate_trait1_pval))
        geom_point(data = manData %>% dplyr::filter(Trait == trait1 & P <= annotate_trait1_pval),
                   color = annotate_trait1_color, size = 1)
    } +
    {
      if(!is.logical(annotate_trait1_pval))
        ggrepel::geom_label_repel(
          data = manData %>% dplyr::filter(Trait == trait1 & P <= annotate_trait1_pval),
          label.size = 0.1, size = 3, color = "black")
    } +

    # trait2 annotations
    {
      if(!is.logical(annotate_trait2_pval))
        geom_point(data = manData %>% dplyr::filter(Trait == trait2 & P <= annotate_trait2_pval),
                   color = annotate_trait2_color, size = 1)
    } +
    {
      if(!is.logical(annotate_trait2_pval))
        ggrepel::geom_label_repel(
          data = manData %>% dplyr::filter(Trait == trait2 & P <= annotate_trait2_pval),
          label.size = 0.1, size = 3, color = "black")
    } +

    # annotate by SNP ID/rsid
    {
      if(!is.null(annotateSNP))
        geom_point(data = annotatedSNPS, color = annotateSNPcolor, size = 1)
    } +

    {
      if(!is.null(annotateSNP))
        ggrepel::geom_label_repel(data = annotatedSNPS, label.size = 0.1, size = 2, color = "black")
    } +

    # highlight SNPs
    {
      if(!is.null(highlight))
        geom_point(data = highlightedSNPs, color = highlightcolor, size = 1)
    } +

    scale_x_continuous(labels = axis_label$CHR, breaks = axis_label$center) +
    scale_y_continuous(breaks = pretty(manData$log), labels = abs(pretty(manData$log))) +
    scale_colour_identity() +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = x_lab, y = bquote(-log[10]~italic((p)))) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 7, vjust = 0.5),
      axis.text.y = element_text(size = 7)
    ) +
    ggtitle(mantitle) +
    theme(plot.title = element_text(hjust=0.5)) +
    geom_hline(yintercept=0, lwd=0.5, colour="black") +
    geom_label(data=annotations, aes(x=x, y=y, label=label),
               color="black",
               fill=c(trait1_chromCols[1], trait2_chromCols[1]),
               size=4, angle=45, fontface="bold")
  manplot
}
