#' Creates and modifies the colour file
#'
#' @param colourADMX file name for the colour file
#' @param Kseq Vector of K to be plotted
#' @keywords Admixture_colour
#' @export
#' @examples
#' \dontrun{Admixture_ModPlot("Q_file", "MAP_file.map", "Sort_file.txt", 2:2)}

ColourFile <- function(colourADMX, Kseq){
  
  # check if colour file exist
  colFile <- NULL
  
  if(!file.exists(colourADMX)){
    
    colMatrix <- matrix(NA, max(Kseq), max(Kseq))
    
    for(r in 1:max(Kseq)) colMatrix[r, 1:r] = 1:r
    
    fwrite(as.data.frame(colMatrix), colourADMX, sep = '\t', col.names = F)
  }
  
  # now read it
  colFile <- fread(colourADMX, header = F)
  
  # if K is now larger
  if (max(Kseq) > NROW(colFile)){
    
    diff = max(Kseq) - NROW(colFile)
    
    # add the columns
    colFile <- cbind(colFile, matrix(NA, NROW(colFile), diff))
    
    # create the extra rows
    colMatrix <- matrix(NA, diff, max(Kseq))
    
    for(r in 1:diff) colMatrix[r, 1:(NROW(colFile) + r)] = 1:(NROW(colFile) + r)
    
    colFile <- as.data.frame(rbind(as.matrix(colFile), colMatrix))
    
    fwrite(colFile, colourADMX, sep = '\t', col.names = F)
  }
}


#' Funtion to generate the Admixture PLot
#'
#' @param Q_is .Q file stem
#' @param Map_is .map/fam file
#' @param Sort_is Order in which the populations should be plotted
#' @param Col_is Filename for the Colour file
#' @param Kseq Vector of K to be plotted
#' @param toPDF If TRUE sends the plot to pdf
#' @param colorPal The colour palette that will be used to plot
#' @param margin defines marigns in plot
#' @param modding If TRUE enters plot modification
#' @param KtoMod Integer: selects the K to modfy
#' @keywords Admixture_plot
#' @export
#' @examples
#' \dontrun{Admixture_ModPlot("Q_file", "MAP_file.map", "Sort_file.txt", 2:2)}

plotAdmixture_MB2 <- function(
  Q_is = "",
  Map_is = "",
  Sort_is = "",
  Col_is = NULL,
  Kseq = 2:4,
  toPDF = F,
  colorPal = NULL,
  margin = .05,
  modding = F,
  KtoMod = 0
){
  
  if(toPDF) pdf(paste(Q_is, ".ADMX.pdf", sep =""), paper = "a4r" )
  
  if(nchar(Map_is) == 0) MAP_is <- paste(Q_is, ".fam", sep = "")
  
  Map <- fread(Map_is, select = 1, col.names = "IDs")
  SortL <- fread(Sort_is, header = F)
  
  nR = NROW(Kseq)
  
  # routine
  MD <- ifelse(modding == F | toPDF == T, 1, 2)
  
  layout(mat = matrix(1:(nR + MD)), heights = c(rep(1, each = nR), rep(.5, MD))) #;layout.show(nR+1)
  par(mai = c(0, 0, margin, 0))
  
  col = 0
  
  # sort POPS and plot bars
  for(K in Kseq){
    Q <- fread(paste(Q_is, ".", K, ".Q", sep = ""))
    Q <- cbind(Map$IDs, Q)
    colnames(Q)[1] <- "IDs"
    Q <- Q[order(factor(Q$IDs, levels = SortL$V1)),]
    
    famChange <- c(0)
    
    for(m in 2:NROW(Q$IDs)) if(Q$IDs[m] != Q$IDs[m-1]) famChange <- c(famChange, m-1)
    
    barplot(t(Q[,-"IDs"]), width = 1, border = NA, col = colorPal[unlist(Col_is[K, 1:K])], axes = F, space = 0)
    text(0, 0.5, K, font = 2, pos = 2)
    abline(v = famChange)
  }
  
  plot(NA, xlim = c(0, NROW(Q)), ylim = c(0,1), axes = F)
  
  MapU <- data.frame(IDs = unique(Q$IDs), Start = famChange)
  
  for(S in 1:NROW(MapU)){
    
    oldV = MapU$Start[S] + (ifelse(S != NROW(MapU), MapU$Start[S+1], NROW(Q)) -  MapU$Start[S])/2
    
    is.odd <- if(S %% 2) 1 else 1.1
    
    text(oldV, 1, adj = c(is.odd, is.odd), MapU$IDs[S], srt = 60, cex = 1-(1-1/nR^0.05), xpd = T)
  }
  
  # ADD EXTRA LINE IF MODDING IS ON
  if(modding == T){
    
    plot(NA, xlim = c(0, KtoMod), ylim = c(0,1), axes = F, xlab = "", ylab = "")
    text(0, 0.5, KtoMod, font = 2, pos = 2)
    #rect(0:(KtoMod-1), 0, 1:KtoMod, 1, col = ColourPalette[unlist(colFile[KtoMod, 1:KtoMod])])
    rect(0:(KtoMod-1), 0, 1:KtoMod, 1, col = colorPal[1:KtoMod])
  }
  
  if(toPDF) dev.off()
}


#' Main funtion to generate and modify the Admixture PLot
#'
#' @param Q_file .Q file stem
#' @param Map_file .map/fam file
#' @param Sort_file Order in which the populations should be plotted
#' @param Kseq Vector of K to be plotted
#' @param ColourPalette A colour palette
#' @param mod_plot Logical. TRUE to
#' @param KtoMod = 3,
#' @param ToPDF = F
#' @keywords Admixture_mod
#' @export
#' @examples
#' \dontrun{Admixture_ModPlot("Q_file", "MAP_file.map", "Sort_file.txt", 2:2)}

Admixture_ModPlot <- function(
  Q_file = "",
  Map_file = "",
  Sort_file = "",
  Kseq = 2:2,
  ColourPalette = NULL,
  mod_plot = F,
  KtoMod = 3,
  ToPDF = F
){
  
  stopifnot(file.exists(Q_file) & file.exists(Map_file) & file.exists(Sort_file))
  
  ColFile.name <- paste(Q_file, "ADMXcolors", sep = ".")
  # checks or creates colour file
  ColourFile(ColFile.name, Kseq)
  # readS it
  ColFile <- fread(ColFile.name, header = F)
  
  # generate color palette if none is provided
  if(is.null(ColourPalette)){
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    ColourPalette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  }
  
  # plot
  plotAdmixture_MB2(Q_is = Q_file, Map_is = Map_file, Sort_is = Sort_file,
                    Col_is = ColFile, colorPal = ColourPalette,
                    Kseq, modding = mod_plot, KtoMod = KtoMod, toPDF = ToPDF & !mod_plot)
  
  # modding?
  if(mod_plot){
    
    # colour to mod
    cols <- as.character(t(as.vector(fread(ColFile.name)[KtoMod, 1:KtoMod])))
    
    # select two colours by mouse clicks
    cat("Select colours to swap!\n")
    coords <- locator(2)
    
    # get colours to swap
    if(between(coords$y[1], 0, 1, T) & between(coords$y[2], 0, 1, T) &
       between(coords$x[1], 0, KtoMod, T) & between(coords$x[2], 0, KtoMod, T)) {
      
      col1 = as.character(ceiling(coords$x[1]))
      col2 = as.character(ceiling(coords$x[2]))
      
      if(col1 != col2){
        cat('swap:', ColourPalette[as.numeric(col1)], "(",col1,")",
            '&', ColourPalette[as.numeric(col2)], "(",col2,")\n")
        
        # swap colours in the colour set and plot again
        cols <- replace(cols, c(which(cols == col1), which(cols == col2)), c(col2, col1))
        
        for(v in 1:KtoMod) ColFile[KtoMod, v] = as.numeric(cols[v])
        
        fwrite(ColFile, ColFile.name, sep = '\t', col.names = F)
        
        plotAdmixture_MB2(Q_is = Q_file, Map_is = Map_file, Sort_is = Sort_file,
                          Col_is = ColFile, colorPal = ColourPalette,
                          Kseq, modding = mod_plot, KtoMod = KtoMod, toPDF = ToPDF)
        
      } else {
        cat("WARNING: same colour selected, no need to swap.\n")
      }
    } else {
      cat("WARNING: coordinates out of reach\n")
    }
  }
}
