

#------------ Yoon: creating block_to_num and block_to_num_genepix ------------#
linpart_block <- c("A1", "B1", "C1", "D1", "E1", "F1")

block_to_num_genepix <- list()
for (i in 1:length(linpart_block)){
  idf <- filter(linpart_genepixdf, A1 == unique(linpart_genepixdf$A1)[i])
  idf["A1"] <- i
  block_to_num_genepix <- rbind(block_to_num_genepix, idf)
}
block_to_num <- block_to_num_genepix
df_inna_genepix <- list()
for (i in 1:length(unique(block_to_num_genepix$ID))){
  temp_prot <- list()
  unique_setting_count <- block_to_num %>% filter(ID == unique(block_to_num$ID)[i]) %>%
    ungroup() %>% select(exposure_setting) %>% unique() %>% nrow()
  for (j in 1:6){
    max_min_df <- block_to_num %>% filter(ID == unique(block_to_num$ID)[i]) %>% arrange(ID, exposure_setting) %>%
      filter(exposure_setting == unique(block_to_num$exposure_setting[j]))
    prot_dil_min <- min(max_min_df$A1)
    prot_dil_max <- max(max_min_df$A1)
    max_min_row <- max_min_df[1, c("exposure_setting", "ID")]
    max_min_row["min"] <- prot_dil_min
    max_min_row["max"] <- prot_dil_max
    temp_prot <- rbind(temp_prot, max_min_row)
  }
  df_inna <- rbind(df_inna, temp_prot)
}
save(block_to_num_genepix, df_inna_genepix, file = "camgenepix_new.RData")

#==============================================================================#
#------------------------------------------------------------------------------#
#                              protein ranking                                 #   
#                        layout for plots per protein                          #
#------------------------------------------------------------------------------#
#==============================================================================#

#load("~/Downloads/gbl_dataframes3.RData")  # df_inna_gbl, linpart_gbldf, 
#load("~/Downloads/genepix_dataframes3.RData")  # df_inna2, df_inna
#save(linpart_gbldf, linpart_genepixdf, file = "~/Downloads/cams.RData")

# dilution for NoDNA to compare - 800

#------------------------------------------------------------------------------#
#                                     Data                                     #
#------------------------------------------------------------------------------#

setwd("~/Documents/Research/Malaria/")
source("functionsProt.R")
library(RColorBrewer)
cols <- c(brewer.pal(9, "Set1")[-6], brewer.pal(8, "Dark2"))
#cols <- rainbow(20)  # not as nice but shows that settings are ordinal
pardef <- par(no.readonly = TRUE)

bwid <- 0.4
nset <- 6    # number of settings in a camera
ndil <- 6    # number of dilutions
nprot <- 250

#load("cams.RData")  # linpart_gbldf, linpart_genepixdf
fdirgbl  <- "datasets/ArrayCAM/"  
fdirgpix <- "datasets/Genepix/"
#gbl  <- dfNsdGBL(fdirgbl)                           # third dilution for NoDNA
#gbl  <- dfNsdGBL(fdirgbl, match.dil = TRUE)         # dilution matches protein's
gbl  <- merge(dfNsdGBL(fdirgbl), dfNsdGBL(fdirgbl, match.dil = TRUE), 
              by = c("block", "ID", "setting", "res"))
gbl  <- gbl[order(gbl$ID, gbl$setting, gbl$block), ]
#gpix <- dfNsdGenepix(fdirgpix)                      # matches Median results      
#gpix <- dfNsdGenepix(fdirgpix, match.dil = TRUE)
gpix <- merge(dfNsdGenepix(fdirgpix), dfNsdGenepix(fdirgpix, match.dil = TRUE), 
              by = c("block", "ID", "setting", "res"))
gpix <- gpix[order(gpix$ID, gpix$setting, gpix$block), ]

#df1gbl  <- dfLinDet( gbl,  nstd = 2, tol = 2.5)
#df2gbl  <- dfDetLin( gbl,  nstd = 2, tol = 2.5)
df1gbl  <- dfLinDet2(gbl,  nstd = 2, tol = 2.5)
df2gbl  <- dfDet2Lin(gbl,  nstd = 2, tol = 2.5)
#df1gpix <- dfLinDet( gpix, nstd = 2, tol = 2.5)
#df2gpix <- dfDetLin( gpix, nstd = 2, tol = 2.5)
df1gpix <- dfLinDet2(gpix, nstd = 2, tol = 2.5)
df2gpix <- dfDet2Lin(gpix, nstd = 2, tol = 2.5)

dilutions <- c(50, 200, 800, 3200, 12800, 51200)
settings1 <- sort(unique(df1gbl$setting))
settings2 <- sort(unique(df1gpix$setting))

gbl$res <- log(gbl$res)               # res replaced by logs!
if (ncol(gbl) == 5) {                 # to accommodate both versions in plots
  gbl$nstd.x <- gbl$nstd
  gbl$nstd.y <- gbl$nstd
}
dlist1 <- split(gbl, gbl$ID)
gpix$res <- log(gpix$res)
if (ncol(gbl) == 5) {                 # to accommodate both versions in plots
  gpix$nstd.x <- gpix$nstd
  gpix$nstd.y <- gpix$nstd
}
dlist2 <- split(gpix, gpix$ID)
  
#==============================================================================#
#------------------------------------------------------------------------------#
#           Plots: detection on nstd for both versions (1 and 2 col)           #
#------------------------------------------------------------------------------#
#==============================================================================#

#------------------------------------------------------------------------------#
#                       Plots: linpart, then detection                         #
#------------------------------------------------------------------------------#

laymat <- matrix(c(1:12, 0, rep(13, 2), rep(14, 2), 0), 3, byrow = TRUE)
layout(laymat)
nstd <- 2
tol  <- 2.5

plist1 <- split(df1gbl,  df1gbl$ID)
plist2 <- split(df1gpix, df1gpix$ID)
sumord <- sapply(plist1, function(df) sum(df$end - df$start))
maxord <- sapply(plist1, function(df) max(df$end - df$start))  
ord    <- order(sumord, maxord, decreasing = TRUE)  # max first, then sum
prots <- names(plist1[ord])

par(ask = TRUE)
for (prot in prots) {
  cam1 <- split(dlist1[[prot]], dlist1[[prot]]$setting)
  cam2 <- split(dlist2[[prot]], dlist2[[prot]]$setting)
  cam1 <- cam1[order(as.numeric(names(cam1)))]
  cam2 <- cam2[order(as.numeric(names(cam2)))]
  for (j in 1:length(settings1)) {
    par(mar = c(0.5, 2, 0.5, 0))
    plot(1:6, cam1[[j]]$res, type = "b", pch = 19, xaxt = "n",
         col = "deepskyblue4", xlab = "dilution", ylab = "", lwd = 1.5)
    ilin <- linpart(cam1[[j]][1:6, "res"], tol)  # index linear
    lines((1:6)[ilin], cam1[[j]][ilin, "res"], col = "firebrick3", lwd = 1.5)
    idet <- cam1[[j]]$nstd.x >= nstd & cam1[[j]]$nstd.y >- nstd
    if (any(!idet)) idet[which(!idet)[1]:length(idet)] <- FALSE 
    points((1:6)[!idet], cam1[[j]]$res[!idet], pch = 16, lty = 2, col = "white", 
           cex = 1, type = "b")
    legend("topright", bty = "n", 
           legend = c("ArrayCAM", paste("setting", names(cam1[j]))))
  }
  for (j in 1:length(settings2)) {
    par(mar = c(0.5, 2, 0.5, 0))
    plot(1:6, cam2[[j]]$res, type = "b", pch = 19, xaxt = "n",
         col = "deepskyblue4", xlab = "dilution", ylab = "", lwd = 1.5)
    ilin <- linpart(cam2[[j]][1:6, "res"], tol)  # index linear
    lines((1:6)[ilin], cam2[[j]][ilin, "res"], col = "firebrick3", lwd = 1.5)
    idet <- cam2[[j]]$nstd.x >= nstd & cam2[[j]]$nstd.y >- nstd
    if (any(!idet)) idet[which(!idet)[1]:length(idet)] <- FALSE 
    points((1:6)[!idet], cam2[[j]]$res[!idet], pch = 16, lty = 2, col = "white", 
           cex = 1, type = "b")
    legend("topright", bty = "n", 
           legend = c("Genepix", paste("setting", names(cam2[j]))))
  }
  
  par(mar = c(2, 4, 2.5, 1.5))
  df <- plist1[[prot]]
  df <- df[order(df$setting), ]
  plot(NULL, xlim = c(1, 6), ylim = c(1 - bwid, nset + bwid), 
       main = paste( "ArrayCAM,", df$ID[1]), 
       xlab = "dilution", ylab = "exposure setting", xaxt = "n", yaxt = "n")
  axis(2, at = 1:nset, labels = settings1, las = 2)
  axis(1, at = 1:ndil, labels = dilutions)
  for (k in 1:nset) {
    polygon(c(df[k, "start"], df[k, "end"], df[k, "end"], df[k, "start"]),
            rep(k + c(-1, 1)*bwid, each = 2), col = cols[k], border = cols[k])
  }  
  df <- plist2[[prot]]
  df <- df[order(df$setting), ]
  plot(NULL, xlim = c(1, 6), ylim = c(1 - bwid, nset + bwid), 
       main = paste("Genepix,", df$ID[1]), 
       xlab = "dilution", ylab = "exposure setting", xaxt = "n", yaxt = "n")
  axis(2, at = 1:nset, labels = settings2, las = 2)
  axis(1, at = 1:ndil, labels = dilutions)
  for (k in 1:nset) {
    polygon(c(df[k, "start"], df[k, "end"], df[k, "end"], df[k, "start"]),
            rep(k + c(-1, 1)*bwid, each = 2), col = cols[k], border = cols[k])
  }  
}
par(pardef)

#------------------------------------------------------------------------------#
#                       Plots: detection, then linpart                         #
#------------------------------------------------------------------------------#

laymat <- matrix(c(1:12, 0, rep(13, 2), rep(14, 2), 0), 3, byrow = TRUE)
layout(laymat)
nstd <- 2
tol  <- 2.5

plist1 <- split(df2gbl,  df2gbl$ID)
plist2 <- split(df2gpix, df2gpix$ID)
sumord <- sapply(plist1, function(df) sum(df$end - df$start))
maxord <- sapply(plist1, function(df) max(df$end - df$start))  
ord    <- order(sumord, maxord, decreasing = TRUE)  # max first, then sum
prots <- names(plist1[ord])

par(ask = TRUE)
for (prot in prots) {
  cam1 <- split(dlist1[[prot]], dlist1[[prot]]$setting)
  cam2 <- split(dlist2[[prot]], dlist2[[prot]]$setting)
  cam1 <- cam1[order(as.numeric(names(cam1)))]
  cam2 <- cam2[order(as.numeric(names(cam2)))]
  for (j in 1:length(settings1)) {
    par(mar = c(0.5, 2, 0.5, 0))
    plot(1:6, cam1[[j]]$res, type = "b", pch = 19, xaxt = "n",
         col = "deepskyblue4", xlab = "dilution", ylab = "", lwd = 1.5)
    idet <- cam1[[j]]$nstd.x >= nstd & cam1[[j]]$nstd.y >= nstd
    if (any(!idet)) idet[which(!idet)[1]:length(idet)] <- FALSE  
    points((1:6)[!idet], cam1[[j]]$res[!idet], pch = 16, lty = 2, col = "white", 
           cex = 1, type = "b")
    ilin <- linpart(cam1[[j]][idet, "res"], tol)  # index linear
    ilin <- which(idet)[ilin]
    lines((1:6)[ilin], cam1[[j]][ilin, "res"], col = "firebrick3", lwd = 1.5)
    legend("topright", bty = "n", 
           legend = c("ArrayCAM", paste("setting", names(cam1[j]))))
  }
  
  for (j in 1:length(settings2)) {
    par(mar = c(0.5, 2, 0.5, 0))
    plot(1:6, cam2[[j]]$res, type = "b", pch = 19, xaxt = "n",
         col = "deepskyblue4", xlab = "dilution", ylab = "", lwd = 1.5)
    idet <- cam2[[j]]$nstd.x >= nstd & cam2[[j]]$nstd.y >= nstd
    if (any(!idet)) idet[which(!idet)[1]:length(idet)] <- FALSE  
    points((1:6)[!idet], cam2[[j]]$res[!idet], pch = 16, lty = 2, col = "white", 
           cex = 1, type = "b")
    ilin <- linpart(cam2[[j]][idet, "res"], tol)  # index linear
    ilin <- which(idet)[ilin]
    lines((1:6)[ilin], cam2[[j]][ilin, "res"], col = "firebrick3", lwd = 1.5)
    legend("topright", bty = "n", 
           legend = c("Genepix", paste("setting", names(cam2[j]))))
  }
  
  par(mar = c(2, 4, 2.5, 1.5))
  df <- plist1[[prot]]
  df <- df[order(df$setting), ]
  plot(NULL, xlim = c(1, 6), ylim = c(1 - bwid, nset + bwid), 
       main = paste( "ArrayCAM,", df$ID[1]), 
       xlab = "dilution", ylab = "exposure setting", xaxt = "n", yaxt = "n")
  axis(2, at = 1:nset, labels = settings1, las = 2)
  axis(1, at = 1:ndil, labels = dilutions)
  for (k in 1:nset) {
    polygon(c(df[k, "start"], df[k, "end"], df[k, "end"], df[k, "start"]),
            rep(k + c(-1, 1)*bwid, each = 2), col = cols[k], border = cols[k])
  }  
  df <- plist2[[prot]]
  df <- df[order(df$setting), ]
  plot(NULL, xlim = c(1, 6), ylim = c(1 - bwid, nset + bwid), 
       main = paste("Genepix,", df$ID[1]), 
       xlab = "dilution", ylab = "exposure setting", xaxt = "n", yaxt = "n")
  axis(2, at = 1:nset, labels = settings2, las = 2)
  axis(1, at = 1:ndil, labels = dilutions)
  for (k in 1:nset) {
    polygon(c(df[k, "start"], df[k, "end"], df[k, "end"], df[k, "start"]),
            rep(k + c(-1, 1)*bwid, each = 2), col = cols[k], border = cols[k])
  }  
}
par(pardef)

# prot <- "PF14_0326e2s4"  # all Genepix below their own NoDNA
# everything checks out; signal in ArrayCAM, none in Genepix (and not because of variability) 
a <- read.csv("~/Downloads/Genepix/500.csv", as.is = TRUE)
a <- a[c("Block", "ID", "F635.Median...B635")]
names(a) <- c("block", "ID", "res")
a$block <- (a$block + 1) %/% 2
nodna <- a[a$ID == "NoDNA", ]
j <- 1
nd <- nodna$res[nodna$block == j]
nd
a[a$ID == prot & a$block == j, "res"]
c(mean(nd), sd(nd))
j <- j + 1

b <- read.csv("~/Downloads/ArrayCAM/200.csv", as.is = TRUE)
b <- b[c("A1", "ID", "Act_Result")]
names(b) <- c("block", "ID", "res")
for (i in 1:6) {
  b$block[grep(LETTERS[i], b$block)] <- i
}
nodna <- b[b$ID == "NoDNA", ]
j <- 1
nd <- as.numeric(nodna$res[nodna$block == j])
nd
b[b$ID == prot & b$block == j, "res"]
c(mean(nd, na.rm = TRUE), sd(nd, na.rm = TRUE))
j <- j + 1

# should be doing two sample t-test (will be easier with 4 replicates)
#   instead of comparing the mean of Protein to all the NoDNA
# might throw in random effect of the block as well

#==============================================================================#
#------------------------------------------------------------------------------#
#                          Plots: summary by camera                            #
#------------------------------------------------------------------------------#
#==============================================================================#

#------------------------------------------------------------------------------#
#                     proportion of proteins per interval                      #
#------------------------------------------------------------------------------#

par(mfrow = c(1, 2))
#slist <- split(df1gbl, df1gbl$setting)  # linpart, then detection
slist <- split(df2gbl, df2gbl$setting)  # detection, then linpart
plot(NULL, xlim = c(1, 6), ylim = c(1 - bwid, nset + bwid), xaxt = "n", 
     yaxt = "n", main = "ArrayCAM", xlab = "dilution", ylab = "camera setting")
axis(2, at = 1:nset, labels = settings1, las = 2)
axis(1, at = 1:ndil, labels = dilutions)
for (k in 1:nset) {
  vv <- NULL
  for (i in 1:nrow(slist[[k]])) {
    vv <- c(vv, slist[[k]][i, "start"]:(slist[[k]][i, "end"] - 1))
  }
  vv <- vv[vv > 0]
  tbl <- table(vv)
  dil <- as.numeric(names(tbl))
  for (j in 1:length(tbl)) {
    polygon(c(dil[j], dil[j] + 1, dil[j] + 1, dil[j]),
            rep(k + c(-1, 1)*bwid, each = 2), 
            col    = adjustcolor(cols[k], alpha.f = tbl[j]/nprot*1), 
            border = adjustcolor(cols[k], alpha.f = tbl[j]/nprot*1))
    text(dil[j] + 0.5, k, tbl[j]/nprot, cex = 0.8)
  }
}

#slist <- split(df1gpix, df1gpix$setting)  # linpart, then detection
slist <- split(df2gpix, df2gpix$setting)  # detection, then linpart
plot(NULL, xlim = c(1, 6), ylim = c(1 - bwid, nset + bwid), xaxt = "n", 
     yaxt = "n", main = "Genepix camera", 
     xlab = "dilution", ylab = "camera setting")
axis(2, at = 1:nset, labels = settings2, las = 2)
axis(1, at = 1:ndil, labels = dilutions)
for (k in 1:nset) {
  vv <- NULL
  for (i in 1:nrow(slist[[k]])) {
    vv <- c(vv, slist[[k]][i, "start"]:(slist[[k]][i, "end"] - 1))
  }
  vv <- vv[vv > 0]
  tbl <- table(vv)
  dil <- as.numeric(names(tbl))
  for (j in 1:length(tbl)) {
    polygon(c(dil[j], dil[j] + 1, dil[j] + 1, dil[j]),
            rep(k + c(-1, 1)*bwid, each = 2), 
            col    = adjustcolor(cols[k], alpha.f = tbl[j]/nprot*1), 
            border = adjustcolor(cols[k], alpha.f = tbl[j]/nprot*1))
    text(dil[j] + 0.5, k, tbl[j]/nprot, cex = 0.8)
  }
}
par(pardef)

#------------------------------------------------------------------------------#
#                      All combinations of linear portions                     #
#------------------------------------------------------------------------------#

#----------------------- Number of proteins - transparency --------------------#

nset <- length(settings1)               # number of settings
ndil <- length(dilutions)               # number of dilutions
lwd <- 7

slist <- split(df1gbl, df1gbl$setting)  # linpart, then detection
#slist <- split(df2gbl, df2gbl$setting)  # detection, then linpart
nmax  <- 0
nintA  <- nintB  <- numeric(nset)
tlistA <- tlistB <- list(nset)
for (k in 1:nset) {
  t2 <- table(slist[[k]]$start, slist[[k]]$end)[-1, -1]
  if (!is.matrix(t2)) t2 <- matrix(t2, nrow = 1)
  diag(t2) <- 0
  nintA[k]    <- sum(t2 > 0)
  nmax        <- max(nmax, t2)
  tlistA[[k]] <- t2
}
slist <- split(df1gpix, df2gpix$setting)  # linpart, then detection
#slist <- split(df2gpix, df2gpix$setting)  # detection, then linpart
for (k in 1:nset) {
  t2 <- table(slist[[k]]$start, slist[[k]]$end)[-1, -1]
  if (!is.matrix(t2)) t2 <- matrix(t2, nrow = 1)
  diag(t2) <- 0
  nintB[k]    <- sum(t2 > 0)
  nmax        <- max(nmax, t2)
  tlistB[[k]] <- t2
}

par(mfrow = c(1, 2))
plot(NULL, xlim = c(1, 6), ylim = c(1, nset + 1), xaxt = "n", yaxt = "n", 
     main = "ArrayCAM", xlab = "dilution", ylab = "camera setting")
axis(2, at = 1.5:(nset + 0.5), labels = settings1, las = 2, tick = FALSE)
axis(1, at = 1:ndil, labels = dilutions)
for (k in 1:nset) {  
  t2 <- tlistA[[k]]
  count <- 0
  for (i in 1:nrow(t2)) {
    for (j in 1:ncol(t2)) {
      if (t2[i, j] == 0) next 
      lines(c(i, j), rep(k + count/nintA[k], 2), lwd = lwd,  
            col = adjustcolor(cols[k], alpha.f = t2[i,j]/nmax))
      loc <- count/nintA[k]
      text(ndil - 1 + loc, k + loc, t2[i, j], col = cols[k], cex = 0.8)
      count <- count + 1
    }
  }
}
plot(NULL, xlim = c(1, 6), ylim = c(1, nset + 1), xaxt = "n", yaxt = "n", 
     main = "Genepix", xlab = "dilution", ylab = "camera setting")
axis(2, at = 1.5:(nset + 0.5), labels = settings2, las = 2, tick = FALSE)
axis(1, at = 1:ndil, labels = dilutions)
for (k in 1:nset) {  
  t2 <- tlistB[[k]]
  count <- 0
  for (i in 1:nrow(t2)) {
    for (j in 1:ncol(t2)) {
      if (t2[i, j] == 0) next 
      lines(c(i, j), rep(k + count/nintB[k], 2), lwd = lwd,  
            col = adjustcolor(cols[k], alpha.f = t2[i,j]/nmax))
      loc <- count/nintB[k]
      text(ndil - 1 + loc, k + loc, t2[i, j], col = cols[k], cex = 0.8)
      count <- count + 1
    }
  }
}
par(pardef)

# ---------------------- Number of proteins - line width  ---------------------#  

nprot <- 250
nset <- length(settings1)               # number of settings
ndil <- length(dilutions)               # number of dilutions
lwd <- 10

slist <- split(df1gbl, df1gbl$setting)  # linpart, then detection
#slist <- split(df2gbl, df2gbl$setting)  # detection, then linpart
nmax  <- 0
nintA  <- nintB  <- numeric(nset)
tlistA <- tlistB <- list(nset)
for (k in 1:nset) {
  t2 <- table(slist[[k]]$start, slist[[k]]$end)[-1, -1]
  if (!is.matrix(t2)) t2 <- matrix(t2, nrow = 1)
  diag(t2) <- 0
  nintA[k]    <- sum(t2 > 0)
  nmax        <- max(nmax, t2)
  tlistA[[k]] <- t2
}
slist <- split(df1gpix, df2gpix$setting)  # linpart, then detection
#slist <- split(df2gpix, df2gpix$setting)  # detection, then linpart
for (k in 1:nset) {
  t2 <- table(slist[[k]]$start, slist[[k]]$end)[-1, -1]
  if (!is.matrix(t2)) t2 <- matrix(t2, nrow = 1)
  diag(t2) <- 0
  nintB[k]    <- sum(t2 > 0)
  nmax        <- max(nmax, t2)
  tlistB[[k]] <- t2
}

par(mfrow = c(1, 2))
plot(NULL, xlim = c(1, ndil), ylim = c(1, nset + 1), xaxt = "n", yaxt = "n", 
     main = "ArrayCAM", xlab = "dilution", ylab = "camera setting")
axis(2, at = 1.4:(nset + 0.4), labels = settings1, las = 2, tick = FALSE)
axis(1, at = 1:ndil, labels = dilutions)
for (k in 1:nset) {  
  t2 <- tlistA[[k]]
  count <- 0
  for (i in 1:nrow(t2)) {
    for (j in 1:ncol(t2)) {
      if (t2[i, j] == 0) next 
      lines(c(i, j), rep(k + count/nintA[k], 2), lwd = lwd*t2[i, j]/nmax, 
            col = adjustcolor(cols[k], alpha.f = 1)) 
      loc <- count/nintA[k]
      text(ndil - 1 + loc, k + loc, t2[i, j], col = cols[k], cex = 0.8)
      count <- count + 1
    }
  }
  text(ndil, k + 0.4, nprot - sum(t2), col = "darkgrey")
}
plot(NULL, xlim = c(1, ndil), ylim = c(1, nset + 1), xaxt = "n", yaxt = "n", 
     main = "Genepix", xlab = "dilution", ylab = "camera setting")
axis(2, at = 1.4:(nset + 0.4), labels = settings2, las = 2, tick = FALSE)
axis(1, at = 1:ndil, labels = dilutions)
for (k in 1:nset) {  
  t2 <- tlistB[[k]]
  count <- 0
  for (i in 1:nrow(t2)) {
    for (j in 1:ncol(t2)) {
      if (t2[i, j] == 0) next 
      lines(c(i, j), rep(k + count/nintB[k], 2), lwd = lwd*t2[i, j]/nmax, 
            col = adjustcolor(cols[k], alpha.f = 1)) 
      loc <- count/nintB[k]
      text(ndil - 1 + loc, k + loc, t2[i, j], col = cols[k], cex = 0.8) 
      count <- count + 1
    }
  }
  text(ndil, k + 0.4, nprot - sum(t2), col = "darkgrey")
}
par(pardef)

#------------------------------------------------------------------------------#
#                         Length of linear portions only                       #
#------------------------------------------------------------------------------#

#----------------------- Number of proteins - line width ----------------------#

nset <- length(settings1)               # number of settings
ndil <- length(dilutions)               # number of dilutions
lwd <- 23

slist <- split(df1gbl, df1gbl$setting)  # linpart, then detection
#slist <- split(df2gbl, df2gbl$setting)  # detection, then linpart
nmax  <- 0
nintA  <- nintB  <- numeric(nset)
tlistA <- tlistB <- list(nset)
for (k in 1:nset) {
  t1 <- table(slist[[k]]$nlin)[-1]
  nintA[k]    <- length(t1)
  nmax        <- max(nmax, t1)
  tlistA[[k]] <- t1
}
slist <- split(df1gpix, df2gpix$setting)  # linpart, then detection
#slist <- split(df2gpix, df2gpix$setting)  # detection, then linpart
for (k in 1:nset) {
  t1 <- table(slist[[k]]$nlin)[-1]
  nintB[k]    <- length(t1)
  nmax        <- max(nmax, t1)          # nmax - out of both cameras!
  tlistB[[k]] <- t1
}

par(mfrow = c(1, 2))
plot(NULL, xlim = c(1, ndil), ylim = c(1, nset + 1), yaxt = "n", 
     main = "ArrayCAM", xlab = "number of points", ylab = "camera setting")
axis(2, at = 1.5:(nset + 0.5), labels = settings1, las = 2, tick = FALSE)
for (k in 1:nset) {  
  t1 <- tlistA[[k]]
  yloc <- k + ((0:(nintA[k] - 1))/nintA[k])
  segments(x0 = rep(1, nintA[k]), y0 = yloc, x1 = as.integer(names(t1)), 
           y1 = yloc, lwd = lwd*t1/nmax, col = cols[k], lend = 1)  # butt caps
  #           col = adjustcolor(cols[k], alpha.f = t1/nmax))  # not vectorized
  text(ndil, yloc, t1, col = cols[k])
}
plot(NULL, xlim = c(1, ndil), ylim = c(1, nset + 1), yaxt = "n", 
     main = "Genepix", xlab = "number of points", ylab = "camera setting")
axis(2, at = 1.5:(nset + 0.5), labels = settings2, las = 2, tick = FALSE)
for (k in 1:nset) {  
  t1 <- tlistB[[k]]
  yloc <- k + ((0:(nintB[k] - 1))/nintB[k])
  segments(x0 = rep(1, nintB[k]), y0 = yloc, x1 = as.integer(names(t1)), 
           y1 = yloc, lwd = lwd*t1/nmax, col = cols[k], lend = 1)
  #           col = adjustcolor(cols[k], alpha.f = t1/nmax))  # not vectorized
  text(ndil, yloc, t1, col = cols[k])
}
par(pardef)

#------------------------------------------------------------------------------#
#                  Check the overall sum of linear segments                    #
#------------------------------------------------------------------------------#

# across all settings
c(sum(df1gbl$nlin - 1), sum(df1gpix$nlin - 1))
c(sum(df2gbl$nlin - 1), sum(df2gpix$nlin - 1))

# for each setting separately
cbind(GBL  = tapply(df1gbl$nlin  - 1, df1gbl$setting,  sum), 
      Gpix = tapply(df1gpix$nlin - 1, df1gpix$setting, sum))
cbind(GBL  = tapply(df2gbl$nlin  - 1, df2gbl$setting,  sum), 
      Gpix = tapply(df2gpix$nlin - 1, df2gpix$setting, sum))

#==============================================================================#
#------------------------------------------------------------------------------#
#                   NoDNA distribution and signal detection                    #   
#------------------------------------------------------------------------------#
#==============================================================================#

ndblock <- "C"  # dilution 800

fdir <- "~/Downloads/ArrayCAM/"  
flist <- list.files(fdir)
par(mfrow = c(2, 2), ask = TRUE)
for (i in 1:length(flist)) {
  dat <- read.csv(paste(fdir, flist[i], sep = ""), as.is = TRUE)
  dat$res <- suppressWarnings(as.numeric(dat$Act_Result))
  # Act_Result is a factor because of those empty lines in excel!!
  nodna <- dat[tolower(dat$ID) == "nodna" & grepl(ndblock, dat$A1), c("A1", "res")]
  mnd <- mean(nodna$res)
  snd <- sd(nodna$res)
  nodna$pair <- as.integer(substr(nodna$A1, 2, 2))
  plot(nodna$res, rep(1, nrow(nodna)), pch = 16, col = (3:4)[nodna$pair])
  hist(nodna$res)
  hist(nodna$res[nodna$pair == 1], col = adjustcolor(3, alpha.f = 0.2)) 
  hist(nodna$res[nodna$pair == 2], col = adjustcolor(4, alpha.f = 0.2)) 
}  
par(pardef)

# start with these proteins:
strange <- c("MAL8P1.4.2o2", "PF10_0143e1s2", "PFB0305c-e1", "PF14_0114.1o1",
             "PF14_0316e1s1", "PFE0055ce3s1", "PFC0440ce1s3", "PF14_0465.1o1", 
             "PF10_0124-e1s2", "PFB0915we2s2", "PF14_0649-e2s1")
i <- 1                                   # block C is the one from which nodna is taken
df <- gbl[gbl$ID == strange[i], ]
split(df, df$setting)
i <- i + 1

#==============================================================================#  
#------------------------------------------------------------------------------#
#                                     R^2                                      #   
#------------------------------------------------------------------------------#
#==============================================================================#

S <- 1e3; N <- 20
r2 <- matrix(nrow = S, ncol = N)
for (n in 1:N) {
  for (s in 1:S) {
    x <- 1:n
    y <- x + rnorm(n, sd = 0.5)
    fit <- lm(y ~ x)
    r2[s, n] <- summary(fit)$r.squared
  }
}
apply(r2, 2, mean)  # from n = 3 goes up
apply(r2, 2, var)   # from n = 3 goes down

#------------------------------------ plots -----------------------------------#

#lgbl  <- split(df1gbl,  df1gbl$setting)  # linear part, then detected
lgbl  <- split(df2gbl,  df2gbl$setting)  # detected, then linear part
#lgpix <- split(df1gpix, df1gpix$setting)
lgpix <- split(df2gpix, df2gpix$setting)

r2gbl <- lapply(lgbl, function(df) {
  df <- df[df$nlin > 2, ]
  tapply(df$r2, df$nlin, mean)
})
r2gpix <- lapply(lgpix, function(df) {
  df <- df[df$nlin > 2, ]
  tapply(df$r2, df$nlin, mean)
})

#ylim <- range(unlist(r2gbl, r2gpix))
ylim <- c(0.94, 1)
plot(NULL, xlim = c(1, 6), ylim = ylim, xlab = "settings",
     ylab = expression(R^2), main = "Fit of linear part")
for (i in 1:length(r2gbl)) {
  x1 <- r2gbl[[i]]
  x2 <- r2gpix[[i]]
  points(rep(i, length(x2)), x2, pch = (21:25)[as.integer(names(x2)) - 2],
         bg = "orchid", cex = 2)
  points(rep(i, length(x1)), x1, pch = (21:25)[as.integer(names(x1)) - 2],
         bg = "deepskyblue", cex = 2)
}
legend("bottomleft", pch = 21:24, legend = paste(3:6, "points in linear part"),
       bty = "n")
legend("bottomright", fill = c("deepskyblue", "orchid"), 
       legend = c("ArrayCAM", "Genepix"), bty = "n")


















