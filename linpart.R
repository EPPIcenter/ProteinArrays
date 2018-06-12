

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
#                              protein ranking                                 #                             #                        layout for plots per protein                          #
#------------------------------------------------------------------------------#
#==============================================================================#

#------------------------------------------------------------------------------#
#                                  functions                                   #
#------------------------------------------------------------------------------#

linpart <- function(means, tol = 2.5) {
  if (length(means) < 2) return(0:0)
  if (any(!is.finite(means))) {
    cutout <- which(!is.finite(means))[1] - 1
    if (cutout < 2) return(0:0)
    means <- means[1:cutout]
  }
  sl <- diff(means)
  min.sl <- min(sl)
  if (min.sl > 0) return(0)
  ratio <- min.sl/sl
  ipass <- which(0 < ratio & ratio < tol)
  return(min(ipass):(max(ipass) + 1))
}

dfLinDet <- function(dat, nstd = 2, tol = 2.5) {
  iprot <- grep("PF|MAL", dat$ID)
  dat <- dat[iprot, c("ID", "Act_Result", "A1", "exposure_setting", "num_stdev_away")]
  dat <- dat[order(dat$ID),]
  dat$res <- log(dat$Act_Result)
  dlist <- split(dat, dat[c("ID", "exposure_setting")])
  df <- data.frame()
  for (i in 1:length(dlist)) {
    ilin <- linpart(dlist[[i]]$res, tol)  # index linear
    ilin1 <- rep(FALSE, 6)
    ilin1[ilin] <- TRUE
    idet <- dlist[[i]]$num_stdev_away >= nstd
    ikeep <- which(ilin1 & idet)
    if (length(ikeep) == 0) ikeep = 0
    dfrow <- data.frame(ID = dlist[[i]]$ID[1],
                        setting = dlist[[i]]$exposure_setting[1],
                        start = min(ikeep), end = max(ikeep))
    df <- rbind(df, dfrow)
  }
  return(df)
}

dfDetLin <- function(dat, nstd = 2, tol = 2.5) {
  iprot <- grep("PF|MAL", dat$ID)
  dat <- dat[iprot, c("ID", "Act_Result", "A1", "exposure_setting", "num_stdev_away")]
  dat <- dat[order(dat$ID),]
  dat$res <- log(dat$Act_Result)
  dlist <- split(dat, dat[c("ID", "exposure_setting")])
  df <- data.frame()
  for (i in 1:length(dlist)) {
    idet <- dlist[[i]]$num_stdev_away >= nstd
    if (any(!idet)) 
      idet[which(!idet)[1]:length(idet)] <- FALSE  # points to the right of 1st F
    ilin <- linpart(dlist[[i]][idet, "res"], tol)  # index linear
    ikeep <- ifelse(ilin == 0, 0, which(idet)[ilin])
    dfrow <- data.frame(ID = dlist[[i]]$ID[1], 
                        setting = dlist[[i]]$exposure_setting[1],
                        start = min(ikeep), end = max(ikeep))
    df <- rbind(df, dfrow)
  }
  return(df)
}

#load("~/Downloads/gbl_dataframes3.RData")  # df_inna_gbl, linpart_gbldf, 
#load("~/Downloads/genepix_dataframes3.RData")  # df_inna2, df_inna
#save(linpart_gbldf, linpart_genepixdf, file = "~/Downloads/cams.RData")

#------------------------------------------------------------------------------#
#                                     Data                                     #
#------------------------------------------------------------------------------#

library(RColorBrewer)
cols <- c(brewer.pal(9, "Set1")[-6], brewer.pal(8, "Dark2"))
#cols <- rainbow(20)  # not as nice but shows that settings are not categorical
pardef <- par(no.readonly = TRUE)

bwid <- 0.4
nset <- 6    # number of settings in a camera
ndil <- 6    # number of dilutions
nprot <- 250

load("~/Downloads/cams.RData")  # linpart_gbldf, linpart_genepixdf
df1gbl     <- dfLinDet(linpart_gbldf,     nstd = 2, tol = 2.5)
df2gbl     <- dfDetLin(linpart_gbldf,     nstd = 2, tol = 2.5)
df1genepix <- dfLinDet(linpart_genepixdf, nstd = 2, tol = 2.5)
df2genepix <- dfDetLin(linpart_genepixdf, nstd = 2, tol = 2.5)

dilutions <- c(50, 200, 800, 3200, 12800, 51200)
settings1 <- sort(unique(df1gbl$setting))
settings2 <- sort(unique(df1genepix$setting))

dat1 <- linpart_gbldf
dat1 <- dat1[c("ID", "Act_Result", "A1", "exposure_setting", "num_stdev_away")]
dat1 <- dat1[order(dat1$ID),]
dat1$res <- log(dat1$Act_Result)
dlist1 <- split(dat1, dat1$ID)

dat2 <- linpart_genepixdf
dat2 <- dat2[c("ID", "Act_Result", "A1", "exposure_setting", "num_stdev_away")]
dat2 <- dat2[order(dat2$ID),]
dat2$res <- log(dat2$Act_Result)
dlist2 <- split(dat2, dat2$ID)

#------------------------------------------------------------------------------#
#                       Plots: linpart, then detection                         #
#------------------------------------------------------------------------------#

laymat <- matrix(c(1:12, 0, rep(13, 2), rep(14, 2), 0), 3, byrow = TRUE)
layout(laymat)
nstd <- 2
tol  <- 2.5

plist1 <- split(df1gbl,     df1gbl$ID)
plist2 <- split(df1genepix, df1genepix$ID)
sumord <- sapply(plist1, function(df) sum(df$end - df$start))
maxord <- sapply(plist1, function(df) max(df$end - df$start))  
ord    <- order(sumord, maxord, decreasing = TRUE)  # max first, then sum
prots <- names(plist1[ord])

par(ask = TRUE)
for (prot in prots) {
  cam1 <- split(dlist1[[prot]], dlist1[[prot]]$exposure_setting)
  cam2 <- split(dlist2[[prot]], dlist2[[prot]]$exposure_setting)
  cam1 <- cam1[order(as.numeric(names(cam1)))]
  cam2 <- cam2[order(as.numeric(names(cam2)))]
  for (j in 1:length(settings1)) {
    par(mar = c(0.5, 2, 0.5, 0))
    plot(1:6, cam1[[j]]$res, type = "b", pch = 19, xaxt = "n",
         col = "deepskyblue4", xlab = "dilution", ylab = "", lwd = 1.5)
    ilin <- linpart(cam1[[j]][1:6, "res"], tol)  # index linear
    lines((1:6)[ilin], cam1[[j]][ilin, "res"], col = "firebrick3", lwd = 1.5)
    idet <- cam1[[j]]$num_stdev_away >= nstd
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
    idet <- cam2[[j]]$num_stdev_away >= nstd
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

plist1 <- split(df2gbl,     df2gbl$ID)
plist2 <- split(df2genepix, df2genepix$ID)
sumord <- sapply(plist1, function(df) sum(df$end - df$start))
maxord <- sapply(plist1, function(df) max(df$end - df$start))  
ord    <- order(sumord, maxord, decreasing = TRUE)  # max first, then sum
prots <- names(plist1[ord])

par(ask = TRUE)
for (prot in prots) {
  cam1 <- split(dlist1[[prot]], dlist1[[prot]]$exposure_setting)
  cam2 <- split(dlist2[[prot]], dlist2[[prot]]$exposure_setting)
  cam1 <- cam1[order(as.numeric(names(cam1)))]
  cam2 <- cam2[order(as.numeric(names(cam2)))]
  for (j in 1:length(settings1)) {
    par(mar = c(0.5, 2, 0.5, 0))
    plot(1:6, cam1[[j]]$res, type = "b", pch = 19, xaxt = "n",
         col = "deepskyblue4", xlab = "dilution", ylab = "", lwd = 1.5)
    idet <- cam1[[j]]$num_stdev_away >= nstd
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
    idet <- cam2[[j]]$num_stdev_away >= nstd
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

#------------------------------------------------------------------------------#
#                         Plots: summary by camera                             #
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

#slist <- split(df1genepix, df1genepix$setting)  # linpart, then detection
slist <- split(df2genepix, df2genepix$setting)  # detection, then linpart
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

