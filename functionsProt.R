
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
#  dat <- dat[iprot, c("ID", "Act_Result", "A1", "exposure_setting", "num_stdev_away")]
  dat <- dat[iprot, c("ID", "res", "block", "setting", "nstd")]
  dat <- dat[order(dat$ID),]
  dat$res <- log(dat$res)
  dlist <- split(dat, dat[c("ID", "setting")])
  df <- data.frame()
  for (i in 1:length(dlist)) {
    ilin <- linpart(dlist[[i]]$res, tol)  # index linear
    ilin1 <- rep(FALSE, 6)
    ilin1[ilin] <- TRUE
    idet <- dlist[[i]]$nstd >= nstd
    if (any(!idet)) 
      idet[which(!idet)[1]:length(idet)] <- FALSE  # points to the right of 1st F
    ikeep <- which(ilin1 & idet)
    if (length(ikeep) == 0) ikeep = 0
    dfrow <- data.frame(ID = dlist[[i]]$ID[1],
                        setting = dlist[[i]]$setting[1],
                        start = min(ikeep), end = max(ikeep))
    df <- rbind(df, dfrow)
  }
  return(df)
}

dfDetLin <- function(dat, nstd = 2, tol = 2.5) {
  iprot <- grep("PF|MAL", dat$ID)
#  dat <- dat[iprot, c("ID", "Act_Result", "A1", "exposure_setting", "num_stdev_away")]
  dat <- dat[iprot, c("ID", "res", "block", "setting", "nstd")]
  dat <- dat[order(dat$ID),]
  dat$res <- log(dat$res)
  dlist <- split(dat, dat[c("ID", "setting")])
  df <- data.frame()
  for (i in 1:length(dlist)) {
    idet <- dlist[[i]]$nstd >= nstd
    if (any(!idet)) 
      idet[which(!idet)[1]:length(idet)] <- FALSE  # points to the right of 1st F
    ilin <- linpart(dlist[[i]][idet, "res"], tol)  # index linear
    ikeep <- ifelse(ilin == 0, 0, which(idet)[ilin])
    dfrow <- data.frame(ID = dlist[[i]]$ID[1], 
                        setting = dlist[[i]]$setting[1],
                        start = min(ikeep), end = max(ikeep))
    df <- rbind(df, dfrow)
  }
  return(df)
}

# names(dat) = c("block", "ID", "setting", "res", "nstd.x", "nstd.y")
dfLinDet2 <- function(dat, nstd = 2, tol = 2.5) {  
  dat$res <- log(dat$res)
  dlist <- split(dat, dat[c("ID", "setting")])
  dfout <- data.frame()
  for (i in 1:length(dlist)) {
    df <- dlist[[i]][order(dlist[[i]]$block), ]
    ilin <- linpart(df$res, tol)  # index linear
    ilin1 <- rep(FALSE, 6)
    ilin1[ilin] <- TRUE
    idet <- df$nstd.x >= nstd & df$nstd.y >= nstd
    if (any(!idet)) 
      idet[which(!idet)[1]:length(idet)] <- FALSE  # points to the right of 1st F
    ikeep <- which(ilin1 & idet)
    if (length(ikeep) == 0) ikeep = 0
    dfrow <- data.frame(ID = df$ID[1], setting = df$setting[1],
                        start = min(ikeep), end = max(ikeep))
    dfout <- rbind(dfout, dfrow)
  }
  return(dfout)
}

dfLinDet2 <- function(dat, nstd = 2, tol = 2.5) {  
  ndil <- length(unique(dat$block))
  dat$res <- log(dat$res)
  dlist <- split(dat, dat[c("ID", "setting")])
  dfout <- data.frame()
  for (i in 1:length(dlist)) {
    df <- dlist[[i]][order(dlist[[i]]$block), ]
    ilin <- linpart(df$res, tol)  # index linear
    ilin1 <- rep(FALSE, ndil)
    ilin1[ilin] <- TRUE
    idet <- df$nstd.x >= nstd & df$nstd.y >= nstd
    if (any(!idet)) 
      idet[which(!idet)[1]:length(idet)] <- FALSE  # points to the right of 1st F
    ikeep <- which(ilin1 & idet)
    if (length(ikeep) == 0) ikeep = 0
    r2 <- ifelse(length(ikeep) > 2, 
                 summary(lm(df$res[ikeep] ~ (1:ndil)[ikeep]))$r.squared, NA)
    dfrow <- data.frame(ID = df$ID[1], setting = df$setting[1],
                        start = min(ikeep), end = max(ikeep), r2 = r2)
    dfout <- rbind(dfout, dfrow)
  }
  dfout$nlin <- with(dfout, end - start + 1)
  return(dfout)
}

# names(dat) = c("block", "ID", "setting", "res", "nstd.x", "nstd.y")
dfDet2Lin <- function(dat, nstd = 2, tol = 2.5) {
  dat$res <- log(dat$res)
  dlist <- split(dat, dat[c("ID", "setting")])
  dfout <- data.frame()
  for (i in 1:length(dlist)) {
    df <- dlist[[i]][order(dlist[[i]]$block), ]
    idet <- df$nstd.x >= nstd & df$nstd.y >= nstd
    if (any(!idet)) 
      idet[which(!idet)[1]:length(idet)] <- FALSE  # points to the right of 1st F
    ilin <- linpart(df[idet, "res"], tol)  # index linear
    ikeep <- ifelse(ilin == 0, 0, which(idet)[ilin])
    r2 <- ifelse(length(ikeep) > 2, 
                 summary(lm(df$res[ikeep] ~ (1:ndil)[ikeep]))$r.squared, NA)
    dfrow <- data.frame(ID = df$ID[1], setting = df$setting[1],
                        start = min(ikeep), end = max(ikeep), r2 = r2)
    dfout <- rbind(dfout, dfrow)
  }
  dfout$nlin <- with(dfout, end - start + 1)
  return(dfout)
}

dfNsdGBL <- function(fdir, ndblock = "C", match.dil = FALSE) {
  flist <- list.files(fdir)
  eset  <- as.integer(gsub(".csv", "", flist))  # exposure settings
  df <- data.frame()
  for (i in 1:length(flist)) {
    dat <- read.csv(paste(fdir, flist[i], sep = ""), as.is = TRUE)
    dat$res <- suppressWarnings(as.numeric(dat$Act_Result))
    if (match.dil) {
      for (j in 1:6) {
        block <- LETTERS[j]
        nodna <- dat[tolower(dat$ID) == "nodna" & grepl(block, dat$A1), "res"]
        mnd <- mean(nodna)
        snd <- sd(  nodna)
        prot <- dat[grepl("PF|MAL", dat$ID) & grepl(block, dat$A1), c("ID", "A1", "res")]
        prot <- aggregate(prot$res, by = list(gsub("[[:digit:]]", "", prot$A1),
                                              prot$ID), FUN = mean)
        names(prot) <- c("block", "ID", "res")
        prot$setting <- eset[i]   
        prot$nstd    <- (prot$res - mnd)/snd
        df <- rbind(df, prot)
      }
    } else {
      nodna <- dat[tolower(dat$ID) == "nodna" & grepl(ndblock, dat$A1), "res"]
      mnd <- mean(nodna)
      snd <- sd(nodna)
      prot <- dat[grepl("PF|MAL", dat$ID) & grepl("[A-F]", dat$A1), c("ID", "A1", "res")]
      prot <- aggregate(prot$res, by = list(gsub("[[:digit:]]", "", prot$A1),
                                            prot$ID), FUN = mean)
      names(prot) <- c("block", "ID", "res")
      prot$setting <- eset[i]   
      prot$nstd    <- (prot$res - mnd)/snd
      df <- rbind(df, prot)
    }
  }  
  return(df[order(df$ID, df$setting, df$block), ])
}

dfNsdGenepix <- function(fdir, ndblock = 3, match.dil = FALSE) {
  flist <- list.files(fdir)
  eset  <- as.integer(gsub(".csv", "", flist))  # exposure settings
  df <- data.frame()
  for (i in 1:length(flist)) {
    dat <- read.csv(paste(fdir, flist[i], sep = ""), as.is = TRUE, check.names = FALSE)
#    dat$res <- suppressWarnings(as.numeric(dat[, "F635 Mean - B635"]))
    dat$res <- suppressWarnings(as.numeric(dat[, "F635 Median - B635"]))  
    dat <- dat[c("Block", "ID", "res")]
    dat$Block <- (dat$Block + 1) %/% 2
    if (match.dil) {
      for (j in 1:6) {
        nodna <- dat[tolower(dat$ID) == "nodna" & dat$Block == j, "res"]
        mnd <- mean(nodna)
        snd <- sd(  nodna)
        prot <- dat[grepl("PF|MAL", dat$ID) & dat$Block == j, c("ID", "Block", "res")]
        prot <- aggregate(prot$res, by = list(prot$Block, prot$ID), FUN = mean)
        names(prot) <- c("block", "ID", "res")
        prot$setting <- eset[i]   
        prot$nstd    <- (prot$res - mnd)/snd
        df <- rbind(df, prot)
      }
    } else {
      nodna <- dat[tolower(dat$ID) == "nodna" & dat$Block == ndblock, "res"]
      mnd <- mean(nodna)
      snd <- sd(nodna)
      prot <- dat[grepl("PF|MAL", dat$ID) & dat$Block %in% 1:6, c("ID", "Block", "res")]
      prot <- aggregate(prot$res, by = list(prot$Block, prot$ID), FUN = mean)
      names(prot) <- c("block", "ID", "res")
      prot$setting <- eset[i]   
      prot$nstd    <- (prot$res - mnd)/snd
      df <- rbind(df, prot)
    }
  }  
  return(df[order(df$ID, df$setting, df$block), ])
}




