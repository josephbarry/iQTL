getCellFtrsThyroid <- function(cal, seg, scales=c(1, 10, 100)) {
    extractFeat <- function(ref, seg, type, scales) {
        r <- ref[, , 1]
        g <- ref[, , 2]
        b <- ref[, , 3]
        gm <- computeFeatures.shape(seg)
        mr <- computeFeatures.moment(seg, r)
        mg <- computeFeatures.moment(seg, g)
        mb <- computeFeatures.moment(seg, b)
        colnames(mr) <- paste0("r.", colnames(mr))
        colnames(mg) <- paste0("g.", colnames(mg))
        colnames(mb) <- paste0("b.", colnames(mb))
        br <- computeFeatures.basic(seg, r)
        bg <- computeFeatures.basic(seg, g)
        bb <- computeFeatures.basic(seg, b)
        colnames(br) <- paste0("r.", colnames(br))
        colnames(bg) <- paste0("g.", colnames(bg))
        colnames(bb) <- paste0("b.", colnames(bb))
        hr <- computeFeatures.haralick(seg, r, haralick.scales=scales)
        hg <- computeFeatures.haralick(seg, g, haralick.scales=scales)
        hb <- computeFeatures.haralick(seg, b, haralick.scales=scales)
        colnames(hr) <- paste0("r.", colnames(hr))
        colnames(hg) <- paste0("g.", colnames(hg))
        colnames(hb) <- paste0("b.", colnames(hb))
        ft <- cbind(gm, mr, mg, mb, br, bg, bb, hr, hg, hb)
	ft <- data.frame(spot=1, object=type, id=seq_len(nrow(ft)), ft, stringsAsFactors=FALSE)
    }
    ftPieces <- extractFeat(ref=cal, seg=seg, type="pieces", scales=scales)
    ftPieces$pieceID <- as.integer(rownames(ftPieces))
    return(ftPieces)
}


genPC2Tile <- function(img2, L=401, ft.mean, ft.sd, loadings) {
    dim0 <- dim(img2)
    seq.x <- seq(1, dim0[1], by=L)
    seq.y <- seq(1, dim0[2], by=L)
    L.x <- length(seq.x)
    L.y <- length(seq.y)
    img2.d <- imageData(img2)
    dummySeg <- array(1, dim=c(L, L))
    tiles <- array(NA, dim=c(L.x-1, L.y-1))
    dimnames(tiles)<- list(seq_len(L.x-1), seq_len(L.y-1))
    for (i in seq_len(L.x-1))
        for (j in seq_len(L.y-1)) {
    	r.x <- seq.x[i]:(seq.x[i+1]-1)
    	r.y <- seq.y[j]:(seq.y[j+1]-1)
	ref0 <- img2.d[r.x, r.y,]
    	fts.tile <- getCellFtrsThyroid(cal=ref0, seg=dummySeg)
        fts.transformed <- sapply(names(ft.mean), function(ft) {
            (log2(fts.tile[, ft])-ft.mean[ft])/ft.sd[ft]        
        })
	names(fts.transformed) <- names(ft.mean)
	fts.transformed[is.na(fts.transformed)] <- 0
	fts.transformed[is.infinite(fts.transformed)] <- 0
	pc2 <- 0
	for (ft in names(fts.transformed)) pc2 <- pc2+fts.transformed[ft]*loadings[ft]
	tiles[i, j] <- pc2    
    }
    img.tile <- as.Image(tiles)
    return(img.tile)
}
