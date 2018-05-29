doTestGSEA <- function(bm, tgt, genesAll) {
    bms <- filter(bm, ensembl_gene_id %in% genesAll)
    bms <- filter(bms, go_id != "")
    bck <- genesAll[! genesAll %in% tgt]
    goIds <- unique(bms$go_id)
    bmsTgt <- filter(bms, ensembl_gene_id %in% tgt)
    bmsBck <- filter(bms, ensembl_gene_id %in% bck)
    print(paste0("Testing ", length(goIds), "GO terms"))
    res <- lapply(goIds, function(gid) {
        tab <- matrix(NA, nrow=2, ncol=2)
        rownames(tab) <- c("tgt", "bck")
        colnames(tab) <- c("inSet", "notInSet")
        tab["tgt", "inSet"] <- sum(bmsTgt$go_id %in% gid)
        tab["bck", "inSet"] <- sum(bmsBck$go_id %in% gid)
        tab["tgt", "notInSet"] <- nrow(bmsTgt)-tab["tgt", "inSet"]
        tab["bck", "notInSet"] <- nrow(bmsBck)-tab["bck", "inSet"]
        fish <- fisher.test(tab)
        res <- data.frame(
            goId=gid,
            pval=fish$p.value,
            odds=fish$estimate,
	    N=tab["tgt", "inSet"], 
	    n=tab["tgt", "notInSet"],
	    M=tab["bck", "inSet"],
	    m=tab["bck", "notInSet"],
	    stringsAsFactors=FALSE)
        return(res)
    })
    res <- bind_rows(res)
    res$pvalAdj <- p.adjust(res$pval, method="BH")
    res <- res %>% arrange(pval)
    ind <- match(res$goId, bms$go_id)
    res$desc <- bms$name_1006[ind]
    res$type <- bms$namespace_1003[ind]
    return(res)
}
