# Functions for Difference detechtion

load_mat <- function(filepaths,threads = 40){
    registerDoParallel(threads)
    mat <- foreach(filepath = filepaths ,.combine = "cbind",.errorhandling = "remove") %dopar%{
        example <- fread(filepath,header = FALSE,tmpdir = "./",col.names = c("chrom","pos1","chrom2","pos2","distance"))
        example <- example %>% mutate(band = pos2-pos1) %>% group_by(band) %>% 
            #mutate(distance_adj = scale(distance))
            #mutate(aveDistance = ave(distance), distance_adj = distance - aveDistanc)
            mutate(aveDistance = ave(distance), distance_adj = distance)
        example <- bintable %>% left_join(example) %>% select(chrom,pos1,pos2,distance_adj)
        example[is.na(example)] <- 10
        example <- cbind(example %>% head(dim(example)[1] / 2), example %>% tail(dim(example)[1] / 2) %>% select(distance_adj)) 
        names(example)[5] <- "distance_adj2"
        gc()
        example %>% select(distance_adj,distance_adj2)
    }
    return(as.matrix(mat))
}

convert_mat_to_contacts <- function(mat){
    mat[mat <= 3] <- 1
    mat[mat > 3] <- 0
    mode(mat) <- "integer"
    return(mat)
}

project1D <- function(bedpe){
    #accept a bedpe file file and project it to 1d-genome ,return bed file
    left <- bedpe%>% select(1:3)
    right <- bedpe %>% select(4:6)
    names(left) <- c("chrom","start","end")
    names(right) <- c("chrom","start","end")
    all <- rbind(left,right) %>%bt.sort() %>%  unique()
    return(all)
}

d3Dtest <- function(x,y,method = "chi-square"){
    # wrapper for test on single loci
    dfy <- as.data.frame(table(y))
    names(dfy) <- c("type","Freq")
    dfx <- as.data.frame(table(x))
    names(dfx) <- c("type","Freq")

    mat <- data.frame(type = factor(c(0,1))) %>% left_join(dfx,by="type") %>% left_join(dfy,by="type")
    mat <- mat[,2:3]
    mat[is.na(mat)] <- 0
    if (method == "chi-square"){
        p <- chisq.test(mat)$p.value %>% suppressWarnings()
    }
    else{
        p <- fisher.test(mat)$p.value %>% suppressWarnings()
    }
    return(p)
}

d3D <- function(mat1,mat2,binnames = rownames(mat1),threads = 200,p.adj.method = "BH",fdr_thres = 0.05,test.method = "chi-square",resolution = 20000){
    library(furrr)
    plan(multicore, workers = threads)
    options(future.globals.maxSize = 320000000000)

    diff_raw <- future_map(seq(ncol(mat1)), function(x) d3Dtest(mat1[,x], mat2[,x],method = test.method))

    diff_format <- cbind(binnames,diff_raw %>% unlist() %>% as.numeric() %>% as.data.frame())
    names(diff_format) <- c("pos","pv")
    diff_format <- diff_format %>% mutate(FDR = p.adjust(pv,method = p.adj.method))
    diff <- (colSums(mat1) / dim(mat1)[1]) - (colSums(mat2) / dim(mat2)[1])
    diff_format <- cbind(diff_format,diff)
    sig <- diff_format %>% filter(FDR < fdr_thres) %>% separate(pos, into = c("chrom1","pos1","pos2")) %>% 
        mutate(start1 = as.numeric(pos1) - resolution /2,start2 = as.numeric(pos2) - resolution / 2,
                end1 = start1 + resolution,end2 = start2 + resolution,chrom2 = chrom1) %>% 
        select(chrom1,start1,end1,chrom2,start2,end2,everything())
    sig <- sig %>% select(-pos1,-pos2)
    return(sig)
}