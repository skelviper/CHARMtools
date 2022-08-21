# Functions for Difference detechtion

load_mat <- function(filepaths,bintable,threads = 40){
    registerDoParallel(threads)
    mat <- foreach(filepath = filepaths ,.combine = "cbind",.errorhandling = "remove") %dopar%{
        example <- fread(filepath,header = FALSE,tmpdir = "./",col.names = c("chrom","pos1","chrom2","pos2","distance"))
        #example <- example %>% mutate(band = pos2-pos1) %>% group_by(band) %>% 
            #mutate(distance_adj = scale(distance))
            #mutate(aveDistance = ave(distance), distance_adj = distance - aveDistanc)
            #mutate(aveDistance = ave(distance), distance_adj = distance)
        example <- bintable %>% left_join(example) %>% select(chrom,pos1,pos2,distance)
        example[is.na(example)] <- 10
        example <- cbind(example %>% head(dim(example)[1] / 2), example %>% tail(dim(example)[1] / 2) %>% select(distance)) 
        names(example)[5] <- "distance2"
        gc()
        example %>% select(distance,distance2)
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
    library(valr)
    #accept a bedpe file file and project it to 1d-genome ,return bed file
    left <- bedpe%>% select(1:3)
    right <- bedpe %>% select(4:6)
    names(left) <- c("chrom","start","end")
    names(right) <- c("chrom","start","end")
    all <- rbind(left,right) %>% valr::bed_sort() %>%  unique()
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

# functions for analysis d3D's results

# Usage:
#registerDoParallel(100)
# ecto <- foreach(cellname = celltype1 ,.combine = "cbind",.errorhandling = "stop") %dopar%{
#     print(cellname)
#     pairsPath <- paste0("../../HiC_clean3/clean3_pairs/",cellname,".pairs.gz")
#     count_contacts_in_region(pairsPath,differences %>% select(1:6)) %>% select(7)
# }
# meso <- foreach(cellname = celltype2 ,.combine = "cbind",.errorhandling = "stop") %dopar%{
#     pairsPath <- paste0("../../HiC_clean3/clean3_pairs/",cellname,".pairs.gz")
#     count_contacts_in_region(pairsPath,differences %>% select(1:6)) %>% select(7)
# }

count_contacts_in_region <- function(pairsPath,regions){
    bins <- regions %>% project1D() %>% mutate(id = row_number())
    names(bins) <- c("chrom","start","end","id")
    
    pairs <- readPairs(pairsPath) %>% filter(chrom1==chrom2) %>% 
        mutate(start1 = pos1,end1 = start1 + 1,start2 = pos2,end2 = start2 + 1) %>% 
        select(chrom1, start1,end1,chrom2,start2,end2) %>% mutate(pairID = row_number())
    
    left <- pairs %>% select(1:3,7)
    names(left) <- c("chrom","start","end","pairID")
    right <- pairs %>% select(4:6,7)
    names(right) <- c("chrom","start","end","pairID")
    left_bin <- valr::bed_intersect(left,bins) %>% select(pairID.x,chrom,start.y,end.y)
    right_bin <- valr::bed_intersect(right,bins)%>% select(pairID.x,chrom,start.y,end.y)
    overlap_bin_pairs <- inner_join(left_bin,right_bin,by = "pairID.x")
    names(overlap_bin_pairs) <- c("pairID","chrom1","start1","end1","chrom2","start2","end2")
    overlap_bin_pairs <- overlap_bin_pairs %>% group_by(chrom1,start1,end1,chrom2,start2,end2) %>% summarise(count = length(unique(pairID))) 
    res <- regions %>% left_join(overlap_bin_pairs)#  %>% na.omit()
    res[is.na(res)] <- 0
    
    return(res)
}