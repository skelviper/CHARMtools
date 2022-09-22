library(foreach)
library(doParallel)
library(data.table)
# Functions for Difference detechtion

load_mat <- function(filepaths,bintable,threads = 40,type="center"){
    registerDoParallel(threads)
    names(bintable) <- c("chrom","pos1","chrom2","pos2")
    mat <- foreach(filepath = filepaths ,.combine = "cbind",.errorhandling = "stop") %dopar%{
        example <- fread(filepath,header = FALSE,tmpdir = "./",col.names = c("chrom","pos1","chrom2","pos2","distance"))

        if(type == "center"){
            example <- example %>% mutate(band = pos2-pos1) %>% group_by(chrom,band) %>% 
            #mutate(aveDistance = mean(distance,na.rm = T), distance_adj = distance - aveDistance)
            mutate(distance_adj = as.vector(scale(distance)))
            example <- bintable %>% left_join(example) %>% select(chrom, 
                pos1, pos2, distance_adj)
            #example[is.na(example)] <- 0
            example <- cbind(example %>% head(dim(example)[1]/2), 
                example %>% tail(dim(example)[1]/2) %>% select(distance_adj))
            names(example)[5] <- "distance2"
            gc()
            example %>% select(distance_adj, distance2)
        }
        else{
            example <- bintable %>% left_join(example) %>% select(chrom,pos1,pos2,distance)
            #example[is.na(example)] <- 10
            example <- cbind(example %>% head(dim(example)[1] / 2), example %>% tail(dim(example)[1] / 2) %>% select(distance)) 
            names(example)[5] <- "distance2"
            gc()
            example %>% select(distance,distance2)
        }
        
    }
    return(as.matrix(mat))
}


convert_mat_to_contacts <- function(mat,type="center"){
    if(type == "center"){
        mat[mat > 0] <- 0
        mat[mat < 0] <- 1
    }
    else{
        mat[mat <= 3] <- 1
        mat[mat > 3] <- 0
    }

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

ttest <- function(...) {
    # R original test cant handle data exactily same.
   obj<-try(t.test(...), silent=TRUE)
   #if (is(obj, "try-error")) return(c(0,1)) else return(as.numeric(c(obj$statistic,obj$p.value)))
    if (is(obj, "try-error")) return(1) else return(obj$p.value)
}

d3Dtest <- function(x,y,method = "t"){
    # wrapper for test on single loci
    # supported methods: chi-square, t, fisher.test

    if(method == "t"){
        p <- ttest(x,y)
    }
    else{
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
    }

    return(p)
}

d3D <- function(mat1,mat2,binnames = rownames(mat1),threads = 200,p.adj.method = "BH",fdr_thres = 0.05,test.method = "t",resolution = 20000){
    library(furrr)
    plan(multicore, workers = threads)
    options(future.globals.maxSize = 320000000000)
    options(scipen = 99999)

    diff_raw <- future_map(seq(ncol(mat1)), function(x) d3Dtest(mat1[,x], mat2[,x],method = test.method))

    diff_format <- cbind(binnames,diff_raw %>% unlist() %>% as.numeric() %>% as.data.frame())
    names(diff_format) <- c("pos","pv")
    diff_format <- diff_format %>% mutate(FDR = p.adjust(pv,method = p.adj.method))
    #diff <- colMeans(mat1,na.rm=TRUE) - colMeans(mat2,na.rm=TRUE)
    ave_celltype1 <- colMeans(mat1,na.rm=TRUE)
    ave_celltype2 <- colMeans(mat2,na.rm=TRUE)
    diff <- ave_celltype1 - ave_celltype2
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

# BAD FUNCTION!!!!!!!!!!!!!
# Unknown bug cause output to be wrong
# cluster_difference <- function(differences,resolution = 20000,flank_size = 1){
#     differences_cluster <- differences %>% select(1:6) %>% bt.cluster(d = flank_size*resolution) %>% select(V4,V5,V6,everything()) %>%
#         bt.sort() %>% bt.cluster(d=flank_size*resolution ) %>% group_by(V7,V8) %>% mutate(cluster= cur_group_id()) %>% 
#         ungroup() %>% select(4:6,1:3,cluster)
#     names(differences_cluster) <- c(names(differences %>% select(1:6)),"cluster")
#     differences_cluster <- cbind(differences_cluster,differences %>% select(-c(1:6)) )        
    
#     return(differences_cluster)
# }

# cluster bedpe file using valr based on the first 6 colunms
cluster_bedpe <- function(cluster_temp,flank = 20000){
    names(cluster_temp)[1:3] <- c("chrom","start","end")
    cluster_temp <- cluster_temp %>% valr::bed_cluster(max_dist = flank)
    names(cluster_temp)[1:6] <- c("chrom1","start1","end1","chrom","start","end")
    names(cluster_temp)[length(names(cluster_temp))] <- "idleft"
    cluster_temp <- cluster_temp %>% valr::bed_cluster(max_dist = flank)
    names(cluster_temp)[4:6] <- c("chrom2","start2","end2")
    names(cluster_temp)[length(names(cluster_temp))] <- "idright"
    cluster_temp <- cluster_temp %>% group_by(idleft,idright) %>% mutate(cluster_id = cur_group_id())%>% ungroup() %>% select(-idleft,-idright) 
    return(cluster_temp)
}



flank_bedpe <- function(bedpe, flank = 20000){
    temp <- bedpe %>% select(1:6)
    names(temp) <- c("chrom1","start1","end1","chrom2","start2","end2")
    temp <- temp %>% mutate(start1 = start1 - flank, end1 = end1 + flank, start2 = start2 - flank, end2 = end2 + flank)
    #bedpe <- bedpe %>% mutate(start = start - flank, end = end + flank)
    return(cbind(temp,bedpe %>% select(-c(1:6))))
}