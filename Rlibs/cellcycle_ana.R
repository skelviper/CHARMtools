# functions for cellcycle analysis
library(tidyverse)
library(patchwork)
library(ggpubr)
library(cellAlign)
library(Seurat)
library(circlize)
library(ComplexHeatmap)
library(DescTools)

# functions for align two group of cell and compare cdps, adapt from cellAlign, NatMethod 2018
getCelltypeOrder <- function(obj,ct,sample_number = NULL){
    # order cellcycle using adapted Tanay's method(Nagano et al. 2017)
    # assuming cdps already save in a Seurat object slot named "cdps"
    # ct: celltype but actually is the expression, e.g. 'celltype == "mix late mesenchyme"'
    #     similarly 'stage == "E75"'
    
    mat <- obj[["cdps"]]@data %>% as.matrix() 
    if(is.null(sample_number)){
        tempOrderDF <- obj@meta.data %>% filter(cellcycle_threshold != "G0",cellcycle_threshold != "M") %>% filter(eval(parse(text=ct))) 
    }
    else{
        tempOrderDF <- obj@meta.data %>% filter(cellcycle_threshold != "G0",cellcycle_threshold != "M") %>% filter(eval(parse(text=ct)))%>% sample_n(sample_number)
    }
    
    #dim(tempOrderDF) %>% print()
    tempOrderDF <- tempOrderDF %>% mutate(clusterOrder = str_replace(paste0(str_replace(celltype," ","_"),"_",cluster,"_",sub_k_cluster),"_sub_",""))
    tempOrderDF$cellcycle_threshold <- factor(tempOrderDF$cellcycle_threshold,levels = c("G1","Early-S","Mid-S","Late-S","G2","M"))
    tempMat <- tempOrderDF %>% select(cellname,clusterOrder)%>% left_join(mat %>% t() %>% as.data.frame() %>% rownames_to_column("cellname")) %>% 
        group_by(clusterOrder) %>% select(-cellname) %>% summarise_all(mean) %>% na.omit() %>% ungroup() %>% column_to_rownames("clusterOrder")

    tempOrderDF <- tempOrderDF %>% ungroup() %>%  mutate(varrepli = var(raw_repli_score),varnearp = var(near_p),snearp=scale(near_p),sfar=scale(farAvg))

    tempOrderDF <- tempOrderDF %>% group_by(cellcycle_threshold) %>% 
                 mutate(order = ifelse(cellcycle_threshold%in% c("Early-S"),near_p /varnearp + raw_repli_score / varrepli,
                                    ifelse(cellcycle_threshold %in% c("Mid-S","G2"), near_p /varnearp,
                                      ifelse(cellcycle_threshold %in% c("Late-S"),near_p / varnearp - raw_repli_score / varrepli,
                                            ifelse(cellcycle_threshold %in% c("G1","G0"),snearp + sfar,
                                                  ifelse(cellcycle_threshold %in% c("M","Unknown"),maxbinorder)))))) %>% arrange(cellcycle_threshold,order)
    return(tempOrderDF %>% pull(cellname))
}

align_cellcycle <- function(obj,celltype1,celltype2,numPts = 50,sample_number = NULL){
    traj_celltype1 <- getCelltypeOrder(obj = obj,ct = celltype1,sample_number = sample_number)%>% as.data.frame() %>% mutate(order = row_number() / n()) %>% deframe()
    traj_celltype2 <- getCelltypeOrder(obj = obj,ct = celltype2,sample_number = sample_number)%>% as.data.frame() %>% mutate(order = row_number() / n()) %>% deframe()

    interGlobal_celltype1 = cellAlign::interWeights(expDataBatch = mat[,traj_celltype1%>% names()], trajCond = traj_celltype1,winSz = 0.1, numPts = numPts)
    interGlobal_celltype2 = cellAlign::interWeights(expDataBatch = mat[,traj_celltype2%>% names()], trajCond = traj_celltype2,winSz = 0.1, numPts = numPts)
    interScaledGlobal_celltype1 = cellAlign::scaleInterpolate(interGlobal_celltype1)
    interScaledGlobal_celltype2 = cellAlign::scaleInterpolate(interGlobal_celltype2)

    # here both local and global alignmnet can be use, but for whole cellcycle i dont think local alignment is necessary since 
    # every celltype has G1~G2 and it's contact decay profile look similar.

    #alignment = localAlign(interScaledGlobal_celltype1$scaledData,interScaledGlobal_celltype2$scaledData,threshPercent = Thresh)
    #
    #mapping = mapRealDataLocal(alignment,intTrajQuery = interScaledGlobal_celltype1$traj, realTrajQuery = traj_celltype1,
    #                            intTrajRef = interScaledGlobal_celltype2$traj, realTrajRef = traj_celltype2)

    alignment = globalAlign(interScaledGlobal_celltype1$scaledData,interScaledGlobal_celltype2$scaledData,scores = list(query = interScaledGlobalML$traj, 
                                                     ref = interScaledGlobalEM$traj),sigCalc = F, numPerm = 20)
    mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobal_celltype1$traj, realTrajQuery = traj_celltype1,
                                intTrajRef = interScaledGlobal_celltype2$traj, realTrajRef = traj_celltype2)

    #code for generate alignment plot
    ref_annotation <- mapping$refAssign %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("refAssign") %>% unnest(V1) %>% mutate(cellname = V1) %>% 
        select(-V1) %>% left_join(obj[[]] %>% select(cellname,cellcycle_threshold)) %>% group_by(refAssign)%>% 
        mutate(cellcycle_ref = cellcycle_threshold[which.max(n())],metaNodeRef=refAssign) %>% ungroup() %>% select(metaNodeRef,cellcycle_ref) %>% unique() 

    query_annotation <- mapping$queryAssign %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("queryAssign") %>% unnest(V1) %>% mutate(cellname = V1) %>% 
        select(-V1) %>% left_join(obj[[]] %>% select(cellname,cellcycle_threshold)) %>% group_by(queryAssign)%>% 
        mutate(cellcycle_query = cellcycle_threshold[which.max(n())],metaNodeQuery=queryAssign) %>% ungroup() %>% select(metaNodeQuery,cellcycle_query) %>% unique() 

    metaNodePt = mapping$metaNodesPt
    metaNodePt = metaNodePt[order(metaNodePt$ptQuery), ]
    metaNodePt$align = 1:nrow(metaNodePt)
    metaNodePtLong = reshape2::melt(metaNodePt[, c("ptQuery", "ptRef", "align")], id.vars = c("align"))
    metaNodePtLong = reshape2::melt(metaNodePt, id.vars = c("align", "metaNodeQuery", 
            "metaNodeRef"))

    lev2num <- function(str){
        return(which(c('G0','M','G1','Early-S','Mid-S','Late-S','G2') == str))
    }
    temp <- metaNodePtLong %>% filter(variable != "diff") %>% left_join(ref_annotation %>% group_by(cellcycle_ref) %>% mutate(rank_ref = row_number() / n())) %>%
                                                      left_join(query_annotation %>% group_by(cellcycle_query) %>% mutate(rank_query = row_number() / n())) %>% 
         rowwise() %>% mutate(shift = rank_ref + lev2num(cellcycle_ref) - rank_query - lev2num(cellcycle_query)) %>% 
        mutate(cellcycle = ifelse(variable == "ptQuery",as.character(cellcycle_query),as.character(cellcycle_ref))) %>% 
        mutate(variable = ifelse(variable == "ptQuery",str_extract(celltype1,pattern = '(?<=")[A-Za-z0-9 ]+(?=")'),str_extract(celltype2,pattern = '(?<=")[A-Za-z0-9 ]+(?=")')))

    # draw alignment plot
    temp$variable <- factor(temp$variable,levels = rev(c(str_extract(celltype1,pattern = '(?<=")[A-Za-z0-9 ]+(?=")'),str_extract(celltype2,pattern = '(?<=")[A-Za-z0-9 ]+(?=")'))))
    alignment_plot <- ggplot(temp %>% filter(variable != "diff"), aes(x = variable, y = value, group = align)) + 
            geom_line(color = "grey") + theme_bw() + geom_point(aes(color = cellcycle)) + 
            coord_flip() + ggtitle("meta-nodes alignment")  + theme_Publication()#+ scale_color_manual(values = cellcyclecolors)

    # draw histogram
    mean_num_old <- temp %>% rowwise() %>%  mutate(shift = lev2num(cellcycle_ref) - lev2num(cellcycle_query)) %>% pull(shift) %>% mean()

    mean_num <- temp %>% pull(shift) %>% mean()

    histogram_plot <- temp %>% gghistogram(x="shift",bins= 16) + theme_Publication() + xlab(paste0("shift (mean = ",round(mean_num,3),")")) + scale_x_continuous(limits = c(-2,2))
    #print(mean_num_old)
    #print(mean_num)
    p_return <- (alignment_plot | histogram_plot) + plot_layout(widths = c(2, 1))

    return(c(list(p_return),mean_num,list(temp)))
}

# plot cdps heatmap
plotCdpsHeatmap <- function(obj,ct){
    order = getCelltypeOrder(obj,ct = ct)
    tempOrderDF = obj[[]][order,]
    heatmap_mat <- mat[,order] %>% as.data.frame
    heatmap_mat[heatmap_mat > 0.025] <- 0.025

    col_fun = colorRamp2(c(-1, 0, 1), c("#ffffff", "#9c76c2", "#272d74"))
    col_fun1 = colorRamp2(c(-1, 0, 1), c("#ffffff","#7eb5b4","#058786"))

    Annodf <- tempOrderDF %>% select(repli_score,G1S.Score,G2M.Score,sub_k_cluster,cluster,cellname,cellcycle_threshold)
    Annodf <- Annodf %>% group_by(cellcycle_threshold) #%>% mutate(repli_score=mean(repli_score),G1S.Score=mean(G1S.Score),G2M.Score=mean(G2M.Score))
    topAnno <- HeatmapAnnotation(df=Annodf %>% column_to_rownames("cellname")%>% select(repli_score,G1S.Score,G2M.Score,cellcycle_threshold) ,
                                col = list(repli_score = col_fun,G1S.Score = col_fun1,G2M.Score = col_fun1,cellcycle_threshold=cellcyclecolors),show_annotation_name=FALSE)

    p <- heatmap_mat%>% as.matrix() %>%
        Heatmap(cluster_rows = FALSE,cluster_columns = FALSE,show_row_names = FALSE, show_column_names = FALSE,colors,
           #    column_split= tempOrderDF %>% select(cellcycle_threshold,clusterOrder),
                use_raster=TRUE,top_annotation= topAnno,column_title=str_extract(ct,pattern = '(?<=")[A-Za-z0-9 ]+(?=")'),
                column_dend_reorder = TRUE,
                   heatmap_legend_param = list(title = "contacts %"), width = 6)
    return(p)
}

# compartment strength calculation
# 1. Tanay's method and simple method

#example of load requirements
#pairsPaths <- paste0("../../HiC_clean3/clean3_pairs/",dir("../../HiC_clean3/clean3_pairs/"))
#resolution = 1000000
#read_tsv("../pileup_stage_lineage/majorgroup/processed/compartment/Emb.compartment.1m.cis.vecs.tsv") %>% mutate(type = ifelse(E1 > 0 ,"A","B")) %>% na.omit() %>% select(chrom,start,end,type) %>% write_tsv("emb.abcomp.1m.tsv",col_names = FALSE)
#defineAB <- read_table("./emb.abcomp.1m.tsv",col_names=FALSE)
#names(defineAB) <- c("chromosome","start","end","AB")

#example of use this funciton
#tanay_res <- comp_score(pairsPaths,defineAB,method = "tanay",threads = 200)

calcCompScore <- function(pairsPaths,defineAB,method = "tanay",threads = 100，bedtoolsPath = "~/miniconda3/envs/py3/bin/"){
    registerDoParallel(threads)
    result <- foreach(filepath=pairsPaths , .combine="rbind", .packages=c("tidyverse","bedtoolsr")) %dopar% {
        options(bedtools.path = bedtoolsPath)
        #c(filepath,calcCompartmentStrengthTanay(pairsPath = filepath,bulkAB = defineAB,resolution=1000000))
        testdata <- read_table(filepath,comment = "#",col_names = F,col_types = cols(X2 = col_character() ,X4 = col_character())) 
        if(dim(testdata)[1] < 100){
            return(NA)
        }
        testdata <- testdata %>% select(2,3,4,5) %>% mutate(pair_id = row_number())     
        names(testdata) <- c("chr1","pos1","chr2","pos2","pair_id")
        testdata <- testdata %>% mutate(pos1e = pos1 + 1,pos2e = pos2+1)  %>% filter(chr1 == chr2)
        #testData <- cbind(bt.intersect(a=testdata %>% select(chr1,pos1,pos1e),b=defineAB,wa = TRUE,wb=TRUE,loj=TRUE),
        #                 bt.intersect(a=testdata %>% select(chr2,pos2,pos2e),b=defineAB,wa = TRUE,wb=TRUE,loj=TRUE) ) %>% na.omit()
        left <-  testdata %>% select(chr1,pos1,pair_id) %>% mutate(pos11 = floor(pos1 / resolution)*resolution)
        names(left) <- c("chromosome","pos1","pair_id","start")
        left <- left %>% left_join(defineAB) %>% select(-end)

        right <-  testdata %>% select(chr1,pos2,pair_id) %>% mutate(pos21 = floor(pos2 / resolution)*resolution)
        names(right) <- c("chromosome","pos2","pair_id","start")
        right <- right %>% left_join(defineAB) %>% select(-end)
        testdata <- cbind(left,right)

        names(testdata) <- paste0("V",seq(1:10))
        testdata <- testdata #%>% filter(!(V7 - V2 < 2000000))
        testdata <- testdata %>% select(V5,V10)
        names(testdata) <- c("left","right")
        ABcount <- testdata %>% group_by(left,right) %>% summarise(count = n()) %>% filter(left != "." & right != ".") %>% arrange(left,right)
        #compartmentStrength = as.numeric(ABcount[1,3]) * as.numeric(ABcount[4,3]) / ((as.numeric(ABcount[2,3]) + as.numeric(ABcount[3,3]))**2)
        Oaa = as.numeric(ABcount[1,3])
        Oab = as.numeric(ABcount[2,3]) + as.numeric(ABcount[3,3])
        Obb = as.numeric(ABcount[4,3])
        
        if(method == "tanay"){
            T = Oaa + Oab + Obb
            Ascore = (2*Oaa + Oab)/T
            Bscore = (2*Obb + Oab)/T
            CompScore = log2(2*Ascore*Bscore*T/Oab)
            return(c(filepath,CompScore,Ascore,Bscore))
        }
        if(method == "simple"){
            CompScore = (Oaa+Obb)/(2*Oab)
            Ascore = Oaa/Oab
            Bscore = Obb/Oab
            return(c(filepath,CompScore,Ascore,Bscore))
        }     
    }
    result <- result%>% as.data.frame() %>% as_tibble()
    names(result) <- c("pairs","compScore","PA","PB")
    result$compScore <- as.numeric(result$compScore)
    result$PA <- as.numeric(result$PA)
    result$PB <- as.numeric(result$PB)
    return(result)
}

# 2. approximation by calculate gini index of scAB value
# scAB value of each cell can be provided as input:x
# return the gini index of scAB value, and A and B score

calcCompGini <- function(x, n = rep(1, length(x)), unbiased = TRUE, conf.level = NA, R = 1000, type = "bca", na.rm = FALSE){
    x <- as.numeric(x)
    x <- rep(x, n)
    if (na.rm) 
        x <- na.omit(x)
    if (any(is.na(x)) || any(x < 0)) 
        return(NA_real_)
    i.gini <- function(x, unbiased = TRUE) {
        n <- length(x)
        x <- sort(x)
        res <- 2 * sum(x * 1:n)/(n * sum(x)) - 1 - (1/n)
        if (unbiased) 
            res <- n/(n - 1) * res
        return(pmax(0, res))
    }
    if (is.na(conf.level)) {
        res <- i.gini(x, unbiased = unbiased)
    }
    else {
        boot.gini <- boot(x, function(x, d) i.gini(x[d], unbiased = unbiased), 
            R = R)
        ci <- boot.ci(boot.gini, conf = conf.level, type = type)
        res <- c(gini = boot.gini$t0, lwr.ci = ci[[4]][4], upr.ci = ci[[4]][5])
    }
    
    n=length(x)
    lcres <- Lc(x)
    Ascore = 1 - lcres$L[(length(lcres$L)/2+1):length(lcres$L)] %>% sum() / (3*length(lcres$L)/8)
    Bscore = 1 - lcres$L[1:(length(lcres$L)/2)] %>% sum() / (length(lcres$L)/8)
    return(c(res,Ascore,Bscore))
}

