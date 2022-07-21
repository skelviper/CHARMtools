# Code for analysis of compartment dynamics
library(ggrepel)

# functions for discover scAB relationship with RNA expression (normally marker gene have higher scAB )

# Usage:
# genes <- read_table("~/share/Data/public/ref_genome/mouse_ref/M23/raw_data/gene.bed",col_names = c("chr","start","end","id","name","strand")) %>% group_by(name) %>% slice(1)
# resolution =1000000
# Idents(hires) <- hires$celltype
# options(repr.plot.width = 5,repr.plot.height=5,repr.plot.res = 200)
# meso_ecto_data <- getCpGDiffVocanoPlotData(hires,c("neural ectoderm"),c("early mesoderm"),genes,1000000,slot_cpg = "scab")
# meso_ecto_plot <- plotCpGDiffVocanoPlot(meso_ecto_data,"ectoderm vs. mesoderm",text_x_Pos = 20,scale_factor = 1e5)
# meso_ecto_plot+ theme_Publication()

getCpGDiffVocanoPlotData <- function(hires,ident1,ident2,genes,resolution=100000,slot_cpg = "cpg",slot_RNA = "SCT"){
    # calculate cpg diff
    DefaultAssay(hires) <- slot_cpg
    cpg_diff <- FindMarkers(hires,`ident.1` = ident1,`ident.2` = ident2,slot = "data",logfc.threshold = 0,mean.fxn = rowMeans,test.use = "wilcox")
    cpg_diff <- cpg_diff %>% rownames_to_column("region") %>% 
    mutate(chr = str_extract(region,"chr[0-9]+"),pos = as.numeric(str_extract(region,"[0-9.]+e.[0-9]{2}"))) 
    names(cpg_diff)[3] <- "mean_diff"
    # calculate rna diff
    DefaultAssay(hires) <- slot_RNA
    diff_SCT <- FindMarkers(hires,`ident.1` = ident1,`ident.2` = ident2,only.pos = TRUE)
    RNAmarkers <- rownames(diff_SCT %>% arrange(p_val_adj) %>% head(100))
    
    plotData <- genes %>% mutate(pos = floor(((start+end)/2+resolution/2)/resolution)*resolution) %>% right_join(cpg_diff) %>% mutate(log_p_val= -log10(p_val))
    plotData <- plotData %>% full_join(tibble(name=RNAmarkers,type="Marker")) %>% filter(!is.na(chr))
    plotData$type[is.na(plotData$type)] <- "Other"
    
    return(plotData)
}
plotCpGDiffVocanoPlot <- function(plotData, plotTitle,text_x_Pos=0.04,name_show_num =15,scale_factor = 1){
    plotData <- plotData %>% mutate(mean_diff = scale_factor * mean_diff)
    plot <- ggplot(plotData %>% filter(type=="Other")) + geom_point(aes(x=mean_diff,y=log_p_val),color='grey',alpha=0.1) + 
    geom_point(data = plotData %>% filter(type=="Marker"),aes(x=mean_diff,y=log_p_val),shape=1) + 
    ggtitle(plotTitle)  + theme(plot.title = element_text(hjust = 0.5))  + 
    geom_vline(aes(xintercept=mean((plotData%>% filter(type=="Marker"))$mean_diff)),linetype="longdash") + 
    geom_vline(aes(xintercept=0))  + 
    annotate("text",label = paste0("Mean diff: ",as.character( round(mean((plotData%>% filter(type=="Marker"))$mean_diff),digits = 3))), x = text_x_Pos, y = 0,size = 4 ) + 
    geom_text_repel(data = plotData %>% filter(type == "Marker") %>% arrange(desc(log_p_val)) %>% head(name_show_num),
                 aes(x=mean_diff,y=log_p_val,label=name),max.overlaps = getOption("ggrepel.max.overlaps", default = 20))

    
    print(paste0("mean difference of CpG of all gene is ",as.character( mean((plotData%>% filter(type=="Marker"))$mean_diff))))
    plot
}

# functions for analysis pairs file
# Usage:
# bintable is a dataframe with chrom-start-end as the first 3 columns

# library(furrr)
# plan(multicore, workers = 100)
# options(future.rng.onMisuse="ignore")
# intermingle_result <- future_map(colnames(E75_hires),~calcInterRatio(paste0("../../HiC_clean3/clean3_pairs/",.x,".pairs.gz"),cluster_res) %>% select(4)) %>% suppressMessages()
# intermingle_result <- bind_cols(intermingle_result)
# colnames(intermingle_result) <- str_extract(colnames(intermingle_result),"(Gas|Org)[a-zA-Z0-9]+")

calcInterRatio <- function(file,bintable){
    library(valr)
    pairs <- readPairs(file)
    pairs <- pairs %>% mutate(type = ifelse(chrom1 == chrom2 ,"intra","inter"))

    left <- pairs %>% select(chrom1,pos1,type) %>% mutate(end = pos1 + 1) %>% select(chrom1,pos1,end,type)
    names(left) <- c("chrom","start","end","type")
    right <- pairs %>% select(chrom2,pos2,type) %>% mutate(end = pos2 + 1) %>% select(chrom2,pos2,end,type)
    names(right) <- c("chrom","start","end","type")

    # betoolsr is the R wrapper of bedtools and not suitable for parallel computation, use valr instead.
    #intermingle_res <- bt.intersect(bintable %>% select(1:3) ,rbind(left,right),wa = T,wb=T,loj=TRUE) %>%
    #     group_by(V1,V2,V3) %>% summarise(count = sum(V7 == "inter") / n()) %>% ungroup()
    #names(intermingle_res) <- c("chrom","start","end","interRatio")
    intermingle_res <- bed_intersect(bintable %>% select(1:3),rbind(left,right)) %>% 
        group_by(chrom,start.x,end.x) %>% summarise(count = sum(type.y == "inter") / n())
    names(intermingle_res) <- c("chrom","start","end",file)
    intermingle_res <- bintable %>% select(1:3) %>% left_join(intermingle_res)
    return(intermingle_res)
}

