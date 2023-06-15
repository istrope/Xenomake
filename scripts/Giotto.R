library(Giotto)
library(data.table)
library(dpyr)
library(ggplot2)
library(fastcluster)

initialize_gobject <- function(spatial_metadata,raw_expression){
    data_path = '~/Documents/GV_Project/final_processing/'
    workdir = '~/Documents/GV_Project/updated_matrices/corrected/Giotto_QC'
    myinst=createGiottoInstructions(save_plot = T,
                                    show_plot = T,
                                    return_plot = F,
                                    save_dir = workdir,
                                    python_path = '~/anaconda3/envs/giotto_env/bin/python')
    
    #load in count matrix and format
    expr <- data.table::fread(raw_expression,sep='\t',header=T) %>% as.data.frame()
    rownames(expr) <- expr$GENE
    expr <- expr[,2:ncol(expr)]
    #load spatial metadata
    spatial_results <- read.table(spatial_metadata,sep=',')

    colnames(spatial_results) <- c("barcode", "in_tissue", "array_row", "array_col", "col_pxl", "row_pxl")
    spatial_results$barcode <- gsub('-1','',spatial_results$barcode)
    spatial_results <- spatial_results[match(colnames(expr),spatial_results$barcode),]

    #create Giotto Object
    gobject  <- createGiottoObject(raw_exprs = expr,
                                    spatial_locs = spatial_locs,
                                    cell_metadata = spatial_results[,c('in_tissue','array_row','array_col')],
                                    instructions = myinst)
    return(gobject)
}

process_gobject <- function(gobject){
    metadata <- pDataDT(gobject)
    #Cluster Raw Expression Matrix to predict "in_tissue spots"
    #groups <- kmeans(gobject@raw_exprs,centers = 2,iter.max=10,nstart = 5)
    hc <- fastcluster::hclust(dist(gobject@raw_exprs))
    in_tissue <- 

    gobject <- filterGiotto(gobject = gobject,
                                     expression_threshold = 1,
                                     gene_det_in_min_cells = 25,
                                     min_det_genes_per_cell = 250,
                                     expression_values = c('raw'),
                                     verbose = T)
    gobject <- normalizeGiotto(gobject=gobject)
    gobject <- addStatistics(gobject = gobject)
    gobject <- calculateHVG(gobject=gobject,nr_expression_groups = 30, expression_values = c('normalized'))
    featgenes = gene_metadata[hvg == 'yes' & perc_cells >3 & mean_expr_det > 0.2]$gene_ID

    gobject <- runPCA(gobject=gobject ,genes_to_use=featgenes,scale_unit=F,center=T,method='irlba')
    gobject <- runUMAP(gobject,dimensions_to_use = 1:30)

    gobject <- createNearestNetwork(gobject = gobject,dimensions_to_use = 1:30,k = 15)
    gobject <- doLeidenCluster(gobject = gobject, resolution = 0.4, n_iterations = 200)

    return(gobject)
}

qc <- function(gobject){
    nr_genes <- spatPlot(gobject=gobject_new ,cell_color='nr_genes',color_as_factor=F,point_size=2,save_param=c(save_name='3-spatplot'),return_plot =T)
    all_counts <- as.data.frame(unlist(as.list(gobject_new@norm_expr)))
    colnames(all_counts) <- c('counts')
    count_histogram <- ggplot(all_counts,aes(counts))+ geom_histogram(binwidth = 0.1)+xlab('Normalized Count') + ylab('Count') + ggtitle('gobject_spacemake') + theme(plot.title = element_text(hjust = 0.5))

        #Plot Filter Combinations
    filter_combinations <- filterCombinations(gobject,expression_values = 'raw',expression_thresholds = 1,gene_det_in_min_cells = c(25,50,100),min_det_genes_per_cell = c(250,500,1000),return_plot = T)
        #Plot Cell Level
    cell_metadata <- gobject_new@cell_metadata
    cell_metadata$sample <- {snakemake.params}
    p1 <- ggplot(cell_metadata,aes(x=sample,y=nr_genes,color = sample)) + geom_boxplot() + ggtitle('number of Genes per Spot') + xlab(NULL) + ylab('nr_genes')+theme_classic()+theme(plot.title = element_text(hjust=0.5))

        #Plot Gene Level QC
    gene_metadata <- gobject_new@gene_metadata
    gene_metadata$sample <- {snakemake.params}
    p2 <- ggplot(gene_metadata,aes(x=sample,y=nr_cells,color=sample))+ geom_boxplot() + ggtitle('Number of Spots per Gene') + xlab(NULL) + ylab('nr_cell')+theme_classic()+theme(plot.title = element_text(hjust=0.5))

    library(ggpubr)
    p3 <- p1 + p2
    p4 <- nr_genes + filter_combinations
    ggarrange(p3,p4,count_histogram,ncol=1,nrow = 3) -> gobject_qc

    ggsave({snakemake.output},plot=gobject_qc)
}

spatial_cluster <- function(gobject){
    gobject <- readRDS({snakemake.input})
        num_cluster <- length(unique(gobject@cell_metadata$leiden_clus))

        #Create Spatial Network
        gobject <- createSpatialGrid(gobject = gobject , sdimx_stepsize = 400, sdimy_stepsize = 400, minimum_padding = 0)
        gobject <- createSpatialNetwork(gobject = gobject , method = 'kNN',k = 5, maximum_distance_knn = 400, name = 'spatial_network')

        #Find Spatially Variable Genes using rank binarization
        ranktest = binSpect(gobject,bin_method = 'rank',calc_hub = T,hub_min_int = 5,spatial_network_name = 'spatial_network')
        ## silhouette
        spatial_genes=silhouetteRank(gobject,expression_values="normalized")
        ext_spatial_genes = spatial_genes[1:1100,]$gene

        spat_cor_netw_DT = detectSpatialCorGenes(gobject ,
                                         method = 'network',
                                         spatial_network_name = 'spatial_network',
                                         subset_genes = ext_spatial_genes,
                                         network_smoothing=NULL,
                                         min_cells_per_grid = 4,
                                         cor_method = c('pearson'))
        # cluster spatial genes from correlation in previous line
        spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT,name = 'spat_clus',k = 4)

        sample_rate=2
        target=500
        tot=0
        gene_list = list()
        clust = spat_cor_netw_DT$cor_clusters$spat_clus
        for(i in seq(1, num_cluster)){gene_list[[i]] = colnames(t(clust[which(clust==i)]))}
        for(i in seq(1, num_cluster)){
            num_g=length(gene_list[[i]])
            tot = tot+num_g/(num_g^(1/sample_rate))
            }
        factor=target/tot
        num_sample=c()
        for(i in seq(1, num_cluster)){
            num_g=length(gene_list[[i]])
            num_sample[i] = round(num_g/(num_g^(1/sample_rate)) * factor)
            }
        set.seed(10)
        samples=list()
        union_genes = c()
        for(i in seq(1, num_cluster)){
            if(length(gene_list[[i]])<num_sample[i]){
            samples[[i]] = gene_list[[i]]
            }else{
                samples[[i]] = sample(gene_list[[i]], num_sample[i])
                }
            union_genes = union(union_genes, samples[[i]])
                }
        union_genes = unique(union_genes)

        #Run Spatially Aware Clustering Method "HMRF"
        HMRF_spatial_genes = doHMRF(gobject = gobject , expression_values = 'scaled', spatial_genes = union_genes, k = num_genes,
                                    spatial_network_name="spatial_network", betas = c(0, 5, 7),
                                    output_folder = {snakemake.output.hmrf})
        #Add HMRF output to giotto gobject
        gobject = addHMRF(gobject = gobject , HMRFoutput = HMRF_spatial_genes,
                        k = num_genes, betas_to_add = c(0, 5,10,15,20),
                        hmrf_name = 'HMRF')

        saveRDS(gobject,{snakemake.output})
}

find_markers <- function(gobject){
    mast_markers = findMarkers_one_vs_all(gobject = ac50, 
                                        method = 'mast',
                                        expression_values = 'scaled',
                                        cluster_column = 'leiden_clus')
    write.table(mast_markers,file={snakemake.output},sep='\t',row.names=F)
}