library(getopt)
command= matrix(c(
  	'help'   		,"h", 0, "logical","No parameters: help document",
  	'count'  		,"c", 2, "character","Optional parameter: gene expression matrix, row are genes, column are cells",
  	'project_name'		,"a", 2, "character","Optional parameter: Seurat object name, if this parameter is defined, all operations are for this object; default: project",
  	'tax'   		,"t", 2, "character","Optional parameter: species information; default: human; human, mouse and rat can directly calculate the proportion of mitochondria, ribosomes, and red blood cells of each cell by specifying tax and filter them",
  	'nFeature_lowlimit'   	,"n", 2, "integer","Optional parameter: one of the cell filtration criteria: the lower limit of the number of genes expressed in each cell; the default is 200",
  	'nFeature_highlimit'   	,"N", 2, "integer","Optional parameter: one of the cell filtration criteria: the upper limit of the number of genes expressed in each cell: default 10000",
  	'mito_cutoff'		,"m", 2, "double","Optional parameter: one of the cell filtration criteria: the threshold of mitochondrial gene expression content of each cell: default 10",
  	'redcell_cutoff'	,"e", 2, "double","Optional parameter: one of the cell filtration criteria: the upper limit of the expression content of red blood cell marker genes in each cell: default 10",
	'double_file'		,"d", 2,"character","Optional parameter: cell doublets  prediction results",
	'sample_list'		,"l", 1,"character","Required parameter: order of samples; input format A, B, C",
	'group_list'		,"2", 2,"character","Optional parameter: sample grouping; input format case, control, case",
	'step'			,"s", 2,"character","Optional parameter: operation steps; default 1, 2, 3, 4 (1,filter cells; 2,Seurat integration data, PCA dimensionality reduction; 3,Seurat clustering identification of marker genes; 4,cell cycle evaluation)",
	'outdir'                ,"o", 2, "character","Optional parameter: output directory: default current directory",
	'batch'               	,"b", 2, "logical","Optional parameter: whether to re-calibrate the batch: Default: TRUE",
        'project.rdata'         ,"R", 2, "character","Optional parameter: Seurat output result .rdata file",
        'pcs_number'            ,"M", 2, "integer","Optional parameter: PCA dimensionality reduction: default 30",
        'resolution'            ,"r", 2, "character","Optional parameter: input format: 0.2, 0.4, 0.6 default 0.9",
        'top_num'               ,"i", 2, "integer","Optional parameter: draw heat map, scatter map and violin map. Number of genes: Default: 10",
        'cluster_number'        ,"I", 2, "character","Optional parameter: cluster number for re-clustering: input format: 1, 2, 3",
        'barcode_file'          ,"F", 2, "character","Optional parameter: barcode ID for re-clustering: input file: barcode ID number in the first column",
        'onlyplot'          	,"p", 2, "logical","Optional parameter: only draw gene scatter plots, violin plots and bubble plots provided by --markers_file, default: FALSE",
	'onlypos'		,"y", 2, "logical","Optional parameter: only detect positively expressed marker genes, default: FALSE",
	'compare_list'		,"C", 2, "character","Optional parameter: Differential expression analysis group list: compare each cluster's internal case&control comparison, input format: control_vs_case,A_vs_B",
	'compare_list2'		,"x", 2, "character","Optional parameter: Differential expression analysis grouping list: case and control overall comparison, input format: control_vs_case,A_vs_B",
	'FindAllMarkers'        ,"f", 2, "logical","Optional parameter: whether to perform marker gene detection, default: TRUE",
	'print'			,"P", 2, "character","Optional parameter: whether to output the gene cell expression matrix: counts, data, or scale.data; not output by default",
	'markers_file'          ,"g", 2,"character","Commonly used marker genes: input file: tabel segmentation, two columns: the first column is the cell type; the second column is the gene name"


), byrow=TRUE, ncol=5)

args=getopt(command)

if(!is.null(args$help) | (is.null(args$project.rdata) & is.null(args$count))) {
	cat(paste(getopt(command, usage = T), "\n"))
	q()
}

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
yy<-unlist(strsplit(script.name,split="/"))
Bin=paste(yy[1:(length(yy)-1)],collapse="/")

source(paste(Bin,"/function.R",sep=""))

library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(reticulate)

use_python("/public/analysis/hefei/xten/software/Anaconda3/envs/my5/bin/python")
py_module_available(module = 'umap')

if ( is.null(args$tax) ) {args$tax="human"}
if ( is.null(args$outdir) ) {args$outdir=getwd()}
if ( is.null(args$nFeature_lowlimit) ) {args$nFeature_lowlimit=200}
if ( is.null(args$nFeature_highlimit) ) {args$nFeature_highlimit=10000}
if ( is.null(args$mito_cutoff) ) {args$mito_cutoff=10}
if ( is.null(args$redcell_cutoff) ) {args$redcell_cutoff=10}
if ( is.null(args$step) ) {args$step="1,2,3,4"}
if ( is.null(args$pcs_number) ) {args$pcs_number=30}
if ( is.null(args$resolution) ) {args$resolution="0.9"}
if ( is.null(args$batch) ) {args$batch="TRUE"}
if ( is.null(args$top_num) ) {args$top_num=10}
if ( is.null(args$onlyplot) ) {args$onlyplot=FALSE}
if ( is.null(args$onlypos) ) {args$onlypos=FALSE}
if ( is.null(args$FindAllMarkers) ) {args$FindAllMarkers="TRUE"}
if ( is.null(args$print) ) {args$print="FALSE"}
if ( args$onlyplot){args$step="NA"}
if(!is.null(args$group_list)){groups <- as.character(unlist(strsplit(args$group_list,split=",")))}
if(!is.null(args$compare_list)){compares <- as.character(unlist(strsplit(args$compare_list,split=",")))}
if(!is.null(args$compare_list2)){compares2 <- as.character(unlist(strsplit(args$compare_list2,split=",")))}

#cat(args$onlyplot)
if(args$onlyplot){
	if(is.null(args$project.rdata) | is.null(args$markers_file)){
      		cat("if only plots,please provide --project.rdata and --markers_file\n")
      		q()
	}
}


mito_cutoff <- args$mito_cutoff
redcell_cutoff <- args$redcell_cutoff
nFeature_lowlimit <- args$nFeature_lowlimit
nFeature_highlimit <- args$nFeature_highlimit
steps <- as.integer(unlist(strsplit(args$step,split=",")))
samples <- as.character(unlist(strsplit(args$sample_list,split=",")))
outdir <- args$outdir
resolutions <- unlist(strsplit(args$resolution,split=","))
resolutions <- as.numeric(resolutions)
pcs_number <- args$pcs_number
scfvalue <- 10000

if (file.exists(outdir)){
	setwd(outdir)
} else {
    	dir.create(outdir)
    	setwd(outdir)
}

##################
star<-paste(paste(rep("*",100),collapse=""),paste(rep("*",100),collapse=""),paste(rep("*",100),collapse=""),sep="\n")
write.table(paste(star,"Seurat's operating environment variables",sep="\n"),file="Seurat.log",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=FALSE)
begins<-paste("the Seurat step begins at the time of ",date(),sep=" ");
write.table(begins,file="Seurat.log",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)

packages.list<-sessionInfo()
R.version=packages.list$R.version$version.string
packages=NULL
for(i in 1:length(packages.list$otherPkgs)){
        package_name=packages.list$otherPkgs[[i]]$Package
        package_version=packages.list$otherPkgs[[i]]$Version
        packages<-c(packages,paste(package_name,package_version,sep=" "))
}
packages_list<-paste(packages,collapse=",")
packages_list2<-paste("R packages: ",packages_list,sep="")
write.table(R.version,file="Seurat.log",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(packages_list2,file="Seurat.log",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)

parmaters <- NULL
if(!is.null(args$sample_list)){parmaters <-c(parmaters,paste("--sample_list ",args$sample_list,sep=""))}
if(!is.null(args$project_name)){parmaters <-c(parmaters,paste("--project_name ",args$project_name,sep=""))}
if(!is.null(args$cluster_file)){parmaters <-c(parmaters,paste("--cluster_file ",args$cluster_file,sep=""))}
if(!is.null(args$cluster_number)){parmaters <-c(parmaters,paste("--cluster_number ",args$cluster_number,sep=""))}
if(!is.null(args$barcode_file)){parmaters <-c(parmaters,paste("--barcode_file ",args$barcode_file,sep=""))}
if(!is.null(args$gene_file)){parmaters <-c(parmaters,paste("--gene_file ",args$gene_file,sep=""))}
if(!is.null(args$group_list)){parmaters <-c(parmaters,paste("--group_list ",args$group_list,sep=""))}
if(!is.null(args$compare_list)){parmaters <-c(parmaters,paste("--compare_list ",args$compare_list,sep=""))}
if(!is.null(args$compare_list2)){parmaters <-c(parmaters,paste("--compare_list2 ",args$compare_list,sep=""))}
parmaters <-c(parmaters,paste("--mito_cutoff ",mito_cutoff,sep=""),paste("--redcell_cutoff ",redcell_cutoff,sep=""),paste("--nFeature_lowlimit ",nFeature_lowlimit,sep=""),paste("--nFeature_highlimit ",nFeature_highlimit,sep=""),paste("--resolution ",args$resolution,sep=""),paste("--pcs_number ",args$pcs_number,sep=""),paste("--tax ",args$tax,sep=""),paste("--batch ",args$batch,sep=""),"--normalization_method LogNormalize")
parmaters_list<-paste(parmaters,collapse=",")
parmaters_list2<-paste("Seurat parmaters: ",parmaters_list,sep="")
write.table(parmaters_list2,file="Seurat.log",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)


#####################################################################
extDiffGene <- function (project=NULL,diffs=NULL){
        lists<-sort(unique(as.numeric(as.vector(project@active.ident))))
        data_mat_frame <- project@assays$RNA@data
        data_mat_frame <- expm1(data_mat_frame)
        allcluster_mean <- NULL
        groupname <- NULL

        for(i in 1:length(lists)){
                        subdiff <- diffs[diffs$cluster == lists[i],]
                        subdata <- data_mat_frame[subdiff$gene,]
                        cluster1_mean <- Matrix::rowMeans(subdata[,project@active.ident == lists[i]])
                        cluster1_other_mean <- Matrix::rowMeans(subdata[,project@active.ident != lists[i]])
                        cluster_mean <- data.frame(Symbol=names(cluster1_mean),Cluster_1_Mean_UMI_Counts=cluster1_mean,Cluster_1_Other_Mean_UMI_Counts=cluster1_other_mean,Cluster_1_Log2_Fold_Change=log2(cluster1_mean+1)-log2(cluster1_other_mean+1))
                        subdiff_mean <- merge(subdiff, cluster_mean, by.x = "gene", by.y = "Symbol", all.x = TRUE,sort=FALSE)
                        colnames(subdiff_mean)[2] <- "Cluster_1_Pvalue"
                        colnames(subdiff_mean)[6] <- "Cluster_1_AdjustPvalue"
                        subdiff_mean <- subdiff_mean[,c("gene","pct.1","pct.2","Cluster_1_Mean_UMI_Counts","Cluster_1_Other_Mean_UMI_Counts","Cluster_1_Log2_Fold_Change","Cluster_1_Pvalue","Cluster_1_AdjustPvalue")]
                        colnames(subdiff_mean)[1] <- "Symbol"
                        colnames(subdiff_mean) <- c("Symbol","pct.1","pct.2",paste("Cluster_",lists[i],"_Mean_UMI_Counts",sep=""),paste("Cluster_",lists[i],"_Other_Mean_UMI_Counts",sep=""),paste("Cluster_",lists[i],"_Log2_Fold_Change",sep=""),paste("Cluster_",lists[i],"_Pvalue",sep=""),paste("Cluster_",lists[i],"_AdjustPvalue",sep=""))
                        otfile <- paste("Cluster",lists[i],"_diff.xls", sep="")
                        write.table(subdiff_mean, sep="\t", file=otfile, row.names=FALSE, quote=FALSE)
                        n <- dim(subdiff_mean)[2]
                        otfile <- paste("Cluster",lists[i], "_diff_significant.xls", sep="")
                        write.table(subdiff_mean[subdiff_mean[,(n-1)] <= 0.05 & abs(subdiff_mean[,(n-2)]) >= 0.585,], sep="\t", file=otfile, row.names=FALSE, quote=FALSE)  
                        tmpname <- paste("Cluster", lists[i] , sep="")
                        cluster_mean <- Matrix::rowMeans(data_mat_frame[,project@active.ident == lists[i]])
                        allcluster_mean <- cbind(allcluster_mean, cluster_mean)
                        groupname <- c(groupname, tmpname)

                }
                colnames(allcluster_mean) <- groupname
                write.table(data.frame(Symbol=rownames(allcluster_mean),allcluster_mean), file="Cluster_level_mean.xls", sep="\t", row.names=FALSE, quote=FALSE)
                return(allcluster_mean)


}



Compare_Of_Clusters <- function(compares,project){
        data_mat_frame <- project@assays$RNA@data
        data_mat_frame <- expm1(data_mat_frame)
        clusters<-sort(unique(project@active.ident))
        for(m in 1:length(compares)){
                pairs <- as.character(unlist(strsplit(compares[m],split="_vs_")))
                number<-NULL
                for(i in 1:length(clusters)){
                        length1<-length(which(project@active.ident == clusters[i] & project$Group %in% pairs[1]))
                        length2<-length(which(project@active.ident == clusters[i] & project$Group %in% pairs[2]))
                        if(length1<3 | length2<3){number<-c(number,0);next}
                        subdiff <- FindMarkers(project,subset.ident =clusters[i],ident.1=pairs[1],ident.2=pairs[2],min.pct = 0.01,logfc.threshold = 0.1,group.by="Group",
                                test.use = "bimod",assay = "RNA",min.cells.group =1)
                        subdata <- data_mat_frame[rownames(subdiff),]
                        group1_mean <- Matrix::rowMeans(subdata[,project@active.ident == clusters[i] & project$Group %in% pairs[1]])
                        group2_mean <- Matrix::rowMeans(subdata[,project@active.ident == clusters[i] & project$Group %in% pairs[2]])
                        name1 <- paste(pairs[1],"_Mean_UMI_Counts",sep="")
                        name2 <- paste(pairs[2],"_Mean_UMI_Counts",sep="")
                        cluster_mean <- data.frame(Symbol=names(group1_mean), name1=group1_mean,name2=group2_mean,Log2_Fold_Change=log2(group1_mean+1)-log2(group2_mean+1))
                        colnames(cluster_mean)<-c("Symbol",name1,name2,"Log2_Fold_Change")
                        subdiff_mean <- merge(data.frame(Symbol=rownames(subdiff),subdiff), cluster_mean, by.x = "Symbol", by.y = "Symbol", all.x = TRUE,sort=FALSE)
                        colnames(subdiff_mean)[2] <- "Pvalue"
                        colnames(subdiff_mean)[6] <- "AdjustPvalue"
                        subdiff_mean2 <- subdiff_mean[,c("Symbol","pct.1","pct.2",name1,name2,"Log2_Fold_Change","Pvalue","AdjustPvalue")]
                        otfile <- paste("Cluster",clusters[i],"_diff_","between_",pairs[1],"_",pairs[2],".xls", sep="")
                        write.table(subdiff_mean2, sep="\t", file=otfile, row.names=FALSE, quote=FALSE)
                        n <- dim(subdiff_mean2)[2]
                        otfile <- paste("Cluster",clusters[i],"_diff_significant_","between_",pairs[1],"_",pairs[2],".xls", sep="")
                        write.table(subdiff_mean2[subdiff_mean2[,(n-1)] <= 0.05 & abs(subdiff_mean2[,(n-2)]) >= 0.585,], sep="\t", file=otfile, row.names=FALSE, quote=FALSE)
                        number<-c(number,dim(subdiff_mean2[subdiff_mean2[,n] <= 0.05 & abs(subdiff_mean2[,(n-2)]) > 0.5,])[1])

                }
		}
}




if(3 %in% steps){
	cat("Step3:FindClusters and FindAllMarkers\n")
	if(!is.null(args$project.rdata) & !(1 %in% steps) & !(2 %in% steps)){
		load(args$project.rdata)
		if(!is.null(args$project_name)){project<-get(args$project_name)}
	}
	if(args$batch){

        	project <- ScaleData(object = project, verbose = TRUE,features=rownames(project@assays$RNA@data),assay="RNA")
	}

	for(j in 1:length(resolutions)){
		output<-paste("clusters_resolution",resolutions[j],sep="")
		dir.create(output)
		setwd(output)
		if(args$FindAllMarkers){
			if(args$batch){
                        	DefaultAssay(object = project) <- "integrated"
                        	project <- FindClusters(object = project, resolution = resolutions[j])
                	}else{
				if("integrated" %in% names(project@assays)){
					DefaultAssay(object = project) <- "integrated"
				}else{
                        		DefaultAssay(object = project) <- "RNA"
				}
                        	project <- FindClusters(object = project, resolution = resolutions[j])
                	}

                	DefaultAssay(object = project) <- "RNA"
			idents<-as.factor(as.numeric(as.vector(project@active.ident))+1)
			names(idents) <- rownames(project@meta.data)
			project@active.ident <- idents
			printClusterFiles(project=project,sample_list=args$sample_list)
			
			diffs <- tryCatch({
                                diffs <- FindAllMarkers(object = project, logfc.threshold = 0.1, only.pos = args$onlypos, test.use="bimod", min.pct = 0.01,assay="RNA")
                        }, error = function(e) {
                                set.seed(1)
                                random<-sample(1:dim(project)[2],100000)
                                project2<-subset(project,cells=colnames(project[,random]))
                                diffs <- FindAllMarkers(object = project2, logfc.threshold = 0.1, only.pos = args$onlypos, test.use="bimod", min.pct = 0.01,assay="RNA")
				return(diffs)
                        })

			#diffs <- FindAllMarkers(object = project, logfc.threshold = 0.1, only.pos = args$onlypos, test.use="bimod", min.pct = 0.01,assay="RNA")
                	write.table(diffs, sep="\t", file="project.FindAllMarkers.xls",quote=FALSE, row.names=FALSE)

			project[['Cluster']] <- project@active.ident
			save(list=c("project","diffs"),file='project.rdata')

			index<-grep("^RP[L|S]",diffs$gene, ignore.case = FALSE)
			if(length(index)>1){
				diffs2<-diffs[-index,]
			}else{
				diffs2<-diffs
			}
			
			allcluster_mean <- extDiffGene(project=project,diffs=diffs)
			
		}
		cat("perform diff gene analysis\n")	
		if(is.null(args$group_list)){setwd(outdir);next}
		project[["Group"]]<-rep(0,dim(project)[2])
                for(i in 1:length(groups)){
                       	yy<-paste("-",i,"$",sep="")
                       	project[["Group"]][grep(yy,colnames(project)),]<-groups[i]
                }
		if(!is.null(args$compare_list)){Compare_Of_Clusters(compares,project)}
		setwd(outdir)
	
	}
}





