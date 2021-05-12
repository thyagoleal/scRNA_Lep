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

if(1 %in% steps){##################filter cells based on mito,ribo,red cells,double cells,nFeature lowlimit and nFeature highlimit
	if(!is.null(args$count)){

		project_mx <- fread(args$count,sep="\t", header=TRUE,data.table=F)
		colnames(project_mx) <- gsub("\\.", "-", colnames(project_mx))
		rownames(project_mx) <- project_mx[,1]
		project_mx[,1] <- NULL
		project <- CreateSeuratObject(counts = project_mx, min.cells = 1, min.features = 0, project = "project")

	}else if(!is.null(args$project.rdata)){

		load(args$project.rdata)
		if(!is.null(args$project_name)){project<-get(args$project_name)}

	}

	if(is.null(args$cluster_number) & is.null(args$barcode_file)){

		cat("step1:filter cells based on mito,ribo,red,double cells,nFeature lowlimit and nFeature highlimit\n")
		genes_list <- read.table(paste(Bin,"/rRNA.gene.list",sep=""),header=TRUE)
                known_tax <- unique(genes_list$Tax)
		if(args$tax %in% known_tax){
			index <- which(genes_list$Tax %in% args$tax & genes_list$Type %in% "mtRNA")
                        mito_genes <- unlist(strsplit(as.character(genes_list[index,3]),split=","))
                        index <- which(genes_list$Tax %in% args$tax & genes_list$Type %in% "rRNA")
                        ribo_genes <- unlist(strsplit(as.character(genes_list[index,3]),split=","))
                        index <- which(genes_list$Tax %in% args$tax & genes_list$Type %in% "redcell")
                        red_genes <- unlist(strsplit(as.character(genes_list[index,3]),split=","))
                        mito_genes <- mito_genes[which(mito_genes %in% rownames(project))]
			ribo_genes <- ribo_genes[which(ribo_genes %in% rownames(project))]
			red_genes <- red_genes[which(red_genes %in% rownames(project))]
			if(length(mito_genes)>=1){
				project[["percent.mt"]] <- PercentageFeatureSet(project, features=mito_genes)
			}else{
				project[["percent.mt"]] <- rep(0,dim(project)[2])
			}
			if(length(ribo_genes)>=1){
				project[["percent.ribo"]] <- PercentageFeatureSet(project, features=ribo_genes)
			}else{
				project[["percent.ribo"]] <- rep(0,dim(project)[2])
			}
			if(length(red_genes)>=1){
				project[["percent.redcell"]] <- PercentageFeatureSet(project, features=red_genes)
			}else{
				project[["percent.redcell"]] <-rep(0,dim(project)[2])
			}
		}
		project[["Sample"]]<-rep(0,dim(project)[2])
		if (!is.null(args$group_list)){project[["Group"]]<-rep(0,dim(project)[2])}
		for(i in 1:length(samples)){
        		yy<-paste("-",i,"$",sep="")
        		project[["Sample"]][grep(yy,colnames(project)),]<-samples[i]
 			if (!is.null(args$group_list))	{project[["Group"]][grep(yy,colnames(project)),]<-groups[i]}
		}
		barcode_info <- project@meta.data
		barcode_info[,1] <- NULL

		########Filter potential cell doublets
		if(!is.null(args$double_file)){
			cat("filter double cells\n")
                        double1 <- read.table(args$double_file,header=TRUE,sep="\t",row.names=1)
                        barcode_info <- cbind(barcode_info,double1[colnames(project),])
                        index <- which(double1$predicted_doublets %in% 0)
			valid_cells <- colnames(project[,which(colnames(project) %in% rownames(double1[index,]))])
			project <- subset(project,cells=valid_cells)
                }

		write.table(data.frame("Barcode"=rownames(barcode_info),barcode_info,check.names=FALSE),file="Barcode_Infor.xls",quote = FALSE,sep="\t", na = "NA", row.names = FALSE)

		QC_index<-c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.redcell")
		QC_index<-QC_index[which(QC_index %in% colnames(barcode_info))]
		pdf(file="QC_index_by_Violinplot.pdf",width=3*length(QC_index)+2)
		p<-VlnPlot(project, features = QC_index, ncol = length(QC_index),group.by = "Sample")
		print(p)
		dev.off()
		if("percent.ribo" %in% colnames(project@meta.data)){
			project <- subset(project, subset = nFeature_RNA > nFeature_lowlimit & nFeature_RNA < nFeature_highlimit & percent.mt < mito_cutoff & percent.redcell < redcell_cutoff)
		}else{
			project <- subset(project, subset = nFeature_RNA > nFeature_lowlimit & nFeature_RNA < nFeature_highlimit)
		}
		 project <- project[Matrix::rowSums(project)>0,]

	}else if(!is.null(args$project.rdata) & !is.null(args$cluster_number) & is.null(args$barcode_file)){


		cat("step1:filter cells based on parmaters --cluster_number: only recluster cells labeled cluster number ",args$cluster_number,"\n")
        	remain <- unlist(strsplit(args$cluster_number,split=","))
		project <- subset(project,idents=remain)
		project <- project[Matrix::rowSums(project)>0,]

	}else if(!is.null(args$barcode_file)  & is.null(args$cluster_number)){


		cat("step1:filter cells: only recluster cells in the file ",args$barcode_file,sep="")
        	barcode<-read.table(args$barcode_file,header=TRUE,stringsAsFactors=FALSE)
        	colnames(barcode)[1]<-"Barcode"
		project <- subset(project,cells=as.character(barcode$Barcode))
		project <- project[Matrix::rowSums(project)>0,]
	}
	
}
		