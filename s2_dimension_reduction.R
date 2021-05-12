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


if(2 %in% steps){
	 if(!is.null(args$project.rdata) & !(1 %in% steps)){
		load(args$project.rdata)
		if(!is.null(args$project_name)){project<-get(args$project_name)}

	 }
	if(dim(project)[2]>pcs_number){dims_final=pcs_number}else{dims_final=dim(project)[2]-1}
	if(args$batch){
		if(min(table(project$Sample))>200){filter=200}else{filter=min(table(project$Sample))-1}
		if(min(table(project$Sample))>100){weight=100}else{weight=min(table(project$Sample))-1}

		cat("Step2:multiple sample batch effect: IntegrateData\n")
                project.list <- SplitObject(object = project, split.by = "Sample")
		cat("normalization_method is LogNormalize\n")
                hvf.info <-NULL
                for (i in 1:length(x = project.list)){
                       	project.list[[i]] <- NormalizeData(object = project.list[[i]], scale.factor = scfvalue)
                       	project.list[[i]] <- FindVariableFeatures(object = project.list[[i]],selection.method = "vst")
                       	hvf.info <-c(hvf.info,project.list[[i]]@assays$RNA@var.features)
                       	var.features=project.list[[i]]@assays$RNA@var.features
                       	output<-paste(names(project.list[i]),"_variable_features.xls",sep="")
                       	write.table(data.frame(Genes=var.features), sep="\t", file=output,quote=FALSE, row.names=FALSE)
               	}
                hvf.info2<-unique(hvf.info)
                write.table(data.frame(Genes=hvf.info2), sep="\t", file="Allsamples_variable_features.xls",quote=FALSE, row.names=FALSE)
		project <- tryCatch({
			project.anchors <- FindIntegrationAnchors(object.list = project.list, dims = 1:dims_final,anchor.features = hvf.info2,k.filter=filter)
			project <- IntegrateData(anchorset = project.anchors, dims = 1:dims_final,k.weight=weight)
			}, error = function(e) {
			cat("choose a sample as reference to intergrate data\n")
			project.list.cells <- NULL
			for(i in 1:length(x = project.list)){
				project.list.cells[i] <- dim(project.list[[i]])[2]
			}
			ref_id <- which.max(project.list.cells)                
			project.anchors <- FindIntegrationAnchors(object.list = project.list, dims = 1:dims_final,anchor.features = hvf.info2,k.filter=filter,
				reference = ref_id)
			project <- IntegrateData(anchorset = project.anchors, dims = 1:dims_final,k.weight=weight)
			return(project)
		})
                assay_name="integrated"

        }else{
		cat("Step2:no batch effect: NormalizeData and FindVariableFeatures\n")
                project <- NormalizeData(object = project, scale.factor = scfvalue)
                project <- FindVariableFeatures(object = project,selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.1, Inf),dispersion.cutoff=c(1, Inf))
                write.table(data.frame(Genes=project@assays$RNA@var.features), sep="\t", file="Allsamples_variable_features.xls",quote=FALSE, row.names=FALSE)
		assay_name="RNA"
		if("integrated" %in% names(project@assays)){
			assay_name="integrated";
		}
        }
	cat("running ScaleData and  RunPCA\n")
	project <- ScaleData(object = project, verbose = TRUE,features=rownames(project[[assay_name]]@data),assay=assay_name)
	project <- RunPCA(object = project, npcs = dims_final, verbose = TRUE,assay=assay_name,features=VariableFeatures(project,assay=assay_name))

	#######The number of PCA components selection method 1: draw a heat map of the positive and negative contribution of each principal component to the gene that contributes the most.
	
	pdf(file="DimHeatmap_PC.pdf",width=15,height=45)
	DimHeatmap(object = project, dims = 1:dims_final, cells = 500, balanced = TRUE)
	dev.off()

	#######PCA component number selection method 2: Randomly select a subset of the data (default is 1%) to perform PCA dimensionality reduction, construct an empty distribution of feature scores, and repeat this process.

	project <- JackStraw(object=project,reduction = "pca",dims = dims_final,assay = assay_name)
	project <- ScoreJackStraw(object = project, reduction = "pca", dims = 1:dims_final,score.thresh = 1e-05)
	pdf(file="JackStrawPlot_PC.pdf",width=15)
	p <- JackStrawPlot(object = project, dims = 1:dims_final)
	print(p)
	dev.off()

	#######The third method of selecting the number of PCA components: sort the principal components based on the percentage of variance explained by each (ElbowPlot function).

	pdf(file="ElbowPlot_PC.pdf",width=15)
	p <- ElbowPlot(project, ndims = dims_final, reduction = "pca")
	print(p)
	dev.off()
	
	#The FindNeighbors function is used to construct a KNN graph based on the Euclidean distance in the PCA space, and the edge weight between any two units is refined according to the shared overlap (Jaccard similarity) in its local neighborhood. The FindClusters function groups cell iterations together and contains a resolution parameter that is used to set the "granularity" of downstream clusters. Increased values will result in more clusters.

	cat("running FindNeighbors,RunTSNE and RunUMAP\n")
	project <- FindNeighbors(object = project, reduction = "pca", dims = 1:dims_final,assay=assay_name)

	#######TSNE nonlinear dimensionality reduction

	project <- RunTSNE(object = project, reduction = "pca",dims=1:dims_final,assay=assay_name,check_duplicates=FALSE)
	project <- RunUMAP(object = project, reduction = "pca",dims = 1:dims_final,assay=assay_name,umap.method = "umap-learn",metric = "correlation")
}


