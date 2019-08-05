
# Color scale
gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

# Plot continuous values in a UMAP from meta.data column
# x is a seurat object
# returns a ggplot object
plot_umap_c <- function(x, 
						colname="percent.mt", 
						legend_title="MT%", 
						palette="Blues"){
	require(ggplot2)
	require(RColorBrewer)
	datf <- data.frame(feat=x@meta.data[,colname],
					   UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
					   UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
	p <- ggplot(datf, aes(x=UMAP1, y=UMAP2, color=feat)) + 
	geom_point() + theme_classic() +
	theme(axis.text=element_blank(), axis.ticks=element_blank()) +
	scale_colour_distiller(palette = palette, name=legend_title)
	return(p)
}

# Plot discrete values in a UMAP from meta.data column
# x is a seurat object
# returns a ggplot object
plot_umap_d <- function(x, colname="orig.ident", legend_title="orig.ident",
						colors=NULL, legend=TRUE, alpha=1){
	require(ggplot2)
	require(RColorBrewer)
	datf <- data.frame(feat=x@meta.data[,ct_col], 
					   UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
					   UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])

	clusters <- sort(unique(datf$feat))
	# Define colors
	if (is.null(colors)){
		nb.cols <- length(clusters)
		mycolors <- gg_color_hue(nb.cols)
	} else {
		mycolors <- colors
	}
	p <- ggplot(datf, aes(x=UMAP1, y=UMAP2)) + 
	geom_point(aes(color=factor(Cluster)), alpha=alpha) + theme_classic() + scale_color_manual(values = mycolors, name = legend_title) + 
	theme(axis.text=element_blank(), axis.ticks=element_blank(), text=element_text(size=20))
	if (!legend){
		p <- p + theme(legend.position = "none")
	}
	return(p)
}

# Plot gene expression values in a UMAP
# x is a seurat object
# returns a ggplot object
plot_umap_gene <- function(x, 
						   gene, 
						   legend_title="UMI", 
						   palette="Blues", 
						   assay=NULL){
	require(ggplot2)
	require(RColorBrewer)

	expr <- GetAssayData(x, slot="data", assay=assay)
	
	datf <- data.frame(feat=expr[gene,],
					 UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
					 UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
	p <- ggplot(datf, aes(x=UMAP1, y=UMAP2, color=feat)) + 
	
	geom_point() + theme_classic() + ggtitle(gene) + 
	theme(axis.text=element_blank(), axis.ticks=element_blank()) +
	scale_colour_distiller(palette = palette, name=legend_title)
	return(p)
}

# Plot cell types in UMAP with labels
# x is a Seurat object
# returns a ggplot object
plot_umap_labels <- function(x, 
							 colname=NULL, 
							 legend_title="Cell Type",
							 colors=NULL){
	library(ggplot2)
	library(RColorBrewer)

	celltypes <- ifelse(is.null(colname), x@active.ident, x@meta.data[,colname])
	df <- data.frame(Cluster=celltypes, 
					 UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
					 UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])

	# Get median UMAP points for each cluster
	clusters <- sort(unique(df$Cluster))
	meds <- lapply(clusters, function(ct) {
				   umap <- df[df$Cluster == ct, c("UMAP1", "UMAP2")]
				   meds <- apply(umap, 2, median)
				   return(meds)
					 })
	names(meds) <- clusters

	# Define colors
	if (is.null(colors)){
		nb.cols <- length(clusters)
		mycolors <- gg_color_hue(nb.cols)
	} else {
		mycolors <- colors
	}

	p <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
	geom_point(aes(color=factor(Cluster))) + 
	theme_classic() + 
	scale_color_manual(values = mycolors, name = legend_title) + 
	theme(axis.text=element_blank(), axis.ticks=element_blank())
	for (ct in clusters){
		p <- p + annotate(geom="label", x=meds[[ct]][1], y=meds[[ct]][2], label=as.character(ct))
	}
	return(p)
}

plot_pc_labels <- function(x, ct_col="SCT_snn_res.0.8", legend_title="Cell Type"){
	df <- data.frame(Cluster=x@meta.data[,ct_col], 
					 PC1=x@reductions$pca@cell.embeddings[,"PC_1"],
					 PC2=x@reductions$pca@cell.embeddings[,"PC_2"])
	# Get median PC points for each cluster
	clusters <- sort(unique(df$Cluster))
	meds <- lapply(clusters, function(ct) {
				   pcs <- df[df$Cluster == ct, c("PC1", "PC2")]
				   meds <- apply(pcs, 2, median)
				   return(meds)
					 })
	names(meds) <- clusters
	p <- ggplot(df, aes(x=PC1, y=PC2, color=Cluster)) + 
	geom_point() + theme_classic() + scale_color_discrete(name = legend_title) + 
	theme(axis.text=element_blank(), axis.ticks=element_blank())
	for (ct in clusters){
		p <- p + annotate(geom="label", x=meds[[ct]][1], y=meds[[ct]][2], label=as.character(ct))
	}
	return(p)
}

plot_pc_cont <- function(x, md_col="percent.mt", legend_title="MT%"){
	df <- data.frame(feat=x@meta.data[,md_col],
					 PC1=x@reductions$pca@cell.embeddings[,"PC_1"],
					 PC2=x@reductions$pca@cell.embeddings[,"PC_2"])
	p <- ggplot(df, aes(x=PC1, y=PC2, color=feat)) + 
	geom_point() + theme_classic() +
	theme(axis.text=element_blank(), axis.ticks=element_blank()) +
	scale_colour_distiller(palette = "Spectral", name=legend_title)
	return(p)
}

# Boxplot of MT%
boxplot_mt <- function(x, mt_col="percent.mt", ct_col="SCT_snn_res.0.8"){
	require(reshape2)
	df <- data.frame(Mito=x@meta.data[,mt_col], Cluster=x@meta.data[,ct_col])
	df <- reshape2::melt(df)
	p <- ggplot(df, aes(x=Cluster, y=value, fill=Cluster)) +
	geom_boxplot() + theme_classic() + 
	ylab("MT%")
	return(p)
}

# Gene expression heatmap
expr_heatmap <- function(x, de, ct_col, n_genes=10, n_cells=50, colors=NULL, breaks=NA){
	cell_types <- sort(unique(de$cluster))
	topn <- lapply(cell_types, function(i) de[de$cluster == i,"gene"][1:n_genes])
	topn <- unlist(topn)
	barcodes <- lapply(cell_types, function(i) {
					   sample(colnames(x)[x@meta.data[,ct_col] == i], size=n_cells)})
	barcodes <- unlist(barcodes)
	expr <- as.matrix(x@assays$SCT@data[topn, barcodes])
	expr <- t(scale(t(expr)))

	anno_col <- lapply(cell_types, function(i) (rep(i, n_cells)))
	anno_col <- as.data.frame(unlist(anno_col))
	rownames(anno_col) <- colnames(expr)
	colnames(anno_col) <- "Cell Type"

	if (is.null(colors)){
	p <- pheatmap(expr, annotation_col=anno_col, breaks=breaks, 
				  cluster_cols=FALSE, cluster_rows=FALSE,
				  show_rownames=FALSE, show_colnames=FALSE,
				  annotation_names_col = FALSE)
	} else{
		p <- pheatmap(expr, annotation_col=anno_col, breaks=breaks, 
					  cluster_cols=FALSE, cluster_rows=FALSE,
					  show_rownames=FALSE, show_colnames=FALSE,
					  annotation_names_col = FALSE,
					  annotation_colors=colors)
	}
	return(p)
}

expr_heatmap_genes <- function(x, genes, ct_col, n_cells=50, colors=NULL){
	cell_types <- sort(unique(x@meta.data[,ct_col]))
	barcodes <- lapply(cell_types, function(i) {
					   sample(colnames(x)[x@meta.data[,ct_col] == i], size=n_cells)})
	barcodes <- unlist(barcodes)
	expr <- as.matrix(x@assays$SCT@data[genes, barcodes])
	expr <- t(scale(t(expr)))

	anno_col <- lapply(cell_types, function(i) (rep(i, n_cells)))
	anno_col <- as.data.frame(unlist(anno_col))
	rownames(anno_col) <- colnames(expr)
	colnames(anno_col) <- "Cell Type"

	if (is.null(colors)){
	p <- pheatmap(expr, annotation_col=anno_col,
				  cluster_cols=FALSE, cluster_rows=FALSE,
				  show_rownames=FALSE, show_colnames=FALSE,
				  annotation_names_col = FALSE)
	} else{
		p <- pheatmap(expr, annotation_col=anno_col,
					  cluster_cols=FALSE, cluster_rows=FALSE,
					  show_rownames=FALSE, show_colnames=FALSE,
					  annotation_names_col = FALSE,
					  annotation_colors=colors)
	}
	return(p)
}

