
#libraries needed for DGE analysis 
library(R6)
library(limma)
library(edgeR)
library(clusterProfiler)
library(biomaRt)
library(dplyr)

DGEA <- R6Class(
  classname = "DGEA",
  public = list(
    initialize = function(gx_filename = NA, meta_filename = NA, save = FALSE) {
      private$gx_filename <- gx_filename
      private$meta_filename <- meta_filename
      private$save <- save
    },
    
    DGEA_analysis = function(save = private$save) {
      list_of_paths <- c(private$gx_filename, private$meta_filename)
      complete_info <- private$read_files(list_of_paths)
      metadata <- complete_info$sampleInfo
      gene_exp_data <- complete_info$gxData
      
      # Create design matrix considering confounding variables gender and age
      design <- model.matrix(~ 0 + etiology + gender + age, data = metadata)
      
      # Apply limma pipeline 
      fit <- lmFit(gene_exp_data, design)
      
      # Create contrasts matrix
      cont.matrix <- makeContrasts(DCMvsControl = etiologyDCM - etiologyNF, levels = design)
      
      # Fit contrasts
      fit2 <- contrasts.fit(fit, cont.matrix)
      
      # Apply empirical Bayes statistics
      ebFit <- eBayes(fit2, trend = TRUE)
      
      # Get results
      dgeRes <- topTable(ebFit, coef = 'DCMvsControl', number = nrow(gene_exp_data))
      
      # Retrieve gene symbols and gene names based on the provided Ensembl gene identifiers
      ensembl <- useEnsembl("ensembl", mirror = "useast")
      hsapiens <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
      
      gene_ids <- rownames(gene_exp_data)
      
      gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"), 
                         filters = "ensembl_gene_id", 
                         values = gene_ids, 
                         mart = hsapiens)
      
      #Eliminating unwated notation from the gene_info dataframe 
      gene_info$description <- sub("\\[Source:.*\\]", "", gene_info$description)
      
      #Adding ensemble_id as the rownames for gene_info
      rownames(gene_info) <- gene_info$ensembl_gene_id
      gene_info <- subset(gene_info, select = -1)
      
      #Merge gene_info and dgeRes by ensemble_id(rowname)
      dgeRes <- merge(gene_info, dgeRes, by = "row.names", all.y = TRUE)
      rownames(dgeRes) <- dgeRes$Row.names
      dgeRes <- subset(dgeRes, select = -1)
      
      #Select genes that were significant
      de_genes <- dgeRes[abs(dgeRes$logFC) > 1 & dgeRes$adj.P.Val < 0.05, ]
      de_genes <- de_genes %>% arrange(desc(abs(logFC)))
      
      # Save the de_genes as csv (if user chooses to)
      if (save == TRUE) {
        write.csv(de_genes, file = "de_genes.csv", row.names = TRUE)
      }
      
      return(list(gene_info = gene_info , de_genes = de_genes))
    }
  ),
  
  private = list(
    gx_filename = NA,
    meta_filename = NA,
    save = NA,
    
    # Read files based on file path and type. 
    read_files = function(file_path_list) {
      meta_keyword <- "sampledata"
      gx_keyword <- "geneexpression"
      sampleInfo = NA
      gxData = NA
      
      for (path in file_path_list) {
        type_of_file <- tools::file_ext(path)
        
        if (type_of_file == "txt" && grepl(gx_keyword, tolower(path))) {
          gxData <- read.delim(path, as.is = TRUE, row.names = 1) 
        } else if (type_of_file == "csv" && grepl(meta_keyword, tolower(path))) {
          sampleInfo <- read.csv(path, as.is = TRUE, row.names = 1)
        }
      }
      
      return(list(sampleInfo = sampleInfo, gxData = gxData))
    }
  )
)

#Libraries needed for Enrichment analysis
library(rWikiPathways)
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)


Enrich <- R6Class(
  classname = "Enrich",
  
  public <- list(
    ensemble_genes = NULL,
    de_genes = NULL,
    save_enrichment_tables = NULL,
    initialize = function(ensemble_genes = NA , de_genes = NA , save_enrichment_tables = FALSE){
      self$ensemble_genes = ensemble_genes
      self$de_genes = de_genes 
      self$save_enrichment_tables = save_enrichment_tables
    }, 
    
    #Perform enrichment analysis function
    enrich_analysis = function(){
      
      #Download all human genes from the Wikipathways archive
      gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", date = "20231110", format = "gmt")
      wp2gene <- readPathwayGMT(gmt)
      
      #Link the pathway_id with the entrez_id of the relevant gene
      wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
      
      #Link the pathway_id with the name of the pathway 
      wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
      
      #Add ENSEMBLE_ID as a column in both de_genes and ensemble_genes
      all_gene_esmbl_ids <- rownames(self$ensemble_genes)
      de_gene_esmbl_ids <- rownames(self$de_genes)
      rownames(self$ensemble_genes) <-NULL
      rownames(self$de_genes) <-NULL
      self$ensemble_genes$ENSEMBLE_ID <- all_gene_esmbl_ids
      self$de_genes$ENSEMBLE_ID <- de_gene_esmbl_ids
      
      #Change the ENSEMBLE_ID to ENTREZID
      all_genes <- unique(self$ensemble_genes[c("ENSEMBLE_ID","hgnc_symbol")])
      de.genes.entrez <- bitr(self$de_genes$ENSEMBLE_ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
      all.genes.entrez <- bitr(all_genes$ENSEMBLE_ID,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
      
      #Perform enrichment analysis using the WIKIPATHWAYS database
      enrichment_results_wiki <- clusterProfiler::enricher(
        de.genes.entrez$ENTREZID , 
        universe = all.genes.entrez$ENTREZID, 
        pAdjustMethod = "fdr", 
        pvalueCutoff = 0.05,        
        qvalueCutoff = 0.02,
        TERM2GENE = wpid2gene,
        TERM2NAME = wpid2name
      )
      
      #Create a dataframe out of the Enrich object
      enrich_dataframe_wiki = as.data.frame(enrichment_results_wiki)
      
      #Turn de.genes.entrez to vector 
      vector.de.genes.entrez = as.vector(de.genes.entrez)
      
      #Perform enrichment analysis using the REACTOME database
      enrich_results_react <- enrichGO(
        gene = vector.de.genes.entrez$ENTREZID, OrgDb = org.Hs.eg.db,
        ont = "BP",
        pvalueCutoff = 0.05
      )
      
      #Create a dataframe out of the Enrich object
      enrich_dataframe_react = as.data.frame(enrich_results_react)
    
      
      #Perform enrichment analysis using the KEGG database
      enrich_results_kegg <- enrichKEGG(gene = de.genes.entrez$ENTREZID,
                                        organism = "hsa", 
                                        keyType = "kegg")
      
      
      #Create a dataframe out of the Enrich object
      enrich_dataframe_kegg <- as.data.frame(enrich_results_kegg)
      
      #-----------------------------------------------------------------------------#
        #Pathway Enrichment Analysis: Hypergeometric distribution
        #-----------------------------------------------------------------------------#
        # The hypergeometric distribution gives the probability of observing at least x differentially expressed genes in a particular pathway by chance.
        # It measures how likely it is that the given amount of genes per pathway in the PEA might be a random observation.
        
        # n: The total number of genes in the dataset, 20781.
        
        # k: This is the number of success states in the drawn sample. So for us, it is the number of genes identified as interesting from the dataset.
        
        # q: This is the number of successes you are interested in. This is x, the number of the interesting genes annotated in the pathway.
        
      # m: This is the total number of success states in the population. Here, it is the total number of genes in a pathway.
      
      # Retrieve m (total number of genes):
      
      # Split the BgRatio string by "/" and take the first element, which is the number of total genes per pathway
      enrich_dataframe_kegg$Total_Genes <- sapply(strsplit(as.character(enrich_dataframe_kegg$BgRatio), "/"), "[", 1)
      enrich_dataframe_react$Total_Genes <- sapply(strsplit(as.character(enrich_dataframe_react$BgRatio), "/"), "[", 1)
      enrich_dataframe_wiki$Total_Genes <- sapply(strsplit(as.character(enrich_dataframe_wiki$BgRatio), "/"), "[", 1)
      
      # Convert the Total_Genes column to numeric
      enrich_dataframe_kegg$Total_Genes <- as.numeric(enrich_dataframe_kegg$Total_Genes)
      enrich_dataframe_react$Total_Genes <- as.numeric(enrich_dataframe_react$Total_Genes)
      enrich_dataframe_wiki$Total_Genes <- as.numeric(enrich_dataframe_wiki$Total_Genes)
      
      n = nrow(self$ensemble_genes) # should be 20781
      k = nrow(self$de_genes) # should be 1085
      
      # Variables x and m are databse-specific
      ##KEGG
      x_kegg = enrich_dataframe_kegg$Count
      m_kegg = enrich_dataframe_kegg$Total_Genes
      
      ##Reactome
      x_react = enrich_dataframe_react$Count
      m_react = enrich_dataframe_react$Total_Genes
      
      ##Wikipathways
      x_wiki = enrich_dataframe_wiki$Count
      m_wiki = enrich_dataframe_wiki$Total_Genes
      
      # We want to calculate the probability of getting at least x successes.
      # This is why we need to subtract the result of phyper from 1 and use x-1 instead of x. 
      # P(X >= x) = 1 - P(X <= x-1)
      
      # So the phyper function called on all databeses would look like this:
      enrich_dataframe_kegg$Prob <- 1 - phyper(q = x_kegg-1 , m = m_kegg, n = n, k = k)
      enrich_dataframe_react$Prob <- 1 - phyper(q = x_react-1 , m = m_react, n = n, k = k)
      enrich_dataframe_wiki$Prob <- 1 - phyper(q = x_wiki-1 , m = m_wiki, n = n, k = k)
      
      
      #Save the all enrichment dataframes per database if user chooses to. 
      if (self$save_enrichment_tables == TRUE){
        write.table(enrich_dataframe_wiki, file="wikipathway_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)
        write.table(enrich_dataframe_react, file="reactome_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)
        write.table(enrich_dataframe_kegg, file="kegg_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)
      }
      #Add gene names to enrichment objects for visualization (Important step!)
      edox_wikipath <- setReadable(enrichment_results_wiki, 'org.Hs.eg.db', 'ENTREZID')
      edox_kegg <- setReadable(enrich_results_kegg, 'org.Hs.eg.db', 'ENTREZID')
      
      
      #Reactome databases already provides the gene names 
      edox_react <- enrichPathway(gene=vector.de.genes.entrez$ENTREZID, pvalueCutoff=0.05, readable=T)
    
      #Return the enrichment objects for visualization. 
      return(list(edox_wikipath = edox_wikipath, edox_kegg = edox_kegg, edox_react = edox_react))
    }
    
    
  )
)



library(ggplot2)
plot_function = function( save_barplots = FALSE,
                         visualize_wikipathways_cnetplot = FALSE,
                         visualize_reactome_cnetplot = FALSE,
                         visualize_kegg_cnetplot = FALSE, 
                         show_wikipathways_barplot = FALSE,
                         show_reactome_barplot = FALSE,
                         show_kegg_barplot = FALSE, 
                         plots_showing = TRUE

                                                 ) {
if (sum(visualize_wikipathways_cnetplot,visualize_reactome_cnetplot,
        visualize_kegg_cnetplot,
        show_wikipathways_barplot,
        show_reactome_barplot,
        show_kegg_barplot) >=2) {
warning ("YOU ARE ONLY SEEING YOUR 1st GRAPH! IF YOU WANT TO SEE THE OTHERS,
         PICK ONE 'TRUE' STATEMENT AT A TIME")
}
  
#Create an instance of the DGEA class
  dgea = DGEA$new(gx_filename = "MAGNET_GeneExpressionData_CPM_19112020.txt", meta_filename = "MAGNET_SampleData_18112022.csv")
  
  #Get the info on all genes and the differentially expressed genes 
  complete_info = dgea$DGEA_analysis()
  
  #Seperate all_genes from differentially expressed genes in variables 
  all_genes = complete_info$gene_info
  de_genes = complete_info$de_genes
  
  #Create and instance of the Enrichment class
  enrich_instance = Enrich$new(all_genes, de_genes)
  
  #Get the enrichment objects of each database (e.g. Wikipathways,KEGG,Reactome)
  enrichment_objects = enrich_instance$enrich_analysis()
  
  #Save all the enrichment objects in a variable
  wikipathways_enr_object = enrichment_objects$edox_wikipath
  reactome_enr_object = enrichment_objects$edox_react
  kegg_enr_object = enrichment_objects$edox_kegg
  
  #Turn Enrich objects into dataframes 
  wikipathways_enr_dataframe = data.frame(wikipathways_enr_object)
  reactome_enr_dataframe = data.frame(reactome_enr_object)
  kegg_enr_dataframe = data.frame(kegg_enr_object)
   
  #If plots are not shown in the "Plots" window, turn off all graphic devices.
  if(plots_showing == FALSE){
    open_devices = dev.cur()
    device_index = as.integer(open_devices)
    if (device_index > 1) {
    dev.off(as.integer(open_devices))
    }
    else{
    warning("The problem is NOT the amount of graphic devices open")
    }
  }
 wikipathway_barplot <- ggplot(wikipathways_enr_dataframe, aes(x = reorder(Description, -p.adjust), y = Count, fill = (p.adjust))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "blue", high = "red", trans = "reverse") +
    labs(x = "Pathway", y = "Count", fill = "p.adjust") +
    theme(axis.text = element_text(size = 11)) +
    ggtitle("Pathway Enrichment Analysis for All Differentially Expressed Genes")
 if (save_barplots == TRUE){
    ggsave(filename = "wikipathway_barplot.pdf", width = 10, height = 10,units = "in")
 }
    
 #Wikipathways based barplot 
  if (show_wikipathways_barplot == TRUE){
  print(wikipathway_barplot)
  }
  
  reactome_barplot <- ggplot(reactome_enr_dataframe, aes(x = reorder(`Description`, -`p.adjust`),y = Count, fill = (p.adjust))) +
     geom_bar(stat = "identity") +
     coord_flip() +
     scale_fill_gradient(low = "blue", high = "red" , trans = "reverse") +
     labs(x = "Pathway", y = "Count", fill = "p.adjust") +
     theme(axis.text = element_text(size = 11)) +
     ggtitle("Reactome-based Enrichment Analysis for All Differentially Expressed Genes")
  if (save_barplots == TRUE){
     ggsave(filename = "reactome_barplot.pdf", width = 10, height = 10, units = "in")
  }
  #Reactome based barplot   
  if (show_reactome_barplot == TRUE){
  print(reactome_barplot)
  }
 
   
   kegg_barplot <- ggplot(kegg_enr_dataframe, aes(x = reorder(`Description`, -`p.adjust`),y = Count, fill = (p.adjust))) +
     geom_bar(stat = "identity") +
     coord_flip() +
     scale_fill_gradient(low = "blue", high = "red" , trans = "reverse") +
     labs(x = "Pathway", y = "Count", fill = "p.adjust") +
     theme(axis.text = element_text(size = 11)) +
     ggtitle("Kegg-based Enrichment Analysis for All Differentially Expressed Genes") 
   if (save_barplots == TRUE){
     ggsave(filename = "kegg_barplot.pdf", width = 10, height = 10 , units = "in")
   }
  #Kegg based barplot 
  if (show_kegg_barplot == TRUE){
   print(kegg_barplot)
   }
     
     

        
  #Visualize the most significant pathways for each database individually 
  if (visualize_wikipathways_cnetplot ==TRUE) {
  #Visualize wikipathways cnetplot (network graph)
  cnetplot(wikipathways_enr_object, categorySize="pvalue", color.params = list(foldChange = de_genes$logFC), cex_label_category = 0.7, cex_label_gene = 0.7, max.overlaps = 5)
  }
  
  else if (visualize_reactome_cnetplot ==TRUE) {
  #Visualize reactome cnetplot (network graph)
  cnetplot(reactome_enr_object, categorySize="pvalue", color.params = list(foldChange = de_genes$logFC), cex_label_category = 0.7, cex_label_gene = 0.7, max.overlaps = 5)
   }
  
  else if (visualize_kegg_cnetplot ==TRUE) {
  #Visualize kegg cnetplot (network graph)
  cnetplot(kegg_enr_object, categorySize="pvalue", color.params = list(foldChange = de_genes$logFC), cex_label_category = 0.7, cex_label_gene = 0.7, max.overlaps = 5)
  }
  
}


#Plot function 
##Notes: If user wants to save the cnetplots specifcally, they will have to do so manually 
##Notes: User can choose to save Barplots by setting save_barplots = TRUE
plot_function(
             save_barplots = FALSE,
             visualize_wikipathways_cnetplot = TRUE,
             visualize_reactome_cnetplot = TRUE,
             visualize_kegg_cnetplot = FALSE, 
             show_wikipathways_barplot = FALSE,
             show_reactome_barplot = FALSE,
             show_kegg_barplot = FALSE, 
             plots_showing = FALSE)

  