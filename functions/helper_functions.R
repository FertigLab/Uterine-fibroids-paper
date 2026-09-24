#Rownames of de_res are gene names
annotate <- function(de_res, orthomap=FALSE) {
  #cellphoneDB secreted
  secreted <- secreted_genes(orthomap)
  de_res$cpdb.secreted <- rownames(de_res) %in% secreted
  
  #Surfaceome
  surface <- surfaceome(orthomap)
  if(orthomap) {
    de_res$Surfaceome <- rownames(de_res) %in% surface$Mouse.gene.name
    de_res$Surfaceome.Membranome.Almen.main.class <- surface[match(rownames(de_res), surface$Mouse.gene.name), "Membranome.Almen.main.class"]
  } else {
    de_res$Surfaceome <- rownames(de_res) %in% surface$UniProt.gene
    de_res$Surfaceome.Membranome.Almen.main.class <- surface[match(rownames(de_res), surface$UniProt.gene), "Membranome.Almen.main.class"]
  }
  
  #Matrisome
  if(orthomap) {
    matrix <- matrisome(org = "Mouse")
  } else {
    matrix <- matrisome(org = "Human")
  }
  de_res <- de_res %>%
    rownames_to_column() %>%
    left_join(matrix, join_by(rowname == Gene.symbol)) %>%
    column_to_rownames()
  
  return(de_res)
}


secreted_genes <- function(orthomap = FALSE) {
  cpdb_genes <- read.csv('ref/cpdb/v4.0.0/gene_input.csv', stringsAsFactors=F)
  cpdb_prot <- read.csv('ref/cpdb/v4.0.0/protein_input.csv', stringsAsFactors=F)
  sec_prot <- cpdb_prot$uniprot[cpdb_prot$secreted == 'True']
  sec_genes <- cpdb_genes$gene_name[cpdb_genes$uniprot %in% sec_prot]
  sec_genes <- unique(sec_genes)
  if(orthomap) {
    ortho_map <- readRDS("ref/HumanMouseMap/human_mouse_orthologs_mart_export.rds")
    sec_genes <- unique(ortho_map$Mouse.gene.name[ortho_map$Gene.name %in% sec_genes])
  }
  return(sec_genes)
}

surfaceome <- function(orthomap = FALSE) {
  surface <- read.csv("ref/Human_surfaceome/human_surfaceome.csv", header = TRUE, stringsAsFactors = F, skip = 1)
  surfaceome <- surface[,c("UniProt.gene", "Membranome.Almen.main.class")]
  if(orthomap) {
    ortho_map <- readRDS("ref/HumanMouseMap/human_mouse_orthologs_mart_export.rds")
    ortho_map.sub <- ortho_map %>%
      select("Gene.name", "Mouse.gene.name") %>%
      unique()
    ortho_map.sub <- ortho_map.sub[!duplicated(ortho_map.sub$Gene.name), ]
    surface$Mouse.gene.name <- ortho_map.sub[match(surface$UniProt.gene, ortho_map.sub$Gene.name), "Mouse.gene.name"] #Can use join
    surfaceome <- surface %>%
      select("Mouse.gene.name", "Membranome.Almen.main.class") %>%
      filter(!is.na(Mouse.gene.name)) %>%
      unique() 
    surfaceome <- surfaceome %>%
      filter(!(Mouse.gene.name %in% c("Or5ak24", "Frg1") & Membranome.Almen.main.class == "Unclassified"))
  }
  return(surfaceome)
}

matrisome <- function(org) {
  if(org == "Human") {
    m <- read.csv("ref/Matrisome/human.csv", header = TRUE, stringsAsFactors = FALSE)
  }else if(org == "Mouse") {
    m <- read.csv("ref/Matrisome/mouse.csv", header = TRUE, stringsAsFactors = FALSE)
  }
  return(m)
}


get_geneSet <- function(geneSet, org) {
  if(geneSet == "GO:BP") {
    geneSets = msigdbr(species = org, category = "C5", subcategory = "GO:BP")
  }else if(geneSet == "H") {
    geneSets = msigdbr(species = org, category = "H")
  }else if(geneSet == "KEGG") {
    geneSets = msigdbr(species = org, category = "C2", subcategory = "CP:KEGG")
  }else if(geneSet == "REACTOME") {
    geneSets = msigdbr(species = org, category = "C2", subcategory = "CP:REACTOME")
  }
  m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
  return(m_list)
}
