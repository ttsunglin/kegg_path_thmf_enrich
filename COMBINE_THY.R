 ########################
 # START OF CALCULATION #
 # last update 20190619 #
 ########################
#------------------^^LIBRARY^^--------------------------


suppressMessages(library(KEGGREST))
suppressMessages(library(dplyr))
suppressMessages(library(rvest))
suppressMessages(library(readr))


#------------------^^Create_output_folder^^----
dir.create(file.path("output"), showWarnings = FALSE)
#------------------^^CONTROLLER^^-------------------------


a <- 1
#b <- defined later in *constant_mmu_path or *chosen_pathway_input (optional)


#------------------^^SHARED_CONSTANT^^----------------------


#----------------------*constant_mmu_path--------------------

#building mmu pathway dataframe
library(KEGGREST)
path_mmu <- keggList("pathway", "mmu")           #downloading all pathway of mmu 
path_mmu_id <- rownames(data.frame(path_mmu))    #taking only the map ID
path_mmu_amount <- nrow(data.frame(path_mmu_id)) #amount of all mmu pathways
b <- path_mmu_amount
b
#there are 326 mmu pathways 

#----------------------*chosen_pathway_input (optional)---------------

#input chosen pathway kegg ID
#NEED TO BE UPDATE FROM THE KEGG
#https://www.genome.jp/kegg/pathway.html#metabolism

library(readr)
#metabolism_pathway <- suppressMessages(read_csv("metabolism_pathway.csv",
#                                                col_names = FALSE))
#metabolism_pathway2 <- data.frame(data.frame(metabolism_pathway)[-c(1:9),])
#row.names(metabolism_pathway2) <- metabolism_pathway2[,1]
#path_mmu_id3 <- rownames(data.frame(metabolism_pathway2))
path_mmu_id3 <- c(1)
if(length(path_mmu_id3)>1){
  #path_mmu_id4 <- data.frame(path_mmu_id)
  path_mmu_id <- intersect(path_mmu_id, path_mmu_id3)
  path_mmu_amount <- nrow(data.frame(path_mmu_id))
  b <- path_mmu_amount
  message("Run with chosen pathways.")
  b
}else{
  message("Run with all pathways.")
  b
  }
# 85 chosen pathways
#----------------------*constant_num_in_path---------------------

#-------------------------/compound/----------------------------
#calculating compound amount in each pathway
PATHCOMPINFO <- function(x){
  e <- keggGet(x)
  f <- list(e[[1]]$COMPOUND)
  nrow(data.frame(f))
}
path_comp <- apply(data.frame(path_mmu_id[a:b]), 1,PATHCOMPINFO)
comp_in_path_df <- data.frame(y1 = path_mmu_id[a:b], y2 = path_comp)
colnames(comp_in_path_df) <- c("PATHWAY_ID", "PATHWAY_COMP_NUM")
#list of compound amount in each pathway -> comp_in_path_df
#-------------------------/gene/-------------------------------
#calculating gene amount in each pathway
PATHGENEINFO <- function(x){
  e <- keggGet(x)
  f <- list(e[[1]]$GENE)
  nrow(data.frame(f))/2
}
path_gene <- apply(data.frame(path_mmu_id[a:b]), 1,PATHGENEINFO)
gene_in_path_df <- data.frame(y1 = path_mmu_id[a:b], y2 = path_gene)
colnames(gene_in_path_df) <- c("PATHWAY_ID", "PATHWAY_GENE_NUM")
#list of gene amount in each pathway -> gene_in_path_df
#-------------------------/combine/----------------------------
#calculating gene and compound amount in each pathway
combine_in_path_df <- data.frame(y1 = comp_in_path_df[,1],
                                 y2 = comp_in_path_df[,2] + gene_in_path_df[,2])
colnames(combine_in_path_df) <- c("PATHWAY_ID", "PATHWAY_COMBINE_NUM")
#list of COMBINATION amount in each pathway -> COMBINE_in_path_df

#----------------------*constant_without_repeat--------------------

#-------------------------/compound/----------------------------
#calculating all pathway compounds without repeats
PATHCOMPINFO2 <- function(x){
  e <- keggGet(x)
  f <- list(e[[1]]$COMPOUND)
  f
}
path_comp2 <- apply(data.frame(path_mmu_id[a:b]), 1,PATHCOMPINFO2)
path_comp_df2 <- data.frame(unlist(path_comp2))
comp_in_path_df_total_union <- nrow(unique(path_comp_df2))
comp_in_path_df_total_union
#there are 3385 different compounds in the total 326 pathways 

PATHCOMPINFO3 <- function(x){
  rownames(data.frame(unlist(PATHCOMPINFO2(x))))
}
path_comp3 <- apply(data.frame(path_mmu_id[a:b]), 1,PATHCOMPINFO3)
path_comp_df3 <- data.frame(unlist(path_comp3))
comp_in_path_df_total_union <- nrow(unique(path_comp_df3))
comp_in_path_df_total_union
#for matching related compounds, path_comp_df3 contains all compounds in path (with repeats)

#-------------------------/gene/-------------------------------
#calculating all pathway genes without repeats
PATHGENEINFO2 <- function(x){
  e <- keggGet(x)
  f <- data.frame(e[[1]]$GENE)
  ff <- f[!row(data.frame(f))%%2 == 0,] 
  ff
}
path_gene2 <- apply(data.frame(path_mmu_id[a:b]), 1,PATHGENEINFO2)
path_gene_df2 <- data.frame(unlist(path_gene2))
gene_in_path_df_total_union <- nrow(unique(path_gene_df2))
gene_in_path_df_total_union
#there are 8585 different genes in the total 326 pathways

#-------------------------/combine/----------------------------
#calculating all pathway genes without repeats
combine_in_path_df_total_union <- gene_in_path_df_total_union + comp_in_path_df_total_union
combine_in_path_df_total_union
#there are 11970 different combine in the total 326 pathways

#----------------------*output_compound_gene_list-----
write.csv(unique(path_comp_df3), file = paste(getwd(), "/output/compound_list.csv", sep = ""))
write.csv(unique(path_gene_df2), file = paste(getwd(), "/output/gene_list.csv", sep = ""))
compound_list <- suppressMessages(read_csv(paste(getwd(), "/output/compound_list.csv", sep = "")))
gene_list <- suppressMessages(read_csv(paste(getwd(), "/output/gene_list.csv", sep = "")))
compound_list <- compound_list[,-1]
gene_list <- gene_list[,-1]

#-----------------^^START_CALCULATION^^--------------------

#----------------------*input-----------------------------
library(readr)
THY_CSID <- read_csv("20190524 comp_MF_THY_FC1.5_P0.05.csv", 
                               col_names = FALSE)
THY_CSID <- unique(THY_CSID)
gene_duplicate_check <- read_csv("gene_duplicate_check.csv")
gene_duplicate_check2 <- gene_duplicate_check[,-1]
gene_duplicate_check3 <- gene_duplicate_check2[!(gene_duplicate_check2$dup == TRUE &
                                                 gene_duplicate_check2$`Input Type` != "current symbol"),]
gene_duplicate_check4 <- gene_duplicate_check3[gene_duplicate_check3$p_value <= 0.05,]

#input data table
library(readr)
CSID_ID_THY <- data.frame(THY_CSID$X2 )#COMPOUND IN CSID (ONLY NUMBER) p<0.05 fc>1.5 
THY_KEGGID <- data.frame(gene_duplicate_check4$`Entrez Gene ID`) #GENE IN KEGG ID (ONLY NUMBER)
#Gene KEGG ID from GENE symbol matched on MGI under "all symbols " 
#mutiple hit dealed by checking, kickout non-matched gene. Gene expression set criteria thy>pc fc>3 p<0.05
#----------------------*convert_ID (CSID only)-------------------

#this part processing with too many data might cause computer shutdown
#convert compound CSID to KEGG ID 
#conneting chemspider url
library(rvest)
c <- 1                                                             #from a to b row, CSID with only numbers
d <- 100#nrow(CSID_ID_THY)

URL <- function(x){                                               #combining the url, note that the url is backstage of Chemspider
  url <- paste("http://www.chemspider.com/ibcontent.ashx?csid=",
               x,
               "&type=ds&ds_types=",
               sep = ""
  )
  doc <- read_html(url)                                           #node of ID in interest, found with SelectorGaget in Chrome
  node <- '//*[(@id = "DS5")]//a'
  html_node <- html_nodes(doc, xpath = node)                      #note the url 
  html_catch <- html_text(html_node)
  html_df <- data.frame(html_catch)
  kegg_row <- which(html_df$html_catch == "KEGG" )                #looking for KEGG ID under "KEGG"
  c1 <- grepl("C",substring(as.character(html_df[ kegg_row+1, ]),1,1))
  c2 <- grepl("C",substring(as.character(html_df[ kegg_row+2, ]),1,1))
  d1 <- grepl("D",substring(as.character(html_df[ kegg_row+1, ]),1,1))
  d2 <- grepl("D",substring(as.character(html_df[ kegg_row+2, ]),1,1))
  if(
    c1&&c2 == TRUE                                                #multiple matches
  ){"Match error**"
  }else if(                                                     
    c1 == TRUE                                                    #check the first row under "KEGG"
  ){as.character(html_df[ kegg_row+1, ])
  }else if(
    c2 == TRUE                                                    #check for second row
  ){as.character(html_df[ kegg_row+2, ])
  }else if(
    d1||d2 == TRUE                                                #if no matches, check if ony KEGG ID start with "D", which is short for "Drug"
  ){"Drug only"                                                   #This step depletes KEGG ID with "D" start
  }else{
    NA                                                            #if CSID does not match with KEGG ID
  }
}

#may crash if calculating all in the same apply function, so separated
kegg_id_thy_1 <- apply(data.frame(CSID_ID_THY[c:d,1]), 1,URL)       
kegg_id_thy_2 <- apply(data.frame(CSID_ID_THY[(c+100):(d+100),1]), 1,URL)
kegg_id_thy_3 <- apply(data.frame(CSID_ID_THY[(c+200):(d+200),1]), 1,URL)
kegg_id_thy_4 <- apply(data.frame(CSID_ID_THY[(c+300):(d+300),1]), 1,URL)
kegg_id_thy_5 <- apply(data.frame(CSID_ID_THY[(c+400):(d+400),1]), 1,URL)
kegg_id_thy_6 <- apply(data.frame(CSID_ID_THY[(c+500):(d+500),1]), 1,URL)
kegg_id_thy_7 <- apply(data.frame(CSID_ID_THY[(c+600):(d+600),1]), 1,URL)
kegg_id_thy_8 <- apply(data.frame(CSID_ID_THY[(c+700):(d+700),1]), 1,URL)
kegg_id_thy_9 <- apply(data.frame(CSID_ID_THY[(c+800):(d+800),1]), 1,URL)
kegg_id_thy_10 <- apply(data.frame(CSID_ID_THY[(c+900):(d+900),1]), 1,URL)
kegg_id_thy_11 <- apply(data.frame(CSID_ID_THY[(c+1000):nrow(CSID_ID_THY),1]), 1,URL)

kegg_id_thy <- rbind(
  data.frame(keggid = kegg_id_thy_1,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_2,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_3,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_4,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_5,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_6,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_7,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_8,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_9,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_10,stringsAsFactors = FALSE),
  data.frame(keggid = kegg_id_thy_11,stringsAsFactors = FALSE)
  )

kegg_id_df_thy <- data.frame(y1 = kegg_id_thy[,1], y2 = CSID_ID_THY[c:nrow(CSID_ID_THY),1], stringsAsFactors = FALSE)

colnames(kegg_id_df_thy) <- c("KEGG_ID", "CSID")

View(kegg_id_df_thy)
kegg_id_df_thy[kegg_id_df_thy$KEGG_ID == "Match error**",] #finding for match error

 ###########################
 #   Caution: pause here   #
 # Correct error maunually #
 ###########################

#------------------------/replace_match_error/-------------------
#Check and replace the match error data ID
#THY

#BROWSECS(834)
kegg_id_df_thy[120,1] <- "C03844"
kegg_id_df_thy[125,1] <- "C03365"
kegg_id_df_thy[198,1] <- "C03844"
kegg_id_df_thy[436,1] <- "C00158"
kegg_id_df_thy[505,1] <- "C01015"
kegg_id_df_thy[583,1] <- "C01231"
kegg_id_df_thy[654,1] <- "C03906"
kegg_id_df_thy[663,1] <- "C06473"
kegg_id_df_thy[756,1] <- "C06118"
kegg_id_df_thy[992,1] <- "C16439"

#kegg_id_df_thy[7,1] <- "C00158"
#kegg_id_df_thy[94,1] <- "C03906"
#kegg_id_df_thy[103,1] <- "C06473"
#kegg_id_df_thy[450,1] <- "C01231"
#kegg_id_df_thy[132,1] <- "C16439"

#kegg_id_df_thy[321,1] <- "C00204"
#kegg_id_df_thy[396,1] <- "C01015"
#kegg_id_df_thy[489,1] <- "C03844"
#kegg_id_df_thy[494,1] <- "C06353"
#kegg_id_df_thy[864,1] <- "C06118"

if(nrow(kegg_id_df_thy[kegg_id_df_thy$KEGG_ID == "Match error**",])>0){
  warning("Match error remaining!!")
}
#------------------------/replace_repeats/-----------------------

#removing repeated compounds
kegg_id_df_thy3 <- kegg_id_df_thy[!duplicated(kegg_id_df_thy$KEGG_ID),]
#removing "Drug only"
kegg_id_df_thy2 <- kegg_id_df_thy3[!kegg_id_df_thy3$KEGG_ID == "Drug only",]
write.csv(kegg_id_df_thy2, file = paste(getwd(), "/output/THY_COMP_KEGG.csv", sep = ""))

#----------------------*merging_my_list_to_total_list------------------------
#my compound list -> kegg_id_df_thy2
#all compound list -> compound_list
#my gene list -> THY_KEGGID
#all gene list -> gene_list

#------------------------/compound/----------------------------
colnames(compound_list) <- "KEGG_ID"
THY_COMP_KEGG3 <- data.frame(kegg_id_df_thy2$KEGG_ID)
colnames(THY_COMP_KEGG3) <- "KEGG_ID"
compound_list <- data.frame(compound_list)
test_comp <- merge(THY_COMP_KEGG3, compound_list)
#------------------------/gene/-------------------------------
colnames(THY_KEGGID) <- "EZ_ID"
THY_GENE_KEGG3 <- data.frame(THY_KEGGID)
colnames(gene_list) <- "EZ_ID"
gene_list <- data.frame(gene_list)
test_gene <- merge(THY_GENE_KEGG3, gene_list)
#------------------------/redirect/-------------------------------
kegg_id_df_thy2 <- test_comp
THY_KEGGID <- test_gene


#----------------------*count_list_num------------------------

#------------------------/compound/----------------------------
#calculation of my compounds
comp_num_thy <- nrow(kegg_id_df_thy2)
comp_num_thy
#184 compounds
#------------------------/gene/-------------------------------
#Calculating gene amount of my list
gene_num_thy <- nrow(THY_KEGGID) 
gene_num_thy
#My gene amount -> thy_gene_num
#584 genes
#------------------------/combine/-----------------------------
combine_num_thy <- comp_num_thy + gene_num_thy
combine_num_thy
#768 combine
#-----------------^^FISHER'S_EXACT_TEST^^-------------------


#----------------------*matching_database----------------------

library(KEGGREST)
library(dplyr)

#------------------------/compound/----------------------------
#matching all the compounds of each pathways to my compound list
match_y <- kegg_id_df_thy2$KEGG_ID
MATCHCOMPTOPATH <- function(x){
  match_x <- rownames(data.frame(unlist(PATHCOMPINFO2(x))))
  nrowoutput <- nrow(data.frame(intersect(match_x, match_y)))
  if(nrowoutput == 0){
    0
  }else{
    nrowoutput
  }
}
path_mscomp_match_num_thy <- apply(data.frame(path_mmu_id[a:b]), 1, MATCHCOMPTOPATH)
path_mscomp_match_num_df_thy <- data.frame(path_mscomp_match_num_thy)
#list of compound match number in each pathway -> path_mscomp_match_num_df_thy
#------------------------/gene/-------------------------------
match_z <- data.frame(THY_KEGGID)
match_z[,1] <- as.character(match_z[,1])
colnames(match_z) <- c("GENE")
MATCHGENETOPATH <- function(x){
  match_x <- data.frame(unlist(PATHGENEINFO2(x)))
  if(ncol(match_x) > 0){
    colnames(match_x) <- c("GENE")
    match_x[,1] <- as.character(match_x[,1])
    nrowoutput <- nrow(data.frame(intersect(match_x, match_z)))
    nrowoutput
  }else{
    0
  }
}
path_gene_match_num_thy <- apply(data.frame(path_mmu_id[a:b]), 1, MATCHGENETOPATH)
path_gene_match_num_df_thy <- data.frame(path_gene_match_num_thy)
#list of GENE match number in each pathway -> path_gene_match_num_df_thy
#------------------------/combine/-----------------------------

path_combine_match_num_df_thy <- data.frame(path_mscomp_match_num_df_thy + path_gene_match_num_df_thy)

#----------------------*fisher's_exact_test---------------------

#fisher's exact test calculation
#p_value
#------------------------/compound/----------------------------
#dhyper(x, m, n, k, log = FALSE)
#https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
FETCOMP <- function(x){
  X <- path_mscomp_match_num_df_thy[x,1]
  M <- comp_in_path_df[x,2]
  N <- comp_in_path_df_total_union - M
  K <- comp_num_thy
  A1 <- X
  A2 <- K-X
  B1 <- M-X
  B2 <- N-(K-X)
  T1 <- data.frame( y1 = c(A1, A2), y2 = c(B1, B2))
  T2 <- fisher.test(T1, alternative = "greater")
  fet_p_value <- T2$p.value
  fet_p_value
}
x_comp_thy <- data.frame(c(a:b))
p_value_path_thy1 <- format(data.frame(y1 = data.frame(path_mmu_id[a:b]),
                                      y2 = data.frame(path_mmu[a:b]),
                                      y3 = data.frame(comp_in_path_df[a:b,2]),
                                      y4 = path_mscomp_match_num_df_thy,
                                      y5 = apply(x_comp_thy, 1, FETCOMP),
                                      stringsAsFactors = FALSE),
                           scientific = FALSE)
colnames(p_value_path_thy1) <- c("PATHWAY_ID", "PATHWAY","TOTAL_COMP", "COMP_MATCH", "COMP_P_VALUE")
rownames(p_value_path_thy1) <- c(a:b)
#------------------------/gene/-------------------------------
FETGENE <- function(x){
  X <- path_gene_match_num_df_thy[x,1]
  M <- gene_in_path_df[x,2]
  N <- gene_in_path_df_total_union - M
  K <- gene_num_thy
  A1 <- X
  A2 <- K-X
  B1 <- M-X
  B2 <- N-(K-X)
  T1 <- data.frame( y1 = c(A1, A2), y2 = c(B1, B2))
  T2 <- fisher.test(T1, alternative = "greater")
  fet_p_value <- T2$p.value
  fet_p_value
}
x_gene_thy <- data.frame(c(a:b))
p_value_path_thy2 <- format(data.frame(y1 = data.frame(gene_in_path_df[a:b,2]),
                                       y2 = path_gene_match_num_df_thy,
                                       y3 = apply(x_gene_thy, 1, FETGENE),
                                       stringsAsFactors = FALSE),
                           scientific = FALSE)
colnames(p_value_path_thy2) <- c("TOTAL_GENE","GENE_MATCH", "GENE_P_VALUE")
rownames(p_value_path_thy2) <- c(a:b)
#------------------------/combine/-----------------------------
FETCOMBINE <- function(x){
  X <- path_combine_match_num_df_thy[x,1]
  M <- combine_in_path_df[x,2]
  N <- combine_in_path_df_total_union - M
  K <- combine_num_thy
  A1 <- X
  A2 <- K-X
  B1 <- M-X
  B2 <- N-(K-X)
  T1 <- data.frame( y1 = c(A1, A2), y2 = c(B1, B2))
  T2 <- fisher.test(T1, alternative = "greater")
  fet_p_value <- T2$p.value
  fet_p_value
}
x_combine_thy <- data.frame(c(a:b))
p_value_path_thy3 <- format(data.frame(y1 = data.frame(combine_in_path_df[a:b,2]),
                                       y2 = path_combine_match_num_df_thy,
                                       y3 = apply(x_combine_thy, 1, FETCOMBINE),
                                       stringsAsFactors = FALSE),
                            scientific = FALSE)
colnames(p_value_path_thy3) <- c("TOTAL_COMPONENT","COMPONENT_MATCH", "COMPONENT_P_VALUE")
rownames(p_value_path_thy3) <- c(a:b)


#-----------------^^OUTPUT^^----------------------------


OUTPUT <- data.frame(p_value_path_thy1, p_value_path_thy2, p_value_path_thy3)
View(OUTPUT)
write.csv(OUTPUT, file = paste(getwd(), "/output/THY_COMBINE.csv", sep = ""))
write.csv(kegg_id_df_thy2, file = paste(getwd(), "/output/THY_my_comp.csv", sep = ""))
write.csv(THY_KEGGID, file = paste(getwd(), "/output/THY_my_gene.csv", sep = ""))

 ########################
 #  END OF CALCULATION  #
 ########################
#-----------------^^Advanced_checking^^-----------------

#----------------------*input---------------------------
#------------------------/gene/-------------------


#------------------------/comp/-------------------
library(readr)
compound_measurements <- read_csv("compound measurements.csv")
compound_identification <- read_csv("compound identification.csv")
comp_feature <- compound_measurements[,c(1,3,5,8,14,34:51,55,57:59,61:63,65:67,69,71:72)]
comp_ident <- compound_identification[,1:2]
comp_total_list_feture_indent <- merge(comp_ident,comp_feature,by = "Compound")
comp_total_list_feture_indent <- comp_total_list_feture_indent[comp_total_list_feture_indent$`Maximum Abundance` >= 100,]
comp_total_list_feture_indent$CSID <- sub("*...D","",comp_total_list_feture_indent$`Compound ID`)

##check on chemspider for kegg id (only rerun when necessary)
#comp_feture_indent_1 <- apply(data.frame(comp_total_list_feture_indent$CSID[1:1000]), 1,URL) 
#comp_feture_indent_2 <- apply(data.frame(comp_total_list_feture_indent$CSID[1001:2000]), 1,URL) 
#comp_feture_indent_3 <- apply(data.frame(comp_total_list_feture_indent$CSID[2001:3000]), 1,URL) 
#comp_feture_indent_4 <- apply(data.frame(comp_total_list_feture_indent$CSID[3001:4000]), 1,URL) 
#comp_feture_indent_5 <- apply(data.frame(comp_total_list_feture_indent$CSID[4001:5000]), 1,URL) 
#comp_feture_indent_6 <- apply(data.frame(comp_total_list_feture_indent$CSID[5001:length(comp_total_list_feture_indent$CSID)]), 1,URL) 

comp_feature_indent <- rbind(
  data.frame(keggid = comp_feture_indent_1,stringsAsFactors = FALSE),
  data.frame(keggid = comp_feture_indent_2,stringsAsFactors = FALSE),
  data.frame(keggid = comp_feture_indent_3,stringsAsFactors = FALSE),
  data.frame(keggid = comp_feture_indent_4,stringsAsFactors = FALSE),
  data.frame(keggid = comp_feture_indent_5,stringsAsFactors = FALSE),
  data.frame(keggid = comp_feture_indent_6,stringsAsFactors = FALSE))

comp_total_list_feture_indent$KEGG_ID <- comp_feature_indent
rownames(comp_total_list_feture_indent) <- NULL

comp_error <- comp_total_list_feture_indent[comp_total_list_feture_indent$KEGG_ID == "Match error**",38:39]
comp_error
BROWSECS <- function(x){                                               #combining the url, note that the url is backstage of Chemspider
  url <- paste("http://www.chemspider.com/ibcontent.ashx?csid=",
               x,
               "&type=ds&ds_types=",
               sep = "")
  browseURL(url)
}

##for replacing metching error, extreme popout, rerun only when checking is updated
#for(n in 26:nrow(comp_error)){
#  BROWSECS(comp_error$CSID[n])
#  n <- n+1
#}

comp_total_list_feture_indent[21,39] <- "C01231"
comp_total_list_feture_indent[42,39] <- "C00158"
comp_total_list_feture_indent[56,39] <- "C03844"
comp_total_list_feture_indent[61,39] <- "C03365"
comp_total_list_feture_indent[144,39] <- "C00711"

comp_total_list_feture_indent[169,39] <- "C03906"
comp_total_list_feture_indent[178,39] <- "C06473"
comp_total_list_feture_indent[186,39] <- "C00663"
comp_total_list_feture_indent[204,39] <- "C01094"
comp_total_list_feture_indent[226,39] <- "C00103"

comp_total_list_feture_indent[295,39] <- "C11909"
comp_total_list_feture_indent[296,39] <- "C11922"
comp_total_list_feture_indent[337,39] <- "C00357"
comp_total_list_feture_indent[418,39] <- "C00204"
comp_total_list_feture_indent[611,39] <- "C05411"

comp_total_list_feture_indent[623,39] <- "C00257"
comp_total_list_feture_indent[631,39] <- "exclude"
comp_total_list_feture_indent[722,39] <- "C00341"
comp_total_list_feture_indent[723,39] <- "C02569"
comp_total_list_feture_indent[1044,39] <- "C06035"

comp_total_list_feture_indent[1072,39] <- "C00897"
comp_total_list_feture_indent[1112,39] <- "C06869"
comp_total_list_feture_indent[1121,39] <- "C00897"
comp_total_list_feture_indent[1307,39] <- "C15586"
comp_total_list_feture_indent[1322,39] <- "C16439"

comp_total_list_feture_indent[1366,39] <- "C00897"
comp_total_list_feture_indent[1960,39] <- "C20322"
comp_total_list_feture_indent[2048,39] <- "C15583"
comp_total_list_feture_indent[2061,39] <- "C00805"
comp_total_list_feture_indent[2107,39] <- "C06677"

comp_total_list_feture_indent[2176,39] <- "C02061"
comp_total_list_feture_indent[2404,39] <- "C03862"
comp_total_list_feture_indent[2892,39] <- "C05665"
comp_total_list_feture_indent[3135,39] <- "C06341"
comp_total_list_feture_indent[3191,39] <- "C04483"

comp_total_list_feture_indent[3523,39] <- "C14826"
comp_total_list_feture_indent[3744,39] <- "C18131"
comp_total_list_feture_indent[3777,39] <- "C14826"
comp_total_list_feture_indent[3934,39] <- "exclude"
comp_total_list_feture_indent[3961,39] <- "C00059"

comp_total_list_feture_indent[4197,39] <- "C06123"
comp_total_list_feture_indent[4327,39] <- "C00059"
comp_total_list_feture_indent[4539,39] <- "C01530"
comp_total_list_feture_indent[4569,39] <- "C01530"
comp_total_list_feture_indent[4644,39] <- "C00897"

comp_total_list_feture_indent[4953,39] <- "C00180"
comp_total_list_feture_indent[5044,39] <- "C01530"
comp_total_list_feture_indent[5305,39] <- "exclude"
comp_total_list_feture_indent[5354,39] <- "C17435"
comp_total_list_feture_indent[5548,39] <- "C06033"

comp_total_list_feture_indent[5578,39] <- "C00059"
comp_total_list_feture_indent[5665,39] <- "C02704"
comp_total_list_feture_indent[5673,39] <- "C06677"
comp_total_list_feture_indent[5730,39] <- "C19208"

comp_error <- comp_total_list_feture_indent[comp_total_list_feture_indent$KEGG_ID == "Match error**",38:39]
if(nrow(comp_error) > 0 ){
  warning("match error remaining!")
}
compound_feature_kegg <- comp_total_list_feture_indent[comp_total_list_feture_indent$KEGG_ID != "Drug only",]
rownames(compound_feature_kegg) <- NULL

chosen_compouund_ident <- read_csv("chosen compouund identification.csv")
input_already_compound <- read_csv("20190524 comp_MF_THY_FC1.5_P0.05.csv", 
                                   col_names = FALSE)
colnames(input_already_compound)[1] <- "CSID"

my_comp_info_1 <- merge(chosen_compouund_ident,input_already_compound)
colnames(my_comp_info_1) <- c("Compound ID","Compound","CSID")
my_comp_info_2 <- merge(kegg_id_df_thy,my_comp_info_1)
my_comp_info_3 <- compound_feature_kegg

#----------------------*using funtion-----------------------------
#------------------------/gene/-------------------
gene_duplicate_check5 <- gene_duplicate_check4
colnames(gene_duplicate_check5)[16] <- "EZ_ID"
my_gene_info <- merge(THY_KEGGID, gene_duplicate_check5, by = "EZ_ID")
#for gene advance information, check my_gene_info
#------------------------/comp/-------------------

match_a2 <- my_comp_info_2  
match_a <- my_comp_info_3
colnames(match_a)[39] <- "keggid"
colnames(match_a2)[2] <- "keggid"

#showing all detected features that match to given pathway
MATCHCOMPTOPATH2 <- function(x){
  output <- suppressWarnings(data.frame(intersect(data.frame(keggid = rownames(data.frame(unlist(PATHCOMPINFO2(x))))), 
                                                  match_a$keggid)))
  
  if(nrow(output) == 0){
    0
  }else{
    output
  }
}
#showing all chosen features and its signal source that match to given pathway
COMPINFO <- function(x){
  output <- merge(data.frame(keggid = rownames(data.frame(unlist(PATHCOMPINFO2(x))))), 
                               data.frame(keggid = match_a2$keggid), all = FALSE)
  output2 <- merge(output,match_a2, all = FALSE)
  output3 <- unique(output2[order(output2$Compound),c(1,2,4)])
  rownames(output3) <- NULL
  return(output3)
  }

pATHCOMPINFO2 <- function(x){
  e <- keggGet(x)
  f <- list(e[[1]]$COMPOUND)
  f
}
mATCHCOMPTOPATH <- function(x){
  match_x <- rownames(data.frame(unlist(pATHCOMPINFO2(x))))
  match_x <- data.frame(match_x)
  colnames(match_x) <- "keggid"
  m <- merge(match_x, match_a)
  if(nrow(m) == 0){
    0
  }else{
    distinct(m)
  }
}
hEATMAPCOMP <- function(x){
  b <- mATCHCOMPTOPATH(x)
  b_0 <- b[,3:14]
  
  if(nrow(b_0) < 2){
    warning("Pathway compound match Null")
  }else{
    row.names(b_0) <- b$Description
    b_0 <- distinct(b_0)
    c <- data.matrix(b_0)
    h <- heatmap(c,
                 Colv = NA,
                 scale = "row",
                 col =  brewer.pal(9,"Blues"),
                 main = paste(x,keggGet(x)[[1]]$NAME, sep = "  "))
    h
    b[,c(1,2,18,16,15)]
  }
}
hEATMAPCOMP2 <- function(x){
  b <- mATCHCOMPTOPATH(x)
  b_0 <- b[,3:14]
  
  if(nrow(b_0) < 2){
    warning("Pathway compound match Null")
  }else{
    b_0 <- distinct(b_0)
    rownames(b_0) <- unique(b$`Compound`)
    c <- data.matrix(b_0)
    h <- heatmap(c,
                 Colv = NA,
                 scale = "row",
                 col =  brewer.pal(9,"Reds"),
                 main = paste(x,keggGet(x)[[1]]$NAME, sep = "  "))
    h
    b[,c(1,2,18,16,15)]
  }
}

#NOt finished, I thought drawing a compound heatmap might not make sense
#HEATMAPCOMP <- function(x){
#  b <- merge(data.frame(keggid = MATCHCOMPTOPATH2(x)), match_a, all = FALSE)
#  
#  if(nrow(b_0) < 2){
#   warning("Pathway compound match Null")
#  }else{
#    row.names(b_0) <- b$Description
#    b_0 <- distinct(b_0)
#    c <- data.matrix(b_0)
#    h <- heatmap(c,
#                Colv = NA,
#                 scale = "row",
#                 col =  brewer.pal(9,"Blues"),
#                 main = paste(x,keggGet(x)[[1]]$NAME, sep = "  "))
#    h
#    b[,c(1,2,18,16,15)]
#  }
#}
#----------------------*Run-----------------------------
message("for advance information, check my_gene_info or use the COMPINFO(single) or CHECKCOMP(multiple) function!")
View(my_gene_info)

sig_path <- OUTPUT_order[which(OUTPUT_order$`COMPONENT_P_VALUE` <= 0.05 
                               & OUTPUT_order$`COMP_MATCH` != " 0" 
                               & OUTPUT_order$`GENE_MATCH` != " 0"), 1] #AsIs data class not dealed
CHECKCOMP <- function(j){
  pathk <- length(j)
  if(pathk <= 25){
    x <-1
    for(x in 1:pathk){
      print(j[x])
      print(COMPINFO(j[x]))
      print("-----------------------------------")
            x <- x + 1
    }
  }else if(pathk > 25){
    warning("Pathway should be less than 25")
  }else{
    warning("Error")
  }
}
CHECKCOMP(sig_path)
 ########################
 #  END OF CALCULATION  #
 ########################
#-----------------^^Visualization^^---------------------------
#----------------------*constant----------------------
match_y <- kegg_id_df_thy2$KEGG_ID
match_z <- data.frame(THY_KEGGID)
match_z[,1] <- as.character(match_z[,1])
colnames(match_z) <- c("GENE")
VISUALCOMPTOPATH <- function(x){
  match_x <- rownames(data.frame(unlist(PATHCOMPINFO2(x))))
  match_xy <- intersect(match_x, match_y)
  match_xy
}
VISUALGENETOPATH <- function(x){
  match_x <- data.frame(unlist(PATHGENEINFO2(x)))
  colnames(match_x) <- c("GENE")
  match_x[,1] <- as.character(match_x[,1])
  match_xx <-intersect(match_x, match_z)
  match_xx[,1]
}
VISUAL <- function(x){
  id <- x
  visual_comp <- VISUALCOMPTOPATH(id)
  visual_gene <- VISUALGENETOPATH(id)
  object <- c(visual_gene, visual_comp)
  url <- suppressMessages(color.pathway.by.objects(id , object,
                                                   fg.color.list = c( rep("red", length(object))),
                                                   bg.color.list = c( rep("black", length(object)))
                                                  )
                          )
  url
}
#----------------------*mark_in_path----------------------
BROWSEPATH <- function(j){
  pathk <- length(j)
  if(pathk <= 25){
  y <- 0
for(x in c(1:pathk)){
    browseURL(VISUAL(j[x]))
    y <- x + y
  }
}else if(pathk > 25){
  warning("Pathway should be less than 25")
}else{
  warning("Error")
  }
}
#----------------------*KEGG_for_significant_pathway_id----------------
#call for significant combine pathways, with GENE_MATCH and COMP_MATCH > 0
OUTPUT_order <- OUTPUT[order(OUTPUT$COMPONENT_P_VALUE),]
sig_path <- OUTPUT_order[which(OUTPUT_order$`COMPONENT_P_VALUE` <= 0.05 
                         & OUTPUT_order$`COMP_MATCH` != " 0" 
                         & OUTPUT_order$`GENE_MATCH` != " 0"), 1] #AsIs data class not dealed
BROWSEPATH(sig_path)
#----------------------*KEGG_for_pathway_id----------------
#pentose phosphate pathway
BROWSEPATH("path:mmu00030")
#pyrimidine metabolism
BROWSEPATH("path:mmu00240")
BROWSEPATH("path:mmu05169")

########################
#  END OF CALCULATION  #
########################