### HPRD database
pathways = read.csv("HPRD/HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt", head=FALSE, sep="\t")
pathways=pathways[,c(1, 7,4)]
colnames(pathways) = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
hprd_dir=G
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="undirected")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
hprd_undir=G

### MINT database

pathways = read.csv("MINT/2013-03-26-mint-human-binary.mitab26.txt",head=TRUE, sep="\t")
pathways=pathways[,c(1, 12,2)]
colnames(pathways) = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
mint_dir=G
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="undirected")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
mint_undir=G

### Intact database

load("allIntactDataFrame.rda")
#pathways = read.csv("Intact/psimitab/intact.txt", head=TRUE, sep="\t")
#pathways=pathways[,c(1, 15,2)]
#node=pathways[,1]
#x=grep(pattern="uniprotkb:", node)
#y=as.character(node[x])
#newintact=sapply(strsplit(y, ":"), "[",2)
#z=c()
#for(i in 1:length(newintact)){
#	if(length(intraIDMapper(newintact[i], species="HOMSA", srcIDType="UNIPROT",destIDType="EG"))!=0)
#		z[i]=intraIDMapper(newintact[i], species="HOMSA", srcIDType="UNIPROT",destIDType="EG")
#	else
#		z[i]="Unknown"
#}


#colnames(pathways)[1:3] = c("GeneA", "INTERACTION_TYPE", "GeneB")
pathways=entIntact
pathways=pathways[complete.cases(pathways),]
pathways <- pathways[!(pathways$GeneA == ""), ]
pathways <- pathways[!(pathways$GeneB == ""), ]

G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
intact_entdir=G

G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="undirected")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
intact_entundir=G


###  Biogrid database
pathways=read.csv("/home/praveen/Paurush/COSBI/interactomes/Alldata/BIOGRID-ORGANISM-Homo_sapiens-3.2.110.tab2.txt", head=TRUE, sep="\t")
pathways = pathways[,c(2, 12, 3)] # For Entrez IDs
colnames(pathways) = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)

biogridEnt_dir=G

G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="undirected")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
biogridEnt_undir=G

pathways=read.csv("/home/praveen/Paurush/COSBI/interactomes/Alldata/homo_sapiens.interactions.txt", head=FALSE, skip=1,sep="\t")
pathways = pathways[,c(8, 12, 9)] #for Symbols
colnames(pathways) = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
biogridSym_dir=G

G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="undirected")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
biogridSym_undir=G



###  Reactome database
pathways=read.csv("/home/praveen/Paurush/COSBI/interactomes/Alldata/homo_sapiens.interactions.txt", head=TRUE, sep="\t")
pathways = pathways[,c(2, 7, 5)] # For ENSEMBL IDs
colnames(pathways) = c("GeneA", "INTERACTION_TYPE", "GeneB")
pathways <- pathways[!(pathways$GeneA == ""), ]
pathways <- pathways[!(pathways$GeneB == ""), ]
genes=sapply(strsplit(as.character(pathways$GeneB), ":"), "[[", 2)
pathways$GeneB=genes
gene=pathways$GeneA
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = gene, mart = mart)
genes=results$hgnc_symbol
names(genes)=results$ensembl_gene_id
gene=genes[pathways$GeneA]
pathways$GeneA=as.character(gene)


G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)

Reactome_entdir=G

G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="undirected")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
Rectome_entundir=G
gg=map_reactome[pathways$GeneA]

> gg=map_reactome[pathways$GeneA]
> length(gg)
[1] 1044417
> dim(pathways)
[1] 1044417       3
> pathways$GeneA=gg
> gg=map_reactome[pathways$GeneB]
> pathways$GeneB=gg


###  String database
pathways=read.csv("/home/praveen/Paurush/COSBI/interactomes/Alldata/STRING/9606.protein.links.v9.1.txt", head=TRUE, sep=" ")
pathways = pathways[,c(1, 3, 2)] # For ENSEMBL IDs
colnames(pathways) = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)

String_Entdir=G

G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="undirected")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
G = removeSelfLoops(G)
String_Entundir=G

##########
> genes=sapply(strsplit(genes, "[.]"), "[[", 2)
> nodes(String_undir)=genes
> genes=nodes(String_dir)
> genes=sapply(strsplit(genes, "[.]"), "[[", 2)
> nodes(String_dir)=genes


library(biomaRt)
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene=as.character(pathways[,"GeneB"])
gene=sapply(strsplit(gene, "[.]"), "[[", 2)
results <- getBM(attributes = c("ensembl_peptide_id", "hgnc_symbol"), filters = "ensembl_peptide_id", values = gene, mart = mart)
g=results$hgnc_symbol
names(g)=results$ensembl_peptide_id
pathways$GeneA=g[pathways$GeneA]
pathways$GeneB=g[pathways$GeneB]
gene
g=results$hgnc_symbol
names(g)=results$ensembl_peptide_id
class(g)
gg=g[gene]
pathways$GeneB=as.character(gg)
head(pathways)

pathways=pathways[complete.cases(pathways),]
pathways <- pathways[!(pathways$GeneA == ""), ]
pathways <- pathways[!(pathways$GeneB == ""), ]


> gg2=gg
> mapname=mapnames[,2]
> names(mapname)=mapnames[,1]
> gg[names(wich(is.na(gg)))]=mapname[names(wich(is.na(gg)))]
Error: could not find function "wich"
> gg[names(which(is.na(gg)))]=mapname[names(which(is.na(gg)))]


String_Entdir=String_Symdir
String_Entundir=String_Symundir
##########
load("allnet.rda")
# [1] "biogridEnt_dir"   "biogridEnt_undir" "biogridSym_dir"   "biogridSym_undir"
# [5] "hprd_dir"         "hprd_undir"       "intact2"          "intact_dir"      
# [9] "intact_undir"     "mint_dir"         "mint_undir"       "Reactome_dir"    
#[13] "Rectome_undir"    "String_dir"       "String_undir"  



library(org.Hs.eg.db)
gg=mget(as.character(genesjan), org.Hs.egSYMBOL, ifnotfound=NA)
gg=as.character(unlist(gg))


library(org.Hs.eg.db)
gg=mget(as.character(nodes(rec)),  org.Hs.egALIAS2EG, ifnotfound=NA)
gg=as.character(unlist(gg))

x=names(which(is.na(gg)))
for(i in 1:length(x)){
gg[x]=map_String[[x]][1]
}

gg[which(is.na(gg))]=map_String[x]
> map_String=gg
> pathways$GeneA=gg[pathways$GeneA]
> pathways$GeneB=gg[pathways$GeneB]
> head(pathways)
  GeneA INTERACTION_TYPE GeneB
1   381              176 51507
2   381              327 79931
3   381              718  3835
4   381              272 57106
5   381              241  1738
6   381              170  9342






> allgraphsEnt=list(BiogridDir=biogridEnt_dir, BiogridUndir=biogridEnt_undir, KeggDir=G.KEGG, HprdDir=hprd_dir_ent, HprdUndir=hprd_undir_ent,IntactDir=intact_entdir,IntactUndir=intact_entundir, ReactomeDir=Reactome_entdir, ReactomeUndir=Rectome_entundir, StringDir=String_Entdir, StringUndir=String_Entundir)
> allgraphsSym=list(BiogridDir=biogridSym_dir, BiogridUndir=biogridSym_undir, KeggDir=G.KEGG_sym, HprdDir=hprd_dir, HprdUndir=hprd_undir,IntactDir=intact_symdir,IntactUndir=intact_symundir, ReactomeDir=Reactome_symdir, ReactomeUndir=Rectome_symundir, StringDir=String_Symdir, StringUndir=String_Symundir)
> nodes(hprd_dir)[1:10]
 [1] "ALDH1A1" "ITGA7"   "PPP1R9A" "SRGN"    "GRB7"    "PAK1"    "DLG4"   
 [8] "PIK3R2"  "PTPN18"  "ERBB2IP"
> save(allgraphsEnt, file="entrezGraphs.rda")
> save(allgraphsSym, file="SymbolGraphs.rda")
> allmaps=list(geneMap_HPRD, map_reactome, map_String)
> save(allmaps, file="nameMaps.rda")


## Get Uniprot IDs

mget(as.character(nodes(graphs$Biogrid)[1:10]),  org.Hs.egUNIPROT, ifnotfound=NA)

allBG_UP=mget(as.character(nodes(graphs$Biogrid)), org.Hs.egUNIPROT, ifnotfound=NA)
allHPRD_UP=mget(as.character(nodes(graphs$HPRD)), org.Hs.egUNIPROT, ifnotfound=NA)
allST_UP=mget(as.character(nodes(graphs$STRING)), org.Hs.egUNIPROT, ifnotfound=NA)
allIN_UP=mget(as.character(nodes(graphs$IntAct)), org.Hs.egUNIPROT, ifnotfound=NA)

graphs$Biogrid
length(allBG_UP)
graphs$HPRD
length(allHPRD_UP)
graphs$STRING
length(allST_UP)
graphs$IntAct
length(allIN_UP)

max(unlist(lapply(allIN_UP, length)))

### Conversion code for entrez graph to uniprot graph
### Converts EID to UProt and joins all the graphs


joinAll = function(myGraphs){
  g = myGraphs[[1]]
  for(i in 2:length(myGraphs)){
    g = join(g, myGraphs[[i]])
  }
  return(g)
}


map2uniprot=function(myGraph){
  allNodes_UP = mget(as.character(nodes(myGraph)),  org.Hs.egUNIPROT, ifnotfound=NA)
  nof = max(unlist(lapply(allNodes_UP, length)))
  listOfGraphs = list()
  for(k in 1:nof){
    listOfGraphs[[k]]=myGraph
  }
  for(i in 1:length(nodes(myGraph))){
    UP_here=allNodes_UP[[i]]
    for(j in 1:length(UP_here)){
      if(!is.na(allNodes_UP[[i]][j])& !(allNodes_UP[[i]][j] %in% nodes(listOfGraphs[j]))){
	nodes(listOfGraphs[[j]])[i]=UP_here[j]
      }else{
	nodes(listOfGraphs[[j]])[i]=nodes(listOfGraphs[[j]])[i]
      }  
    }
  }
  newGraph = joinAll(listOfGraphs)
  return(newGraph)
}

  
#     for(j in 1:length(allNodes_UP[[i]])){
# 	print(nodes(listOfGraphs[[j]])[i])
# 	cat("\n")
# 	print(allNodes_UP[[i]][[j]])
# 	if(!is.na(allNodes_UP[[i]][[j]])){
# 	  nodes(listOfGraphs[[j]])[i] = allNodes_UP[[i]][[j]]
# 	}else{
# 	  nodes(listOfGraphs[[j]])[i]=nodes(listOfGraphs[[j]])[i]
# 	}
#     }
#   }
#   newGraph = joinAll(listOfGraphs)
#   return(newGraph)
# }






for(i in 1:length(nodes(myGraph))){
    UPid=mget(as.character(nodes(myGraph)[i]),  org.Hs.egUNIPROT, ifnotfound=NA)
    for(j in 1:length(unlist()))
  }




ConsMat=function(ids, graphList){
  mat=matrix(0, nrow=length(ids), ncol=length(ids))
  colnames(mat)=rownames(mat)=as.character(ids)
 
}


computeSP=function(v,g){
matSP=matrix(NA, nrow=length(v), ncol=length(v))
colnames(matSP)=rownames(matSP)=v
  for(x in 1:length(v)){
    for(y in 1:length(v)){
	matSP[x,y]=sp.between(g,start=v[x],finish=v[y], detail=TRUE)[[1]]$length
    }
  }
  return(matSP)
}


runTest=function(allgenes2, graphs2){
mat2=list()
for(i in 1:length(allgenes2)){
  mat2[[i]]=list()
  for(j in 1:length(graphs2)){
    nodes_here=as.character(allgenes2[[i]])
    mat2[[i]][[j]]=matrix(NA, ncol=length(nodes_here),nrow=length(nodes_here))
    colnames(mat2[[i]][[j]])=rownames(mat2[[i]][[j]])=nodes_here
    nodes_here=intersect(nodes(graphs2[[j]]), nodes_here)
    cat("Nodes found= ", length(nodes_here))
    cat("\n")
    cat("Graph Running: ")
    cat(names(graphs2)[[j]])
    cat("\n")
    cat("Genes Running: ")
    cat(names(allgenes2)[[i]])
    cat("\n")
    mat2[[i]][[j]][nodes_here, nodes_here]=computeSP(nodes_here, graphs2[[j]])
  }
}
return(mat2)
}


DIR="/home/praveen/Paurush/COSBI/interactomes/Results/EvalInp"
filenames <- list.files(DIR, pattern="*.tab", full.names=TRUE)
setwd(DIR)
allgenes=list()
categ=c()
for(i in 1:length(filenames)){
temp=read.csv(filenames[i], sep="\t", head=TRUE)
tempname=strsplit(filenames[i], ":")[[1]][2]
tempname=strsplit(tempname, "[.]")[[1]][1]
allgenes[[i]]=temp[,1]
categ=c(categ,tempname)
rm(temp, tempname)
}
names(allgenes)=categ






graphs2=graphs[c(2,5,7,11)]
allgenes2=allgenes[c(1,11,10,9)]



# A venn Diagram
library(VennDiagram)
pdf("venn.pdf")
venn.diagram(x = list(Biogrid = nodes(graphs2[[1]]),HPRD = nodes(graphs2[[2]]),IntAct = nodes(graphs2[[3]]),STRING = nodes(graphs2[[4]])),
resolution =500, col = "transparent", fill = c("cornflowerblue","green","yellow","darkorchid1"),
alpha = 0.50, label.col = c("orange", "white", "darkorchid4", "white", "white", 
"white",    "white", "white", "darkblue", "white", "white", "white", "white", 
"darkgreen", "white"), cex = 1.5, fontfamily = "serif", fontface = "bold",
cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5,
cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", rotation.degree = 270,
margin = 0.2, imagetype="svg")
dev.off()

venn.plot = venn.diagram(x = list(Biogrid = nodes(graphs2[[1]]),HPRD = nodes(graphs2[[2]]),IntAct = nodes(graphs2[[3]]),STRING = nodes(graphs2[[4]])), filename = NULL, col = "transparent", fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),alpha = 0.50,c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),cex = 1.5,main="Sharing of proteins across PPI databases",fontfamily = "serif",fontface = "bold",cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),cat.cex = 1.5,cat.pos = 0,cat.dist = 0.07,cat.fontfamily = "serif",rotation.degree = 270,margin = 0.2)

pdf("venn-genes.pdf")
grid.draw(venn.plot)
dev.off()

ppidbInfo=data.frame(names=c("Biogrid","HPRD","IntAct","STRING"), 
Proteins=c(length(nodes(graphs2[[1]])),Nodes=length(nodes(graphs2[[2]])),length(nodes(graphs2[[3]])),length(nodes(graphs2[[4]]))),
Edges=c(129740, 37039, 24311, 2100237))

ppidbInfo=melt(ppidbInfo)
colnames(ppidbInfo)=c("Database", "Property", "Number")
pdf("dbstat.pdf")
ggplot(ppidbInfo, aes(x=Database, y=Number, fill="red"))+geom_bar(position="dodge", stat="identity")+facet_wrap(~Property, scales="free_y")
dev.off()



Nodes found=  29
Graph Running: BiogridUndir
Genes Running: Citric
Nodes found=  23
Graph Running: HprdUndir
Genes Running: Citric
Nodes found=  23
Graph Running: IntactUndir
Genes Running: Citric
Nodes found=  30
Graph Running: StringUndir
Genes Running: Citric


Nodes found=  32
Graph Running: BiogridUndir
Genes Running: RNA
Nodes found=  24
Graph Running: HprdUndir
Genes Running: RNA
Nodes found=  25
Graph Running: IntactUndir
Genes Running: RNA
Nodes found=  30
Graph Running: StringUndir
Genes Running: RNA
Nodes found=  151


Graph Running: BiogridUndir
Genes Running: Ribosomal
Nodes found=  79
Graph Running: HprdUndir
Genes Running: Ribosomal
Nodes found=  139
Graph Running: IntactUndir
Genes Running: Ribosomal
Nodes found=  145
Graph Running: StringUndir
Genes Running: Ribosomal


ed=list()
for(i in 1:length(graphs2)){
f=edgeL(graphs2[[i]])
f=unlist(f)
nam=sapply(strsplit(names(f), "[.]"), "[[",1)
ed[[i]]=paste(nam,f, sep="->")
}


venn.plot2 = venn.diagram(x = list(Biogrid = ed[[1]],HPRD = ed[[2]],IntAct = ed[[3]],STRING = ed[[4]]), filename = NULL, col = "transparent", fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),alpha = 0.50,c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),cex = 1.5,main="Sharing of edges across PPI databases",fontfamily = "serif",fontface = "bold",cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),cat.cex = 1.5,cat.pos = 0,cat.dist = 0.07,cat.fontfamily = "serif",rotation.degree = 270,margin = 0.2)
pdf("venn-ed.pdf")
grid.draw(venn.plot2)
dev.off()



fuzzyCluster=function(X){
    if (missing(X)) 
        stop("The data set must be given")
    if (is.null(X)) 
        stop("The data set X is empty")
    X = as.matrix(X)
    if (any(is.na(X))) 
        stop("The data set X must not contain NA values")
    if (!is.numeric(X)) 
        stop("The data set X is not a numeric data.frame or matrix")
#     cat(" ", fill = TRUE)
#     cat("WELCOME to the interactive Fclust program", fill = TRUE)
#     cat("Warning: If you insert an object of mode CHARACTER when not requested, an error occurs and the program stops!")
#     cat(" ", fill = TRUE)
    n = nrow(X)
    startU = NULL
    k = ceil(n/2)
    m = 2
    ent = 1
    vp = rep(1, k)
    RS = 1
    conv = 1e-09
    maxit = 1e+6
    stand = 1
    delta = NULL
    ana = 1
    ok = 1
    param = 0
    while (ana == 1) {
#        
            clust = FKM(X, k, m, RS, stand, startU, conv, maxit)
#         
        ana = scan("", n = 1)
        if (length(ana) == 0) {
            ana = 0
        }
        if (ana == 1) {
            cat(" ", fill = TRUE)
            cat("If you want to run the same clustering algorithm, specify 1: ", 
                fill = TRUE)
            ok = scan("", n = 1)
            if (length(ok) == 0) {
                ok = 0
            }
            if (ok != 1) {
                cat(" ", fill = TRUE)
                cat("You chose to run a new clustering algorithm", 
                  fill = TRUE)
                cat("If you want to use the same parameter values, specify 1: ", 
                  fill = TRUE)
                param = scan("", n = 1)
                if (length(param) == 0) {
                  param = 0
                }
            }
            else param = 0
#              
        }
    }
    return(clust)
}

cc=fuzzyCluster(Mc[,1:(ncol(Mc)-1)])
cc[[1]]
