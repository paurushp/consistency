# KEGG simulation 1 (A fixed numer of nodes) and 2 (an entire pathway)
sampleGraph = function(n){
set.seed(123456789)
	load("KEGGgraphs.rda")
	pdf(file.path(outputdir, paste("sampledNetworks_",n,"Sgenes.pdf")))
	subnets = list()
	i = 0
	while(i < 10){
		cand.net = sampleKEGGPathway(graphs, nacc, n)
		cand.net.mat = as(cand.net, "matrix")
		exists.net = any(sapply(subnets, function(subn) all(subn == cand.net.mat))) # CAUTION: This is not an exact test for graph isomorphism!!!
		if(!exists.net && (i < 9 || (i == 9 && length(tsort(cand.net)) == 0))){ # at least one graph has to have a cycle			
			i = i + 1			
			plot(cand.net)
			subnets[[i]] = cand.net.mat			
		}
	}
	save(subnets, file=file.path(outputdir, paste("sampledNetworks_",n,"Sgenes.rda")))
	dev.off()	
	return(G.samp=subnets)
}

#############################################################################################
#############################################################################################
#############################################################################################
pathways2 = read.csv("rupi_metacoreCur.txt", sep="\t", header=F)
colnames(pathways2) = c("Source", "Target", "INTERACTION")
pathways2 = pathways2[pathways2[,3] %in% c("Activation", "Inhibition"), ]
G2 = new("graphNEL", nodes=unique(c(as.character(pathways2$Source), as.character(pathways2$Target))), edgemode="directed")
G2 = addEdge(as.character(pathways2$Source), as.character(pathways2$Target), G2)
G2

paths=list()
paths2=list()
source("sp.path.compute.R")
for(i in 1:length(graphs)){
print(i)
paths[[i]] = sp.path.compute(genelist=c(nodes(graphs[[i]])), G=graphs[[i]])
paths2[[i]] = paths[sapply(paths[[i]], function(p) !is.na(p$length))]
}
names(paths)=names(paths2)=names(graphs)
testgraph=sampleGraph(n=c(5,10,20,30,40,50))
for(i in 1:length(tau)){
	G	
	Valid_Lit=validateLit(GraphObj=G, GG=graphs[[i]])
	allNone=sum(Valid_Lit$explained.by=="none")
	allNF=sum(Valid_Lit$explained.by=="NF")
	allEnt=nrow(Valid_Lit)-allNF
	allExp=nrow(Valid_Lit)-(allNone+allNF)
	Frac1=(allExp/allEnt)*100
	cat("Validating Literature against Model \n")
	Valid_Mode=validateMod(GraphObj=G, Allpaths=paths2[[i]])
	allNone=sum(Valid_Mode$explained.by=="none")
	allNF=sum(Valid_Mode$explained.by=="NF")
	allEnt=nrow(Valid_Mode)-allNF
	allExp=nrow(Valid_Mode)-(allNone+allNF)
	Frac2=(allExp/allEnt)*100
	exp=rbind(exp, data.frame(Threshold=tau[i], Model.View=Frac1,Know.View=Frac2))
	write.csv(Valid_Lit, file=fileL)
	write.csv(Valid_Mode, file=fileM)
}

ggplot(exp, aes(x=Threshold, y=))

#########################################

validateLit=function(GraphObj, GG=G2){
gh=GraphObj
em = edgeMatrix(gh)
eval = data.frame()
if(length(em)>0){
for(i in 1:NCOL(em)){
	cand.path = paste(nodes(gh)[em["from", i]], "->",nodes(gh)[em["to", i]], sep="")
	#print(nodes(gh)[em["from", i]])
	#print(nodes(gh)[em["to", i]])
	if((nodes(gh)[em["to", i]] %in% nodes(GG)) & (nodes(gh)[em["from", i]] %in% nodes(GG))){
		path = sp.between(GG, nodes(gh)[em["from", i]], nodes(gh)[em["to", i]])
	if(!is.na(path[[1]]$length))
		eval = rbind(eval, data.frame(edge=cand.path, explained.by=paste(path[[1]]$path_detail, collapse="->")))
	else
		eval = rbind(eval, data.frame(edge=cand.path, explained.by="none"))
	}
	else
		eval = rbind(eval, data.frame(edge=cand.path, explained.by="none"))
}
}
return(eval)
}


#########################################
validateMod=function(GraphObj, Allpaths=paths2){
G=GraphObj
eval2 = data.frame()
for(p in 1:length(Allpaths)){
	mypath = names(Allpaths[[p]])
	mynodes = strsplit(mypath, ":")[[1]]
	#print(mynodes)
	if((mynodes[1] %in% nodes(G)) & (mynodes[2] %in% nodes(G))){
	path = sp.between(G, mynodes[1],  mynodes[2])
	if(!is.na(path[[1]]$length)){
		eval2 = rbind(eval2, data.frame(literature.path=paste(mynodes[1], "->", mynodes[2], sep=""), explained.by=paste(path[[1]]$path_detail, collapse="->")))
	}
	else{
		eval2 = rbind(eval2, data.frame(literature.path=paste(mynodes[1], "->", mynodes[2], sep=""), explained.by="none"))
	}
	
}
else{
		eval2 = rbind(eval2, data.frame(literature.path=paste(mynodes[1], "->", mynodes[2], sep=""), explained.by="NF"))
	}
	}
return(eval2)
}


exp2=data.frame(Threshold=rep(exp2[c(1:10),2],2), exp2[11:30,c(1,2)])

ggplot(data=exp2, aes(x=Threshold, y=Fraction, fill=Measure)) + geom_bar(stat="identity", position=position_dodge())
#############################################################################################
#############################################################################################
#############################################################################################
load("SymbolGraphs.rda")
#graphs=list(HPRD=hprd,INTACT=intactSym, BIOGRID=biogridSym)
graphs=allgraphsSym
# Simulation 2 entire database
testConsis=function(graphs, iter=1, size=c(200,500, 1000)){
	d=list()
	vert=nodes(graphs[[1]])	
	for(i in 2:length(graphs)){
		vert=intersect(vert, nodes(graphs[[i]]))
	}
	for(j in 1:iter){
		testNodes=sample(vert, size, replace=FALSE)
		dt=data.frame()
		for(k in 1:length(graphs)){
			dt = rbind(dt, c(as(subGraph(testNodes, graphs[[k]]), "matrix")))
		}
		dt=t(dt)
		colnames(dt) = names(graphs)
		n1=as(subGraph(testNodes, graphs[[k]]), "matrix")
		n=c()		
		for(a in 1:ncol(n1)){
			for(b in 1:nrow(n1)){
			n=c(n, paste(colnames(n1)[a], rownames(n1)[b], sep="->"))			
			}
		}
		rownames(dt) = n
		dt = dt[which(rowSums(dt) > 0),]
		d[[j]]=dt
		#mar.default <- c(5,4,4,2) + 0.5
		#par(mar = mar.default + c(6, 4, 0, 0)) 
		#barplot(t(dt[1:250,]), las=2, legend = colnames(dt))
			
	}		
return(d)
}

d200=testConsis(graphs,iter=10,size=200)
d500=testConsis(graphs,iter=10,size=500)
d1000=testConsis(graphs,iter=10,size=1000)

> pdf("test.pdf", height=12, width=6)
 mar.default <- c(5,4,4,2) + 1
 par(mar = mar.default + c(6, 5, 1, 1)) 
 barplot(t(d200[[1]][50:100,]), las=2, horiz=TRUE, xlab =NULL, legend = colnames(d200[[1]]))
> dev.off()

mar.default <- c(5,4,4,2) + 0.5
par(mar = mar.default + c(6, 6, 0, 0)) 
barplot(t(dt[1:250,]), las=2,  horiz=TRUE, ylab = " ",legend = colnames(dt))	
barplot(t(d200[[1]])[,1:50], las=2, horiz=TRUE, xlab = " ",legend = colnames(d200[[1]]))
barplot(t(d200[[1]][50:100,]), las=2, horiz=TRUE, xlab = NULL, legend = colnames(d200[[1]]))
source("sp.path.compute.R")

testNodes=sample(vert, size=100, replace=FALSE)
testMap <-function(testNodes, graphs){
#subGraph(testNodes, g)
	names(graphs)=NULL	
	paths=list()
	paths2=list()
	for(k in 1:length(graphs)){
		testNodes=intersect(testNodes, nodes(graphs[[k]]))
		cat(testNodes)
		cat("\n")
		paths[[k]] = sp.path.compute(genelist=testNodes, G=graphs[[k]])
		paths2[[k]] = paths[[k]][sapply(paths[[k]], function(p) !is.na(p$length))]
	}
	#names(paths)=names(paths2)=names(graphs)
	Mod=list()
	Lit=list()
	for(j in 1:length(graphs)){
	g=subGraph(testNodes, graphs[[j]])
	Mod[[j]] = list()
	Lit[[j]] = list()
		for(m in 1:length(graphs)){
			#Mod[[j]][[m]]=validateMod(g, Allpaths=paths2[[m]])
			Lit[[j]][[m]]=validateLit(g, GG=graphs[[m]])	
		}	
	
	}
	return(Lit)
}

vert=nodes(graphs[[1]])
for(i in 2:length(graphs)){
vert=intersect(vert, nodes(graphs[[i]]))
}
vert=intersect(nodes(graphs[[1]]), intersect(nodes(graphs[[3]]),nodes(graphs[[2]])))
testNodes=sample(vert, size=50, replace=FALSE)
tm.50=testMap(testNodes, graphs)
sink("test.txt")
tm.50
sink()
testNodes=sample(vert, size=100, replace=FALSE)
tm.100=testMap(testNodes, graphs)
testNodes=sample(vert, size=200, replace=FALSE)
tm.200=testMap(testNodes, graphs)

### Run on KEGG
library(org.Hs.eg.db)
randWalk = function(g, n){
	mygraph = g
	u_mygraph = ugraph(mygraph)		
	start = sample(nodes(mygraph),1)
	no = start
	mynodes = start
	while(length(mynodes) < n){
		nei = adj(u_mygraph, no)
		no = sample(nei[[1]], 1)[[1]]           
		if(!(no %in% mynodes))
			mynodes = c(mynodes, no)
	}
	S = subGraph(mynodes, mygraph)		
	S = removeSelfLoops(S)	
	S
	geName = mget(nodes(S), org.Hs.egSYMBOL, ifnotfound=NA)
	return(list(Genes=geName, SubGraph=S))
}
simKegg(G.KEGG, n=10)

load("allKEGG.rda")
DIR="/home/praveen/Paurush/COSBI/interactomes/Results/KeggTest/"
e=c(20,50,100,150)
for(q in 1:length(graphs)){
for(j in 1:10){
#simKegg = function(G.KEGG, set){
	file=paste(DIR, "res", "size=", e[q],"_iter#",j, sep="_")	
	cand.net = randWalk(G.KEGG, n=e[q])
	print(cand.net$Genes)
	cand.net$SubGraph
	#plot(cand.net$SubGraph)
	cat("#######\n")	
	runGene = as.character(cand.net$Genes)
	runGene = runGene[which(!is.na(runGene))]
	print(runGene)
	myTM=testMap(runGene, graphs)
	save(myTM, cand.net, file=paste(file,".rda", sep=""))	
	#return(list(testRes=myTM, graph=cand.net))
}
}


DIR="/home/praveen/Paurush/COSBI/interactomes/Results/KeggTest/"
e=c(20,50,100)
for(q in 1:length(e)){
	for(j in 1:5){
		file=paste(DIR, "res", "size=", e,"_iter#",j, sep="_")
		myres=simKegg(G.KEGG, set=e[q])
		save(myres, file=paste(file,".rda", sep=""))
		pdf(paste(file, ".pdf", sep=""))		
		plot(myres$graph)
		dev.off()
	}
}


"6505"

[1] "HLA-DQB1"     "HLA-DRA"      "LOC100507709" "HLA-DQA2"     "HLA-DMA"     
 [6] "CD74"         "CTSL1"        "LGMN"         "HLA-DOA"      "HLA-DRB3"    
[11] "HLA-DQA1"     "LOC100507714" "HLA-DMB"      "HLA-DRB1"     "LOC101060835"
[16] "LOC100509457" "HLA-DRB5"     "HLA-DPA1"     "CD4"          "HLA-DRB4"    

 [1] "HLA-DQB1"     "HLA-DRA"      "LOC100507709" "HLA-DQA2"     "HLA-DMA"     
 [6] "CD74"         "CTSL1"        "LGMN"         "HLA-DOA"      "HLA-DRB3"    
[11] "HLA-DQA1"     "LOC100507714" "HLA-DMB"      "HLA-DRB1"     "LOC101060835"
[16] "LOC100509457" "HLA-DRB5"     "HLA-DPA1"     "CD4"          "HLA-DRB4"    

d=testMap(testNodes, graphs)
validateLit(g, GG=graphs[[m]])





detach("package:RSQLite", unload=TRUE)
detach("package:AnnotationDbi", unload=TRUE)
detach("package:org.Hs.eg.db", unload=TRUE)
detach("package:biomaRt", unload=TRUE)

library(RBGL)
library(org.Hs.eg.db)
load('/home/praveen/Paurush/COSBI/interactomes/sessionConsis.rda')
e=100
DIR="/home/praveen/Paurush/COSBI/interactomes/Results/KeggTest/"
for(j in 1:5){
#simKegg = function(G.KEGG, set){
	file=paste(DIR, "res", "size=", e,"_iter#",j, sep="_")	
	cand.net = randWalk(G.KEGG, n=e)
	print(cand.net$Genes)
	cand.net$SubGraph
	#plot(cand.net$SubGraph)
	cat("#######\n")	
	runGene = as.character(cand.net$Genes)
	runGene = runGene[which(!is.na(runGene))]
	print(runGene)
	myTM=testMap(runGene, graphs)
	save(myTM, cand.net, file=paste(file,".rda", sep=""))	
	#return(list(testRes=myTM, graph=cand.net))
}


















j=7
	file=paste(DIR, "res", "size=", e[q],"_iter#",j, sep="_")	
	cand.net = randWalk(G.KEGG, n=e[q])
	print(cand.net$Genes)
	cand.net$SubGraph
	#plot(cand.net$SubGraph)
	cat("#######\n")	
	runGene = as.character(cand.net$Genes)
	runGene = runGene[which(!is.na(runGene))]
	print(runGene)
	myTM=testMap(runGene, graphs)
	save(myTM, cand.net, file=paste(file,".rda", sep=""))	
j=8
	file=paste(DIR, "res", "size=", e[q],"_iter#",j, sep="_")	
	cand.net = randWalk(G.KEGG, n=e[q])
	print(cand.net$Genes)
	cand.net$SubGraph
	#plot(cand.net$SubGraph)
	cat("#######\n")	
	runGene = as.character(cand.net$Genes)
	runGene = runGene[which(!is.na(runGene))]
	print(runGene)
	myTM=testMap(runGene, graphs)
	save(myTM, cand.net, file=paste(file,".rda", sep=""))	

j=9
	file=paste(DIR, "res", "size=", e[q],"_iter#",j, sep="_")	
	cand.net = randWalk(G.KEGG, n=e[q])
	print(cand.net$Genes)
	cand.net$SubGraph
	#plot(cand.net$SubGraph)
	cat("#######\n")	
	runGene = as.character(cand.net$Genes)
	runGene = runGene[which(!is.na(runGene))]
	print(runGene)
	myTM=testMap(runGene, graphs)
	save(myTM, cand.net, file=paste(file,".rda", sep=""))	
j=10
	file=paste(DIR, "res", "size=", e[q],"_iter#",j, sep="_")	
	cand.net = randWalk(G.KEGG, n=e[q])
	print(cand.net$Genes)
	cand.net$SubGraph
	#plot(cand.net$SubGraph)
	cat("#######\n")	
	runGene = as.character(cand.net$Genes)
	runGene = runGene[which(!is.na(runGene))]
	print(runGene)
	myTM=testMap(runGene, graphs)
	save(myTM, cand.net, file=paste(file,".rda", sep=""))	



#########################################
### Script to compute consistency score # DONE WORKS *******************
#########################################
conScore=function(graphList, mode=c("direct", "path")){
	graphList[[1]]=NULL
	v=unique(unlist(lapply(allGraphs,nodes)))
	cmat=array(0, dim=c(length(v), length(v), length(graphList)))
	dimnames(cmat)[[1]]=dimnames(cmat)[[2]]=v
	dimnames(cmat)[[3]]=names(graphList)
	for(a in 1:nrow(cmat)){
		for(b in 1:ncol(cmat)){
			for(i in 1:length(graphList)){
				if(rownames(cmat)[a] %in% nodes(graphList[[i]])& colnames(cmat)[b] %in% nodes(graphList[[i]])){
					if(mode=="direct"){
						cmat[rownames(cmat)[a],colnames(cmat)[b],i]=as(graphList[[i]], "matrix")[rownames(cmat)[a],rownames(cmat)[b]]
					}					
					if(mode=="path"){
						cmat[rownames(cmat)[a],colnames(cmat)[b],i]=1/(sp.between(graphList[[i]],rownames(cmat)[a],rownames(cmat)[b], detail=TRUE)[[1]]$length)
					}				
				}
				else
					cmat[rownames(cmat)[a],colnames(cmat)[b],i]=NA
			}
		}	
	}
	score=mat_mean_sd(cmat)	
	return(score)
}
mat_mean_sd = function(mat){
mmat=matrix(NA, dim(mat)[1], dim(mat)[2])
smat=mmat
	for(i in 1:dim(mat)[1]){
		for(j in 1:dim(mat)[2]){
			mmat[i,j]=mean(mat[i,j,], na.rm=TRUE)
			smat[i,j]=sd(mat[i,j,], na.rm=TRUE)		
		}	
	}
mmat[is.nan(mmat)]<-NA
dimnames(mmat)=dimnames(smat)=list(dimnames(mat)[[1]], dimnames(mat)[[2]])
return(list(Mean=mmat, SD=smat))
}
scores=conScore(allGraphs, mode="path")

plotConsistencyMap=function(scores, file="scorePlot.pdf"){
	#par(mfrow=c(1,2))
	pdf(file)	
	scores$Mean[is.na(scores$Mean)]<--1
	scores$SD[is.na(scores$SD)]<--1
	myImagePlot(scores$Mean)
	myImagePlot(scores$SD)
	dev.off()
}
#########################################
### GO similarity Score computation	# DONE WORKS *******************
#########################################
library(GOSim)
library(org.Hs.eg.db)
goSimMatrix = function(names.nsgenes, ont){ 
	gg=mget(as.character(names.nsgenes), org.Hs.egALIAS2EG, ifnotfound=NA)
	gg=as.character(unlist(sapply(gg, "[[", 1))) # selects only first Entrez ID
	genesjanENT=as.character(unlist(sapply(gg, "[[", 1))) # selects only first Entrez ID
	setOntology(ont)	
	genes = gg
	SimMat=getGeneSim(genes,similarity="hausdorff", similarityTerm='Lin')#, similarityTerm = terms, verbose = FALSE)
	#colnames(SimMat)=genes
	#rownames(SimMat)=genes
	SimMat=	mergeMatrix(SimMat, gg)
	colnames(SimMat)=rownames(SimMat)=names.nsgenes
	return(SimMat)
	
}
v=unique(unlist(lapply(allGraphs,nodes)))
genesjanENT=as.character(unlist(sapply(gg, "[[", 1))) # selects only first Entrez ID
GOMF=goSimMatrix(names.nsgenes=v, ont="MF")
GOBP=goSimMatrix(names.nsgenes=v, ont="BP")
GOCC=goSimMatrix(names.nsgenes=v, ont="CC")

library(org.Hs.eg.db)
results=mget(gg, org.Hs.egGO)
tab=data.frame()
for(i in 1:length(results)){
	a=lapply(results[[i]], data.frame)
	a=do.call(rbind.data.frame, a)
	a=subset(a, Ontology == "CC", select = GOID)
	terms=c()
	for(j in 1:nrow(a)){
		t=mget(rownames(a)[j], GOTERM, ifnotfound=NA)[[1]]
		if(!is.na(t)){
			t=Term(t)
			terms=c(terms, t)
		}
	}
	tab=rbind(tab, data.frame(gene=v[i], CC=terms))
}
tab[,3]=1
library(reshape2)
f=acast(tab, gene~CC)
f[is.na(f)]<-0
myImagePlot(t(f))

####


mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("hgnc_symbol", "go_cellular_component_id",
           "go_cellular_component__dm_name_1006"), filters = "hgnc_symbol",
           values = v, mart = mart)


mergeMatrix=function(mat,genes){
	myMatrix = matrix(NA, ncol = length(genes), nrow = length(genes))
	colnames(myMatrix) = genes
	rownames(myMatrix) = genes
		for(i in 1:nrow(mat)){
			for(j in 1:ncol(mat)){
				v=mat[i,j]
				row=rownames(mat)[i]
				col=colnames(mat)[j]
				ind1=match(row, genes)
				ind2=match(col, genes)
				myMatrix[ind1, ind2]=v
			}
		}
	return(myMatrix)
 }

#########################################
### Compute transitively closed graph	# DONE WORKS *******************
#########################################
plot(transitive.closure(allGraphs[[3]]))
#########################################
### Script to compute graphs		# DONE WORKS *******************
#########################################
DIR="/home/praveen/Paurush/COSBI/interactomes/Results/KeggTest/20/"
allFiles=list.files(path = DIR, pattern=".rda")
dbases=c("HPRD","INTACT","BIOGRID")
for(i in 1:length(allFiles)){
	load(allFiles[i])
	db=length(myTM)
	for(j in 1:db){
		#refNet=data.frame(do.call('rbind', strsplit(as.character(myTM[[j]][[j]]$edge),'->',fixed=FALSE)))
		#refNet=		
		#colnames(refNet) = c("Source", "Target")
		allNets=lapply(myTM, data.frame)[[j]]
		x=colnames(allNets)
		x=grep("edge",x)
		x=x[-1]	
		allNets=allNets[,-x]
		colnames(allNets)=c("RefNet", dbases)
		allGraphs=apply(allNets[,1:ncol(allNets)],2, getGraph)	
	}
	score1=conScore(allGraph, mode="direct")
	score1=conScore(allGraph, mode="path")
	allSIF=lapply(allGraphs, Graph2DF)
	for(k in 1:length(allSIF)){
		allSIF[[k]]$edgeType=rep(colnames(allNets)[k], nrow(allSIF[[k]]))
	}
}
comNet=do.call("rbind", allSIF)
rownames(comNet)=NULL

#########################################
### Construct graphNEL from a dataframe # DONE WORKS *******************
#########################################

getGraph=function(dframe){
e=strsplit(as.character(dframe),'->',fixed=FALSE)
g1 = new("graphNEL", nodes=unique(unlist(e)), edgemode="directed")
for(i in 1:length(e)){
	for(j in 1:(length(e[[i]])-1)){
	#print(e[[i]][j])
	#print(e[[i]][j+1])
	g1 = addEdge(e[[i]][j],e[[i]][j+1], g1)
	}
}
return(g1)
}

#########################################
### Construct a SIF file with edge type # DONE WORKS *******************
#########################################

Graph2DF = function (g){
    edges.g <- edges(g)
    edge.names = as.character(unlist(sapply(names(edges.g), function (a) {
         bs = edges.g[[a]]; 
         if (length(bs) > 0) paste(a, edges.g[[a]], sep='~') 
         })))
    pairs = strsplit(edge.names, '~')
    a = sapply(pairs, "[", 1)
    b = sapply(pairs, "[", 2)

    if ('edgeType' %in% eda.names(g))
        edgeType = as.character(edgeData(g, from=a, to=b, attr='edgeType'))
    else
        edgeType = rep('unspecified', length (a))

   return(data.frame(source=a, target=b, edgeType=edgeType, stringsAsFactors=FALSE))
}

#########################################################
### 
#########################################################
########################################################################################################
write.table(comNet, quote=FALSE, row.names=FALSE, sep="\t","net.sif")

########################################################################################################
addNode(node, object, edges)
#########################################
### old junk code			# DONE WORKS *******************
#########################################

getGraph = function(dframe, dbase){
	g = new("graphNEL", nodes=unique(c(as.character(dframe[,1]), as.character(dframe[,2]))), edgemode="directed")
	g = addEdge(as.character(dframe[,1]), as.character(dframe[,2]), g)
	net=as.network(as(g, "matrix"), matrix.type = NULL, directed = TRUE,hyper = FALSE, loops = FALSE, multiple = FALSE, bipartite = FALSE, ignore.eval = TRUE, names.eval = NULL, na.rm = FALSE, edge.check = FALSE)
	plot.network(net, label = network.vertex.names(net), edge.col = db, vertex.col = db)
}



###########################################################################################################
###########################################################################################################

sink("tnet.txt")
for(k in 1:10){
	cat(paste(refNet[k,1], refNet[k,2], "\n" ,sep=" "))
}
sink()
sink("qnet.txt")
for(k in 1:10){
	cat(paste(refNet[k,1], refNet[k,2], "\n" ,sep=" "))
}
sink("sim.txt")
for(k in 1:length(nodes(g))){
for(l in 1:length(nodes(g))){
	if(nodes(g)[k]==nodes(g)[l])
		cat(paste(nodes(g)[k], 1, nodes(g)[l], "\n" ,sep=" "))
	else
		cat(paste(nodes(g)[k], 0, nodes(g)[l], "\n" ,sep=" "))	
}
}
###########################################################################################################
###########################################################################################################
       V <- LETTERS[1:4]
       edL1 <- vector("list", length=4)
       names(edL1) <- V
       for(i in 1:4)
          edL1[[i]] <- list(edges=c(2,1,4,3)[i], weights=sqrt(i))
       gR <- graphNEL(nodes=V, edgeL=edL1)
       gX <- addNode("X", gR)
     
     set.seed(123)
     g1 <- randomGraph(letters[1:10], 1:4, p=.3)
     g2 <- addNode("z", g1, edges=list(c("a", "h", "g")))


###########################################################################################################
###########################################################################################################
net.align("qnet.txt", "tnet.txt", "sim.txt", query.type=1,delta.d=1e-10, delta.c=0.0, delta.e=1, delta.s=1,   output="result.txt")
###########################################################################################################
###########################################################################################################

myImagePlot <- function(x, ...){
	min <- min(x)
	max <- max(x)
        yLabels <- rownames(x)
	xLabels <- colnames(x)
        title <-c()

  # check for additional function arguments

  	if( length(list(...)) ){
    		Lst <- list(...)
    		if( !is.null(Lst$zlim) ){
		       min <- Lst$zlim[1]
		       max <- Lst$zlim[2]
		}
		if( !is.null(Lst$yLabels) ){
	               yLabels <- c(Lst$yLabels)
		}
		if(!is.null(Lst$xLabels) ){
	        	xLabels <- c(Lst$xLabels)
        	}
		if( !is.null(Lst$title) ){
		title <- Lst$title
		}
	}

# check for null values
	if( is.null(xLabels) ){
		xLabels <- c(1:ncol(x))
	}
	if( is.null(yLabels) ){
		yLabels <- c(1:nrow(x))
	}
	layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
	ColorRamp <- rgb( seq(0.95,0.99,length=256),  # Red
                   seq(0.95,0.05,length=256),  # Green
                   seq(0.95,0.05,length=256))  # Blue
	ColorLevels <- seq(min, max, length=length(ColorRamp))
 # Reverse Y axis
	reverse <- nrow(x) : 1
	yLabels <- yLabels[reverse]
	x <- x[reverse,]
 # Data Map
	par(mar = c(3,5,2.5,2))
	image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",ylab="", axes=FALSE, zlim=c(min,max))
	if( !is.null(title) ){
		title(main=title)
	}
	axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7, las=3)
	axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,cex.axis=0.7)
 # Color Scale
	par(mar = c(3,2.5,2.5,2))
	image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp,xlab="",ylab="",xaxt="n")
	layout(1)

}


################################################################
## Edge wise clustering and module detection in Graphs
################################################################
# Compute 3D array for graph edge

fillArray=function(graphList){
	for(i in 1:length(graphList)){
	ed=edgesAll(graphList)
	myArray=array(NA, c(length(graphList), length(graphList), ed))
		for(j in 1:length(graphList)){
			if(i!=j){
				v=ComputeEdgeMap(graphList[[i]], graphList[[j]])
				rownames(v)=v[,1]
				
			}
			
		}
	}
}
compute3DArray=function(G1,G2){
	
	a=array()
	
}

# Compute Distance matrix for 3D array

# Perform Clustering


################################################################
## Read Signature File from Laura
################################################################
sigs=read.csv('/home/praveen/Paurush/COSBI/CorradoProjects/NDDsig.csv', head=TRUE, sep="\t" )
library(VennDiagram)

results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = gene, mart = mart)

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene=as.character(pathways[,"GeneB"])
gene=sapply(strsplit(gene, "[.]"), "[[", 2)
results <- getBM(attributes = c("ensembl_peptide_id", "hgnc_symbol"), filters = "hgnc_symbol", values = gene, mart = mart)

gene=as.character(sigs[,1])
library(org.Hs.eg.db)
gg=mget(as.character(gene),  org.Hs.egALIAS2EG, ifnotfound=NA)
gg=as.character(unlist(gg))
library(org.Hs.eg.db)
select(org.Hs.eg.db, cols=c("SYMBOL", "UNIPROT"), keys= as.character(gene), keytype="SYMBOL")


myVennSet=list(as.character(sigs[,1]), as.character(sigs[,2]), as.character(sigs[,3]), as.character(sigs[,6]), as.character(sigs[,5]))

names(myVennSet)=colnames(sigs)[c(1:3,6)]

venn.plot <- venn.diagram(x = myVennSet, filename = NULL, col = "transparent", fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),alpha = 0.50,c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),cex = 1.5,main="Sharing of genes signatures across NDDs",fontfamily = "serif",fontface = "bold",cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),cat.cex = 1.5,cat.pos = 0,cat.dist = 0.07,cat.fontfamily = "serif",rotation.degree = 270,margin = 0.2)

pdf("venn-genes.pdf")
grid.draw(venn.plot)
dev.off()


graphs=list(allgraphsEnt[[2]], allgraphsEnt[[3]], allgraphsEnt[[5]], allgraphsEnt[[7]], allgraphsEnt[[9]], allgraphsEnt[[11]])
 names(graphs)=c("Biogrid", "KEGG", "HPRD", "IntAct", "Reactome", "STRING")
# Simulation 2 entire database
testConsis=function(graphs, iter=1, size=c(200,500, 1000)){
	d=list()
	vert=nodes(graphs[[1]])	
	for(i in 2:length(graphs)){
		vert=intersect(vert, nodes(graphs[[i]]))
	}
	for(j in 1:iter){
		testNodes=sample(vert, size, replace=FALSE)
		dt=data.frame()
		for(k in 1:length(graphs)){
			dt = rbind(dt, c(as(subGraph(testNodes, graphs[[k]]), "matrix")))
		}
		dt=t(dt)
		colnames(dt) = names(graphs)
		n1=as(subGraph(testNodes, graphs[[k]]), "matrix")
		n=c()		
		for(a in 1:ncol(n1)){
			for(b in 1:nrow(n1)){
			n=c(n, paste(colnames(n1)[a], rownames(n1)[b], sep="->"))			
			}
		}
		rownames(dt) = n
		dt = dt[which(rowSums(dt) > 0),]
		d[[j]]=dt
		#mar.default <- c(5,4,4,2) + 0.5
		#par(mar = mar.default + c(6, 4, 0, 0)) 
		#barplot(t(dt[1:250,]), las=2, legend = colnames(dt))
			
	}		
return(d)
}


# connected components

connComp


