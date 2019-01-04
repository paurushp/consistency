##############################################################################
## PROJECT: CONSISTENCY IN PPI GRAPHS
## AUTHOR: PAURUSH PRAVEEN
## DATE: 13 MARCH 2014
## APPROACH: 
## DEPENDENCIES:
## CONTACT: praveen@cosbi.eu
##############################################################################
### GET DATABASES IN graphNEL FORMAT
library(RBGL)

# HPRD to graphNEL

pathways = read.csv("HPRD/HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt", head=FALSE, sep="\t")
pathways=pathways[,c(1, 7,4)]
colnames(pathways) = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)
# > G
# A graphNEL graph with directed edges
# Number of Nodes = 9617 
# Number of Edges = 39184 

G = removeSelfLoops(G)
hprd=G

# > hprd
# A graphNEL graph with directed edges
# Number of Nodes = 9617 
# Number of Edges = 37049 

# MINT to graphNEL

pathways = read.csv("MINT/2013-03-26-mint-human-binary.mitab26.txt",head=TRUE, sep="\t")

pathways=pathways[,c(1, 12,2)]
colnames(pathways) = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathways$GeneA), as.character(pathways$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways$GeneA), as.character(pathways$GeneB), G)

# > G
# A graphNEL graph with directed edges
# Number of Nodes = 9916 
# Number of Edges = 24939 

G = removeSelfLoops(G)
mint=G

# > mint
# A graphNEL graph with directed edges
# Number of Nodes = 9916 
# Number of Edges = 24458 



# BIOGRID to graphNEL


pathways = read.csv("MINT/2013-03-26-mint-human-binary.mitab26.txt",head=TRUE, sep="\t")
pathwaysSYM=pathways[,c(2, 19,3)]

colnames(pathwaysSYM) = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathwaysSYM$GeneA), as.character(pathwaysSYM$GeneB))), edgemode="directed")
G = addEdge(as.character(pathwaysSYM$GeneA), as.character(pathwaysSYM$GeneB), G)
# > G
# A graphNEL graph with directed edges
# Number of Nodes = 18110 
# Number of Edges = 146803 

G = removeSelfLoops(G)
biogridEnt=G

# > biogridSym
# A graphNEL graph with directed edges
# Number of Nodes = 18110 
# Number of Edges = 145036 


pathwaysSym=pathways[,c(8, 19,9)]
colnames(pathwaysSym) = c("GeneA", "INTERACTION_TYPE", "GeneB")
pathwaysSym=pathwaysSym[complete.cases(pathwaysSym),]
G = new("graphNEL", nodes=unique(c(as.character(pathwaysSym$GeneA), as.character(pathwaysSym$GeneB))), edgemode="directed")
G = addEdge(as.character(pathwaysSym$GeneA), as.character(pathwaysSym$GeneB), G)
# > G
# A graphNEL graph with directed edges
# Number of Nodes = 17762 
# Number of Edges = 146632 

G = removeSelfLoops(G)
biogridSym=G
# > biogridSym
# A graphNEL graph with directed edges
# Number of Nodes = 17762 
# Number of Edges = 144860 






# Intact to graphNEL
G


pathways = read.csv("Intact/psimitab/intact.txt", head=TRUE, sep="\t")
pathways=pathways[,c(1, 15,2)]
node=pathways2[,1]
x=grep(pattern="uniprotkb:", node)
y=as.character(node[x])
newintact=sapply(strsplit(y, ":"), "[",2)
z=c()
for(i in 1:length(newintact)){
if(length(intraIDMapper(newintact[i], species = "HOMSA", srcIDType="UNIPROT", destIDType = "EG"))!=0)
z[i]=intraIDMapper(newintact[i], species = "HOMSA", srcIDType="UNIPROT", destIDType = "EG")
else
z[i]="Unknown"
}

intact2=G
> intact2
A graphNEL graph with directed edges
Number of Nodes = 15651 
Number of Edges = 54127 

colnames(pathways2)[1:3] = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(pathways2$GeneA), as.character(pathways2$GeneB))), edgemode="directed")
G = addEdge(as.character(pathways2$GeneA), as.character(pathways2$GeneB), G)
# > G
# A graphNEL graph with directed edges
# Number of Nodes = 57665 
# Number of Edges = 200627 

G = removeSelfLoops(G)
intact=G
# > intact
# A graphNEL graph with directed edges
# Number of Nodes = 57665 
# Number of Edges = 196642 

### NAME NORMALIZATION OF NODES IN GRAPHS (1) ENTREZ



### NAME NORMALIZATION OF NODES IN GRAPHS (2) SYMBOL
# intact database [1] 50897 Uniprots

x=grep(pattern="uniprotkb:", nodes(intact2))
y=nodes(intact2)[x]
newintact=sapply(strsplit(y, ":"), "[",2)
nodes(intact)[x]=newintact
for(i in )
HumanEGs = intraIDMapper(newintact, species="HOMSA", srcIDType="UNIPROT",destIDType="EG")



library('org.Hs.eg')
annotation.col1 <- select(org.Hs.eg.db, keys=c('Q7TNF6','Q53XJ8','P05787','P0CG48','Q96CG1','D3DR86','Q96FS5'), cols=c('UNIPROT', 'SYMBOL', 'ENTREZID'), keytype="UNIPROT")

library('biomaRt')
ensembl <- useMart('ensembl', dataset="hsapiens_gene_ensembl")



annotation <- getBM(attributes=c("uniprot_swissprot_accession", "hgnc_symbol", "entrezgene"), filters="uniprot_swissprot_accession", values=newintact, mart=ensembl)


pathways2=subset(pathways, Taxid.interactor.A=="taxid:9606(human)|taxid:9606(Homo sapiens)"& Taxid.interactor.B=="taxid:9606(human)|taxid:9606(Homo sapiens)", select=c(1,15,2, 10,11))
head(pathways2)

geneB=as.character(pathways2[,3])
w=grep(pattern="uniprotkb:", geneB)
z=geneB[w]
newintact=sapply(strsplit(z, ":"), "[",2)
geneB=z


geneA=pathways2[,2]
a=data.frame()
for(i in 1:length(geneA)){

a=rbind(a, getBM(attributes=c("uniprot_swissprot_accession", "hgnc_symbol", "entrezgene"), filters="uniprot_swissprot_accession", values=geneA[i], mart=ensembl))
}

b=data.frame()
for(i in 1:length(geneB)){

b=rbind(b, getBM(attributes=c("uniprot_swissprot_accession", "hgnc_symbol", "entrezgene"), filters="uniprot_swissprot_accession", values=geneB[i], mart=ensembl))
}



symIntact=data.frame()
entIntact=data.frame()
for(i in 1:nrow(pathways3)){
	A=getBM(attributes=c("uniprot_swissprot_accession", "hgnc_symbol", "entrezgene"), filters="uniprot_swissprot_accession", values=sapply(strsplit(as.character(pathways3[i,1]),":"), "[",2), mart=ensembl)
	entA=A[,3]
	symA=A[,2]
	B=getBM(attributes=c("uniprot_swissprot_accession", "hgnc_symbol", "entrezgene"), filters="uniprot_swissprot_accession", values=sapply(strsplit(as.character(pathways3[i,3]),":"), "[",2), mart=ensembl)
	entB=B[,3]
	symB=B[,2]
	if(length(entB)!=0 & length(entA!=0)){
	for(j in 1:length(entA)){
		for(k in 1:length(entB)){
			entIntact=rbind(entIntact, data.frame(GeneA=entA[j], Interaction=pathways3[i,2], GeneB=entB[k]))
			
		}	
	}
	}
	if(length(symB)!=0 & length(symA!=0)){
for(j in 1:length(symA)){
		for(k in 1:length(symB)){
			symIntact=rbind(symIntact, data.frame(GeneA=symA[j], Interaction=pathways3[i,2], GeneB=symB[k]))
			
		}	
	}
	}
}


symIntact=rbind(symIntact, data.frame(GeneA=A[j], Interaction=pathways3[i,2], GeneB=B[k]))
sapply(strsplit(z, ":"), "[",2)








# Extract human interactions
pathways2=subset(pathways, Taxid.interactor.A=="taxid:9606(human)|taxid:9606(Homo sapiens)"& Taxid.interactor.B=="taxid:9606(human)|taxid:9606(Homo sapiens)", select=c(1,15,2, 10,11))
# Extract PPIs (look for protein IDs in inteacting molecule column)
s=grep(pattern="uniprotkb:", as.character(pathways2[,1]))
t=grep(pattern="uniprotkb:", as.character(pathways2[,3]))
length(intersect(s,t))
u=intersect(s,t)
pathways3=pathways2[u,]
dim(pathways3)
#

> s=grep(pattern="uniprotkb:", as.character(pathways2[,1]))
> t=grep(pattern="uniprotkb:", as.character(pathways2[,3]))
> length(intersect(s,t))
[1] 72198
> dim(pathways2)
[1] 74365     5
> u=intersect(s,t)
> pathways3=pathways2[u,]
> dim(pathways3)
[1] 72198     5
> head(pathways3)
## Convert to graph object
row_sub = apply(symIntact2, 1, function(row) all(row !="" ))
symIntact2=symIntact[complete.cases(symIntact),] # remove rows with NAs 
colnames(symIntact2)[1:3] = c("GeneA", "INTERACTION_TYPE", "GeneB")
G = new("graphNEL", nodes=unique(c(as.character(symIntact2$GeneA), as.character(symIntact2$GeneB))), edgemode="directed")
G = addEdge(as.character(symIntact2$GeneA), as.character(symIntact2$GeneB), G)
# > G
# A graphNEL graph with directed edges
# Number of Nodes = 57665 
# Number of Edges = 200627 
> G
A graphNEL graph with directed edges
Number of Nodes = 7950 
Number of Edges = 29573 

G = removeSelfLoops(G)
intactEnt=G
# A graphNEL graph with directed edges
# Number of Nodes = 7950 
# Number of Edges = 29161 

> G #(Symbol)
# A graphNEL graph with directed edges
# Number of Nodes = 7601 
# Number of Edges = 26944 

G = removeSelfLoops(G)
intactSym=G
# A graphNEL graph with directed edges
# Number of Nodes = 7601 
# Number of Edges = 26541 

save(list=ls(), file="allconsistency.rda")


# input matrices from multiple sources 

scoreCon=function(sources){
  sname = names(sources)
  nlist = length(sources)
  ew =list()
  for(i in 1:nlist){
    mat.here=1/sources[[i]]
    mat.here[is.na(mat.here)] <- 0
    graph.here=as(mat.here,"graphNEL")
    ew[[i]] = unlist(edgeWeights(graph.here))
    ew[[i]][ew[[i]]=="Inf"] <- 0
  }
  all.edges=c()
  for(i in 1:length(ew)){
  all.edges=union(names(ew[[i]]),all.edges)
  }
  compMat=matrix(0, nrow=nlist, ncol=length(all.edges))
  colnames(compMat)=all.edges
  rownames(compMat)=sname
  for(i in 1:nlist){
    compMat[i,names(ew[[i]])]=ew[[i]]
  }
  p=calcscore(compMat)
  return(p)
}


calcscore=function(edgewts){
s1=list(Con=mean(apply(edgewts, 2, mean)), Div=mean(apply(edgewts,2,sd)))
#   ns=nrow(edgewts)
#   ned=ncol(edgewts)
#   n1=ns-1
#   n2=ne-1
#   x=c()
#   for(i in 1:n1){
#     for(j in i+1:ns){
#       for(k in 1:n2){
# 	for(l in k+1:ne){
# 	   if(((compMat[i,k] > compMat[i,l]) & (compMat[j,k] >= compMat[j,l])) | ((compMat[i,k] < compMat[i,l]) & (compMat[j,k] <= compMat[j,l]))){
# 	    count=count+1
# 	}
#       }
#       
#     }
#   }
return(s1)
  }

cm=compMat

# Concordance Probability Method Score
computeConcProb=function(cm){
  p1=p2=0
  rr=nrow(cm)*(nrow(cm)-1)/2
  for(i in 1:(ncol(cm)-1)){
    t=i+1
    for(j in t:ncol(cm)){
      for(s in 1:(nrow(cm)-1)){
	count1=count2=0
	u=s+1
	for(v in u:nrow(cm)){
	  if((cm[s,i] > cm[s,j]) & (cm[v,i] >= cm[v,j])){
	    count1=count1+1
	  }
	  if((cm[s,i] < cm[s,j]) & (cm[v,i] <= cm[v,j])){
	    count2=count2+1
	  }
	}
	if(count1!=0){
	  temp=count1/rr
	 # cat(temp)
	  p1=p1+log(temp)
	}
	if(count2!=0){
	  temp=count2/rr
	  p2=p2+log(temp)
	}
      }
    #print(p1)
    #print(p2)
    }
  }
  #score=exp(p1/(ncol(cm)*(ncol(cm)-1)/2))+exp(p2/(ncol(cm)*(ncol(cm)-1)/2))
  score=(exp((p1)/(ncol(cm)*(ncol(cm)-1)/2))+exp((p2)/(ncol(cm)*(ncol(cm)-1)/2)))/2
  #print(score)
  return(score)
}

computeConcProb(cm)
computeConcProb(cm2)
computeConcProb(cm3)
computeConcProb(cm4)

computeConcProb(cm[c(1,4),])



computeVarScore=function(cm){
  s=mean(apply(cm,2,sd))
  s=1/s
  return(s)
}



computeCorScore=function(ss){
  co=0
  for(i in 1:(nrow(ss)-1)){
    h=i+1
    for(j in h:nrow(ss)){
      co=co+cor( as.numeric(ss[i,]),  as.numeric(ss[j,]))
    }
   }
  s=co/(nrow(ss)*(nrow(ss)-1)/2)
  return(s)
}


generateNoisy=function(g, noise, source){
  ss=data.frame(matrix(0,ncol=length(g), nrow=1))
  colnames(ss)=paste("E", 1:length(g), sep="_")
  #gg=c(as(g, "matrix"))
  for(so in 1:source){
#   print(so)
    gh=g
    n=round(noise*length(g)/100)
    n=sample(length(gh), n)
#   print(n)
    for(m in 1:length(n)){
      b=n[m]
#     print(gh[b])
      if(gh[b]==0){
	gh[b] = 1
      }
      else{
	gh[b] = 0
      }
    }
    temp = matrix(gh, nrow=1, ncol=ncol(ss))
    colnames(temp)=colnames(ss)
    ss=rbind(ss, temp)
  }
  ss=ss[-1,]
  return(ss)
}


##
plength=function(ss){
  ss2=ss
  for(s in 1:nrow(ss)){
    m = matrix(ss[s,],nrow=sqrt(ncol(ss)), ncol=sqrt(ncol(ss)), byrow=TRUE)
    pathL = matrix(0, nrow=nrow(m), ncol=ncol(m))
    G = as(m, "graphNEL")
    nodes(G)=genelist=paste("N", c(1:ncol(m)))
    for(i in 1:length(genelist)){
	for(j in 1:length(genelist)){
		if(i != j){
			path = sp.between(G, genelist[i], genelist[j])
			if(!is.na(path[[1]]$length)){
				pathL[i,j] = length(path[[1]][[3]][[1]])
			}
		}
	}
    }
  ss2[s,]=1/c(pathL)
  }
  ss2=as.matrix(ss2)
  ss2[is.infinite(ss2)]<-0
  rownames(ss2)=rownames(ss)
  return(ss2)
}
 
# ss1=plength(ss)
                   
################### KEGG based simulation
noisyNet=function(adj,noise){
  change=max(1,round(noise*ncol(adj)/100))
  sampRow=sample(ncol(adj), change)
  sampCol=sample(ncol(adj), change)
  #print(sampRow)
  #print(sampCol)
  adjN=adj
  operation=c(0,0,0)
  for(x in 1:length(sampRow)){
    for(y in 1:length(sampCol)){
      l=sampRow[x]
      m=sampCol[y]
      #print(l)
      #print(m)
      if(adj[l,m] == 0)
	operation[1] = 1
      if(adj[l,m] == 1 )
        operation[2] = 1
      if(adj[l,m]!=adj[m,l])
        operation[3] = 1
      op = sample(1:3, 1)        # Random selection of operation
      if(operation[op]==1 && op==1)
        adjN[l,m]= 1 # add edge
      if(operation[op]==1 && op==2)
        adjN[l,m]= 0 # remove edge
      if(operation[op]==1 && op==3){
        tempLV= adjN[l,m]
        adjN[l,m] = adjN[m,l] # swap direction
        adjN[m,l] = tempLV
      }

    }
  }
  return(adjN)
  }


generateSources=function(adj, noise, source){
  ss=data.frame(matrix(0,ncol=(ncol(adj)*nrow(adj)), nrow=1))
  #print(ss)
  colnames(ss)=paste("E", 1:length(adj), sep="_")
	for(k in 1:source){
		#alpha = rnorm(1, 2, noise/100)
		#beta = rnorm(1, 2, noise/100)
		if(noise==0){
		  adj=as(adj, "graphNEL")
		  temp=computeSP(nodes(adj), adj) # Comment for noisy uncomment above
		}
		adjNew=noisyNet(adj, noise)
		adjNew=as(adjNew, "graphNEL")
		temp=computeSP(nodes(adjNew), adjNew)
		temp=c(temp)
		temp[is.na(temp)]<-0
		#print(temp)
		ss=rbind(ss,c(temp))
	}
	ss=ss[-1,]
	return(ss)
}


#### Based on beta distribution
generateSource=function(adj, alpha, beta){
	size=ncol(adj)
	s=matrix(0, size, size)
	for(d in 1:size){
		for(e in 1:size){
			if(adj[d,e]==1)
				s[d,e]=rbeta(1,alpha,1)
			if(adj[d,e]==0)
				s[d,e]=rbeta(1,1,beta)
		}
	}
	#print(s)
#	s=c(s)
	return(s)
}

# Compute shortest path length

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



#ss=generateSources(ag, 5, source=3)

load('/home/praveen/Paurush/Prior/StrPrior/Data/KEGGgraphs.rda')
#########################################
### Sample graph from KEGG ##############
#########################################

sampleKEGGPathway = function(graphs, nacc, n){
	mygraph = sample(graphs, 1)
	mynacc = nacc[[names(mygraph)]]
	mygraph = mygraph[[1]]
	u_mygraph = ugraph(mygraph)		
	start = names(sample(mynacc[mynacc >= n],1))
	mynodes = start
	no = start
	while(length(mynodes) < n){
		nei = adj(u_mygraph, no)
		no = sample(nei[[1]], 1)[[1]]           
		if(!(no %in% mynodes))
			mynodes = c(mynodes, no)
	}
	S = subGraph(mynodes, mygraph)		
	S = removeSelfLoops(S)	
	S
}


#############
set.seed(123)
V = letters[1:10]
M = 1:4
g1 = randomGraph(V, M, 0.2)



# Simulation

# Simulate 
set.seed(123)
simulateTest=function(n=c(10,20,30,40,50), noise=c(5,10,15,20,25,50), sources=c(3,4,5,6), graph.type="random",sgen="edges"){
  resmat=matrix(0, nrow=1, ncol=7)
  resmat=data.frame(resmat)
  ag=list()
  for(sims in 1:simulations){
    if(graph.type=="random"){
	M=1:4
	p=0.2
	ag[sims]=c(as(randomGraph(V, M, 0.2),"matrix"))
    }
    if(graph.type=="power.law"){
      ag[sims]=barabasi.game(n[i])
      ag[sims]=c(as(igraph.to.graphNEL(ag[sims]),"matrix"))
    }
    if(graph.type=="KEGG"){
      ag[sims]=sampleKEGGPathway(graphs,nacc,n)
      ag[sims]=as(ag[sims], "matrix")
    }
  }
      save(ag, file="simnets.RData")
      for()
colnames(resmat)=c("Sources","Size", "Noise", "CS", "CP", "VS", "KC")
src=paste("S", c(1:6), sep="_")
  for(s in 1:length(sources)){
    source=sources[s]
   # cat("Sources \t")
   # print(source)
    #cat("\n")
    for(i in 1:length(n)){
    #  cat("Nodes \t")
    #  print(n[i])
     # cat("\n")
      V=paste("V", c(1:n[i]))
      
      for(j in 1:length(noise)){
	#cat("Noise \t")
	#print(noise[j])
	#cat("\n")
	  if(sgen=="edges"){
	    ss=generateNoisy(ag, noise[j], source)
	    rownames(ss)=src[1:nrow(ss)]
	  }
	  if(sgen=="plength"){
	    ss=generateSources(ag, noise[j], source)
	    rownames(ss)=src[1:nrow(ss)]
	  }
	
	cp=computeConcProb(ss)
	vs=computeVarScore(ss)
	cs=computeCorScore(ss) #0.759788518571958
	kc=computeKendalConc(ss)
# 	print(source)
# 	print(n[i])
# 	print(noise[j])
# 	cat("CS score \n")
# 	print(cs)
# 	cat("CP score \n")
# 	print(cp)
# 	cat("VS score \n")
# 	print(vs)
	#resmat=rbind(resmat, data.frame(Sources=source, Size=n[i] ,Noise=noise[j], CS=cs, CP=cp, VS=vs))
	resmat=rbind(resmat, data.frame(Sources=source, Size=n[i], Noise=noise[j], CS=cs, CP=cp, VS=vs, KC=kc))
	cat(paste("Source =", source, "| Nodes = ", n[i], "| Noise = ", noise[j], "| CS = ", cs, "| CP = ", cp, "| VS = ", vs,"| KC = ", kc, "\n", sep=" "))
      }
    }
  }
  #resmat=resmat[-1,]
  return(resmat)
 }

 
 ##### Proof of principle Simulations
 
 # Simulate for random graph direct edges
simres.rand.ed=list() 
for(a in 1:5){
cat(paste("Running simulation >> ", a, "\n", sep=""))
simres.rand.ed[[a]]=simulateTest(n=c(10,20,30,40,50), noise=c(5,10,15,20,25,50), sources=c(3,4,5,6), graph.type="random", sgen="edges")
}
# Simulate for scale free graph direct edges
simres.sf.edges=list() 
for(a in 1:5){
  cat(paste("Running simulation >> ", a, "\n", sep=""))
  simres.sf.edges[[a]]=simulateTest(n=c(10,20,30,40), noise=c(10,15,20,25,50), sources=c(3,4,5,6), graph.type="power.law", sgen="edges")
}

# Simulate for KEGG graph path length

simres.ppi.plength3=list() 
for(a in 1:5){
  cat(paste("Running simulation >> ", a, "\n", sep=""))
  simres.ppi.plength3[[a]]=simulateTest(n=c(10,20,30,40, 50), noise=c(10,15,20,25,50), sources=c(3,4,5,6), graph.type="power.law", sgen="plength")
}

simres.kegg.plength3=list() 
for(a in 1:5){
  cat(paste("Running simulation >> ", a, "\n", sep=""))
  simres.ppi.plength3=simulate_score_test(n=c(10,20,30,40,50,60), noise=c(5,10,15,20,25,30), num_sources=c(3,4,5,6), evidence="plength", net_type="power.law", trials=5)
}

simres.kegg.plength3=list() 
for(a in 1:5){
  cat(paste("Running simulation >> ", a, "\n", sep=""))
simulate_score_test(n=c(10,20,30,40, 50), noise=c(10,15,20,25,50), num_sources=c(3,4,5,6), net_type="power.law", evidence="plength")
}

simres.kegg.plength.no.noise=list() 
for(a in 1:5){
  cat(paste("Running simulation >> ", a, "\n", sep=""))
  simres.kegg.plength.no.noise[[a]]=simulateTest(n=c(10,20,30,40.50), noise=0, sources=c(3,4,5,6), graph.type="KEGG", sgen="plength")
}
 ##### Running for House keeping and non house keepig stuff
 source('/home/praveen/Paurush/COSBI/CorradoProjects/Functions.R')

 ### Cross talk detector
 install.packages('/home/praveen/Downloads/RSQLite-0.11.4.tar.gz', repos=NULL, type="source") #RSQLite version 1.0.0, which is not backwards-compatible so switch to old one
 
 library(org.Hs.eg.db)
 library(ReactomePA)
 library(graphite)
 g = convertIdentifiers(reactome$"mTOR signalling", type="symbol")
 g = pathwayGraph(pathway=g, edge.types=NULL)
 
 pathwayGraph(pathway, edge.types=NULL) # Biocarta pathway
 viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList) # reactome pathway
 
 ipPath=read.csv('/home/praveen/Paurush/COSBI/interactomes/pathways.csv', head=TRUE, sep="\t")
 seedPathways=unique(ipPath[,3])
 
 pathDB=list(biocarta,spike,nci,kegg, reactome)
 
 seedwd=c("Glycolysis", "Pentose", "Glucose", "TCA", "Oxidative", "Aspartate", "Leucine", "Cholesterol", "Threonine", "Proline", "Biosynthesis", "Arginine", "ketone", "Polyamine", "Creatine", "Purine", "Pyrimidine")
 for(i in 1:5){
  for(j in 1:length(seedwd)){
  x= grep(seedwd[j], names(pathDB[[i]]))
  cat(i)
  cat("\t")
  cat(seedwd[j])
    cat("\t")
  print(x)
  #cat(paste(names(pathDB[i]), seedwd[j], x,"\n" ,sep=" "))
  }
}


gluc.nci=
gluc.kegg=
gluc.reactome=

gluc.nci=convertIdentifiers(nci$"Glucose metabolism", type="symbol")
gluc.nci2=pathwayGraph(pathway=gluc.nci2, edge.types=NULL)
gluc.nci

g1=pathwayGraph(pathway=g1, edge.types=NULL)
g2=convertIdentifiers(kegg[[152]], type="symbol")
g2=pathwayGraph(pathway=g2, edge.types=NULL)
g3=convertIdentifiers(kegg[[153]], type="symbol")
g3=pathwayGraph(pathway=g3, edge.types=NULL)
g4=convertIdentifiers(kegg[[46]], type="symbol")
g4=pathwayGraph(pathway=g4, edge.types=NULL)
g1=convertIdentifiers(reactome[[765]], type="symbol")
g1=pathwayGraph(pathway=g1, edge.types=NULL)
g2=convertIdentifiers(reactome[[441]], type="symbol")
g2=pathwayGraph(pathway=g2, edge.types=NULL)
g3=convertIdentifiers(reactome[[442]], type="symbol")
g3=pathwayGraph(pathway=g3, edge.types=NULL)
g4=convertIdentifiers(reactome[[1106]], type="symbol")
g4=pathwayGraph(pathway=g4, edge.types=NULL)
g5=convertIdentifiers(reactome[[219]], type="symbol")
g5=pathwayGraph(pathway=g5, edge.types=NULL)
g6=convertIdentifiers(reactome[[860]], type="symbol")
g6=pathwayGraph(pathway=g6, edge.types=NULL)
g7=convertIdentifiers(reactome[[1124]], type="symbol")
g7=pathwayGraph(pathway=g7, edge.types=NULL)


compSPmat=function(dbsize, dbs){
  mats1=list()
  for(i in 1:dbsize){
    mats1[[i]] = computeSP(v=nodes(dbs[[i]]),g=dbs[[i]])
  }
  return(mats1)
}

joiner=function(db, index){
mylist=db[index]
G=convertIdentifiers(mylist[[1]], type="symbol")
G=pathwayGraph(pathway=G, edge.types=NULL)
if(length(index)>1)
for(i in 1:length(mylist)){
  G=join(G,pathwayGraph(pathway=convertIdentifiers(mylist[[i]], type="symbol"),edge.types=NULL))
}
return(G)
}

g3=joiner(reactome, reactomenfkb)

msp_nfkb=compSPmat(3,db_nfkb)
msp_insulin=compSPmat(3,db_insulin)
msp_toll=compSPmat(3,db_toll)

> keggnfkb=134
> keggtol=201
> kegginsulin=108
> nciInsulin=c(375:378, 574:577, 664)
> ncinfkb=c(49,83,118,343,729,330, 461,462)
> ncitoll=c(731:736)
> reactomeinsulin=c(518,527:530,896:900,1017,1031,1106)
> reactometoll=c(1127:1146)
> reactomenfkb=c(43,699,700,756,1148,1153)


for(i in 1:dbsize){
  interface = intersect(genes, nodes(dbs[[i]]))
  neib = intersect(intersect(nodes(dbs[[1]]), nodes(dbs[[2]])), nodes(dbs[[3]]))
  ss1 = ss2 = data.frame(matrix(0,ncol=length(neib), nrow=1))
  for(j in 1:length(interface)){
    temp1 = c(1/as.matrix(mats1[[i]][interface[j],neib]))
    temp2 = c(1/as.matrix(mats1[[i]][neib,interface[j]]))
    temp1[is.na(temp1)]<-0
    temp2[is.na(temp2)]<-0
    ss1 = rbind(ss1,temp1)
    ss2 = rbind(ss2,temp2)
}
 }
 
crossTalker2=function(dbmat, dblist, markers){
    neib=colnames(dbmat[[1]])
    for(i in 2:length(dbmat)){
      neib=intersect(neib, colnames(dbmat[[i]]))
    }
    x=1
    gn=nei=trlist1=trlist2=list()
    tab=data.frame(matrix(0, nrow=1,ncol=8))
    colnames(tab)=c("db", "marker", "cs1", "cs2", "pv","rank1", "rank2", "rank3")
    for(j in 1:length(dbmat)){
      hm=names(dblist)[j]
      interface=intersect(markers, colnames(dbmat[[j]]))
      ss1 = ss2 = data.frame(matrix(0,ncol=length(neib), nrow=1))
      
      for(k in 1:length(interface)){
	r=c()
	for(l in 1:length(dbmat)){
	print(paste("####", hm, sep=" "))
	print(interface[k])
	print(names(dblist)[l])
	cat("\n")
	if(interface[k] %in% nodes(dblist[[l]])){
	  r=c(r,betweenness(graph=igraph.from.graphNEL(dblist[[l]]), v = interface[k], directed = TRUE))#r=c(r,page.rank(graph=igraph.from.graphNEL(dblist[[l]]), algo ="prpack",vids = interface[k], directed = TRUE))$vector
	}else{
	  r=c(r,0)
	} 
	if(interface[k] %in% colnames(dbmat[[l]])){
	cat("YES")
	gn[[x]]=interface[k]
	nei[[x]]=neib
	tr1=as.matrix(dbmat[[l]][interface[k],neib])
	tr2=as.matrix(dbmat[[l]][neib,interface[k]])
	temp1 = c(1/tr1)
	temp2 = c(1/tr2)
	}else{
	temp1 =rep(0, length(interface))
	temp1 =rep(0, length(interface))
	}
	temp1[is.na(temp1)]<-0
	temp2[is.na(temp2)]<-0
	ss1 = rbind(ss1,temp1)
	ss2 = rbind(ss2,temp2)
	}
	ss1=ss1[-1,]
	ss2=ss2[-1,]
	trlist1[[x]]=ss1
	trlist2[[x]]=ss2
	x=x+1
	cs1=computeConcProb(ss1)
	cs2=computeConcProb(ss2)
	tab=rbind(tab, data.frame(db=names(dblist)[j], marker=interface[k], cs1=cs1$Score, cs2=cs2$Score, pv=cs1$p.value,rank1=r[1], rank2=r[2], rank3=r[3]))
      }
      
    }
    tab=tab[-1,]
    res=list(table=tab, interface=neib, nb=nei, ge=gn, ilist1=trlist1,ilist2=trlist2)
    return(res)
}
 
allscores_1=list()
allscores_1[[1]]=crossTalker(dbmat=msp_gluc, dblist=db_gluc, markers=genes)
allscores[[2]]=crossTalker(dbmat=msp_insulin, dblist=db_insulin, markers=genes)
allscores[[3]]=crossTalker(dbmat=msp_toll, dblist=db_toll, markers=genes) 
allscores[[4]]=crossTalker(dbmat=msp_nfkb, dblist=db_nfkb, markers=genes) 
names(allscores)=c("Glucose","Insulin","TOLL","NFKB")

restab=data.frame(matrix(0, nrow=1,ncol=13))
    colnames(restab)=c("Pathway", "Marker", "Frequency","CS1","CS2","RNK1", "RNK2", "RNK3", "SDCS1", "SDCS2", "SDRNK1","SDRNK2","SDRNK3")
    for(i in 1:length(allscores)){
      mytab.here=allscores[[i]]$table
      markers=unique(mytab.here$marker)
      print(markers)
      for(j in 1:length(markers)){
	  tab.here=as.data.frame(subset(allscores[[i]]$table, marker==markers[j]))
	  #print(tab.here)
	  #print(markers[j])
	  mn=apply(tab.here[,c(3:7)],2,mean)
	  sdv=apply(tab.here[,c(3:7)],2,sd) 
	  mn=t(data.frame(mn))
	  sdv=t(data.frame(sdv))
	  colnames(mn)=c("CS1","CS2","RNK1", "RNK2", "RNK3")
	  colnames(sdv)=c("SDCS1", "SDCS2", "SDRNK1","SDRNK2","SDRNK3")
	  restab=rbind(restab, data.frame(Pathway=names(allscores)[i],Marker=markers[j], Frequency=nrow(tab.here),cbind(mn,sdv)))
      }
    }
    
 
val=list(x=restab1$Frequency,y=restab1$CS2,z=restab1$CS1)

open3d()
surface3d(x, y, z, color = col, back = "lines")

ggplot(restab1, aes(x=Frequency, y=CS1))+geom_line()+stat_smooth(method="lm")

library(scatterplot3d)
# attach(mtcars)
# scatterplot3d(wt,disp,mpg, pch=16, highlight.3d=TRUE, type="h", main="3D Scatterplot") 
  
col=rep("red",37)
col[which(as.numeric(factor(color1))==1)]="black"
s3d=scatterplot3d(restab1$Frequency,restab1$CS1,restab1$CS2, pch=16, xlab="Frequency of occurance" , ylab="Score for outgoing edge",zlab="Score for incoming edge",highlight.3d=FALSE, color=as.numeric(factor(restab1$RNK1)), type="h", main="Knowledge-consistency space")

s3d.coords <- s3d$xyz.convert(restab1$Frequency,restab1$CS1,restab1$CS2)
#al.char <- toupper(substr(as.character(Allocation), 1, 1))
text(s3d.coords$x, s3d.coords$y, labels = restab1$Marker, pos = 2, offset = 0.50, color=as.numeric(factor(restab1$RNK1)))

layout(cbind(1:2, 1:2), heights = c(7, 1))
s3d=scatterplot3d(restab1$Frequency,restab1$CS1,restab1$CS2, pch=16, xlab="Frequency of occurance" , ylab="Score for outgoing edge",zlab="Score for incoming edge",highlight.3d=FALSE, color=as.numeric(factor(restab1$RNK1)), type="h", main="Knowledge-consistency space")
par(mar=c(5, 3, 0, 1))
plot(seq(min(as.numeric(factor(restab1$RNK1))), max(as.numeric(factor(restab1$RNK1))), length = length(unique(as.numeric(factor(restab1$RNK1))))), rep(0, length(unique(as.numeric(factor(restab1$RNK1))))), pch = 16,axes = FALSE, xlab = "color code of centrality \"Pagerank method\"",ylab = "", col = sort(unique(as.numeric(factor(restab1$RNK1))), decreasing=FALSE))
axis(1, at = 1:length(unique(restab1$RNK1)), labels = sort(unique(restab1$RNK1)))

unique(restab2$RNK1)
  
  
plotTimeNetwork = function(graph, layout=c())
net=graph.adjacency(result.meta$network,mode="directed",weighted=T,diag=FALSE)
E(GlucI)$curved=0.25
V(GlucI)$color= "yellow"
V(GlucI)$color[which(genes%in%V(GlucI)$name)]="red"

V(GlucI)$size=15
plot.igraph(GlucI,vertex.label=V(GlucI)$name,edge.label=NULL,layout=layout.reingold.tilford(g, circular=T))  # plot the graphite

E(net)$weight2= gsub(1,"green",E(net)$weight)  # FAST
E(net)$weight2= gsub(7,"red",E(net)$weight2)   # SLOW
E(net)$weight2= gsub(6,"red",E(net)$weight2)   # SLOW
E(net)$weight2= gsub(5,"red",E(net)$weight2)   # SLOW
E(net)$weight2= gsub(4,"cyan",E(net)$weight2)  # INTERMEDIATE
E(net)$weight2= gsub(3,"cyan",E(net)$weight2)  # INTERMEDIATE
E(net)$weight2= gsub(2,"green",E(net)$weight2) # FAST
E(net)$color=E(net)$weight2

title(main=paste("Network For Murine Stem Cell Development"), sub= "\t\tFast Signalling\t\t", col.sub='green')
title(sub= "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
title(sub= "\t\tIntermediate Signalling\t\t", col.sub='cyan')
title(sub= "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
title(sub= "\t\t Slow Signalling\t", col.sub='red')


Layouts in iGraph
layout.spring
layout.kamada.kawai
layout.lgl

 #Plot Functions

 restab1=subset(restab, Pathway=="Gluc")
  color1=factor(restab1$RNK1)
 restab1=data.frame(restab1, color1=color1, tsize=restab1$Frequency*4)
 bbubc1 <- ggplot(restab1, aes(x = -log(CS1), y = -log(CS2), size = RNK1, label = Marker, colour=RNK1)) + #colour = factor(SDCS1))
    geom_point(aes(shape = factor(Frequency))) + 
    scale_size(range = c(5, sqrt(max(restab1$Frequency)/min(restab1$Frequency)*5^2)), name = "Frequency in \n Pathway databases") +
    geom_text(size = restab1$tsize, colour = "black", vjust = -5, position = position_jitter(width=-1, height=2)) + 
    #scale_colour_manual(values = levels(factor(restab$SDCS1)), name = "Standard deviation of score", label = levels(factor(restab$SDCS1))) + 
    xlab("Score (IN)") +
    ylab("Score (OUT)") + 
    #opts(axis.text.x=theme_text(angle=0, hjust=0, size = 16)) + 
    #opts(axis.text.y=theme_text(angle=0, hjust=0, size = 16)) +
    #opts(axis.title.x=theme_text(size = 16)) + 
    #opts(legend.title = theme_text(size = 16)) + 
    #opts(axis.title.y=theme_text(size = 16, angle = 90)) + 
    geom_vline(colour = I("grey")) + 
    geom_hline(colour = I("grey")) + 
    ylim(c(0,30))+
    xlim(c(0,30))
    bbubc1 
    
 
generateSources=function(adj, noise, source){
  ss=data.frame(matrix(0,ncol=(ncol(adj)*nrow(adj)), nrow=1))
  #print(ss)
  colnames(ss)=paste("E", 1:length(adj), sep="_")
	for(k in 1:source){
		#alpha = rnorm(1, 2, noise/100)
		#beta = rnorm(1, 2, noise/100)
		if(noise==0){
		  adj=as(adj, "graphNEL")
		  temp=computeSP(nodes(adj), adj) # Comment for noisy uncomment above
		}
		adjNew=noisyNet(adj, noise)
		adjNew=as(adjNew, "graphNEL")
		temp=computeSP(nodes(adjNew), adjNew)
		temp=c(temp)
		temp[is.na(temp)]<-0
		#print(temp)
		ss=rbind(ss,c(temp))
	}
	ss=ss[-1,]
	return(ss)
}












# Python code for 3d plot



x_surf=[[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3]]
y_surf=[[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10]]
z_surf=[[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10]]


from mpl_toolkits.mplot3d import *
import matplotlib.pyplot as plt
import numpy as np
from random import random, seed
from matplotlib import cm


z=[1, 3, 2, 2, 2, 2, 2, 1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
x=[8.572998e-01, 4.068936e-01, 4.368407e-02, 1.053834e-02, 1.281430e-02,1.065945e-02, 3.506698e-02, 1.797963e-02, 8.581283e-02, 1.428335e-02,1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 5.267207e-02, 6.484517e-02, 7.281187e-02, 7.840046e-02, 8.244115e-02, 8.540193e-02, 8.757703e-02, 8.916416e-02, 9.030222e-02, 8.224146e-05, 6.977281e-06, 1.195190e-11,2.746540e-11, 4.433214e-11]
y=[5.731617e-01, 3.303881e-01, 1.140550e-02, 7.120484e-03, 1.062121e-02, 6.298342e-03, 6.026997e-03, 4.153908e-04, 8.468981e-02, 1.406555e-02,6.125711e-01, 7.399015e-01, 1.000000e+00, 1.000000e+00, 1.000000e+00,1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,1.000000e+00, 1.000000e+00, 1.000000e+00, 5.137235e-03, 9.218224e-03,1.200756e-02, 1.407228e-02, 1.562898e-02, 1.680907e-02, 1.770226e-02,1.837320e-02, 1.886986e-02, 3.917972e-08, 1.403402e-11, 1.182024e-18,1.181652e-17, 4.204634e-17]

x_surf=[1,2,3]
y_surf=[0.5, 0.25, 0.08]
x_surf, y_surf = np.meshgrid(x_surf, y_surf)
z_surf=[[0.5, 0.27, 0.10],[0.5, 0.27, 0.10],[0.5, 0.27, 0.10]]

x_surf=[1,1,1,1,2,2,2,3,3,3]
y_surf=[1,0.2,0.6,0.5,0.6,0.4,0.5,0.08,0.12,0.10]
x_surf, y_surf = np.meshgrid(x_surf, y_surf)
z_surf=np_array([[0.2,1,0.6,0.5,0.4,0.6,0.5,0.12,0.08,0.10]])
z_surf=np.repeat(z_surf, 10, axis=0)
z_surf = np.sqrt(y_surf.T)
fig = plt.figure()
ax = fig.gca(projection='3d')               # to work in 3d
plt.hold(True)
ax.plot_surface(x_surf, y_surf, z_surf, cmap=cm.hot)

ax.scatter(x, y, z);                        # plot a 3d scatter plot

ax.set_xlabel('x label')
ax.set_ylabel('y label')
ax.set_zlabel('z label')

plt.show()


 ### Read BioPAX files
 



# Result analysis and plots
library(reshape)
library(ggplot2)
library(plotrix)

plot.analysis=function(res, filename){
  for(i in 1:length(res)){
    res[[i]]=res[[i]][-1,]
  }
  filename=paste(filename, ".pdf", sep="")
  res.here=array(0, dim=c(nrow(res[[1]]), ncol(res[[1]]), length(res)))
  dimnames(res.here)[[2]]=c("Sources", "Size", "Noise", "Kendall", "CP", "FCCS", "Lin")
  dimnames(res.here)[[3]]=paste("trials", 1:length(res), sep="")
  for(i in 1:length(res)){
    res.here[,,i]=as.matrix(res[[i]])
  }
 ### Plots
# Sources Vs Scores
#   bplot_SS=data.frame(matrix(0, ncol=4+length(res), nrow=1))
#   scoring=data.frame(Scoring=melt(res[[1]][,4:6])[,1])
#   df=list()
#   for(i in 1:length(res)){
#   df[[i]]=melt(res[[i]][,4:6])[,2]
#   }
#   df=as.data.frame(df)
#   colnames(df)=paste("Trial", 1:length(res))
#   df=data.frame(Source=rep(res[[1]][,1], 3), Size=rep(res[[1]][,2], 3), Noise=rep(res[[1]][,3], 3), Scoring= scoring, df)
  allMeans=t(apply(res.here[,c(4:7),], 1, rowMeans))
  allSD=apply(res.here[,4:7,], c(1,2), sd)
  colnames(allMeans)=colnames(rr[[i]]) 
  colnames(allSD)=paste("SD",colnames(rr[[i]])[4:7],sep=".")
  plotRes=data.frame(allMeans, allSD)
  V=unique(plotRes$Size)
  method=c("Kendall", "CP", "FCCS", "Lin")
    pdf(filename)
#     col=c("red", "green", "blue")
  for(v in 1:length(V)){
    data.here=subset(plotRes, Size==V[v])
    param=data.here[,1:4]
    #dataSel=data.frame(param, Method=rep(method[2], nrow(param)), Score=data.here[,method[2]], SD=data.here[,paste("SD", method[1], sep=".")])

    for(m in 1:length(method)){ #length(method)){
#       dataSel=rbind(dataSel, data.frame(param, Method=rep(method[m], nrow(param)), Score=data.here[,method[m]], SD=data.here[,paste("SD", method[m], sep=".")]))
	title.here=paste("Consistency score", method[m], "#nodes =", V[v], sep=" ")
        dataSel= data.frame(param, Method=rep(method[m], nrow(param)), Score=data.here[,method[m]], SD=data.here[,paste("SD", method[m], sep=".")])
        p=ggplot(dataSel, aes(x=Noise, y=Score))+geom_errorbar(aes(ymin=Score-SD, ymax=Score+SD), width=.1)+geom_point()+stat_smooth(method="lm")+facet_grid(~ Sources)
	p=p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_rect(aes(fill = Sources), xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.1)+ ggtitle(title.here)
	print(p)
    }
   } 
    dev.off()
  }
   plot.analysis(simres.sf.ed, file="simres.plength2")


   
   plotRes2=function(obj, filename){
    res.here=array(NaN, dim=c(nrow(obj)/max(obj[,1]), ncol(obj)-1, max(obj[,1])))
    #dimnames(res.here)[[2]]=c("Sources", "Size", "Noise", "Kendall", "CP", "FCCS", "Lin")
   rr=list()
    for(i in 1:max(obj[,1])){
      rr[[i]]=subset(obj, Trial==i)
      rr[[i]]=rr[[i]][,-1]
      #res.here[,,i]=data.matrix(rr[[i]])
    }
    for(i in 1:length(rr)){
      for(j in 1:ncol(rr[[i]]))
      res.here[,j,i]=as.numeric(as.character((rr[[i]][,j])))
    }
    # dimnames(res.here)[[3]]=paste("trials", 1:max(obj[,1]), sep="")
    allMeans=t(apply(res.here[,,], 1, rowMeans))
    colnames(allMeans)=c("Sources", "Size", "Noise", "Kendall", "CP", "FCCS", "Lin")
    allSD=apply(res.here[,4:7,], c(1,2), sd)
    colnames(allSD)=paste("SD", colnames(allMeans)[4:7], sep=".")
    plotRes=data.frame(allMeans, allSD)
    V=unique(plotRes$Size)
    method=c("Kendall", "CP", "FCCS", "Lin")
    filename=paste(filename, ".pdf", sep="")
    pdf(filename)
    for(v in 1:length(V)){
    data.here=subset(plotRes, Size==V[v])
    param=data.here[,1:4]
    #dataSel=data.frame(param, Method=rep(method[2], nrow(param)), Score=data.here[,method[2]], SD=data.here[,paste("SD", method[1], sep=".")])

    for(m in 1:length(method)){ #length(method)){
#       dataSel=rbind(dataSel, data.frame(param, Method=rep(method[m], nrow(param)), Score=data.here[,method[m]], SD=data.here[,paste("SD", method[m], sep=".")]))
	title.here=paste("Consistency score", method[m], "#nodes =", V[v], sep=" ")
        dataSel= data.frame(param, Method=rep(method[m], nrow(param)), Score=data.here[,method[m]], SD=data.here[,paste("SD", method[m], sep=".")])
        p=ggplot(dataSel, aes(x=Noise, y=Score))+geom_errorbar(aes(ymin=Score-SD, ymax=Score+SD), width=.1)+geom_point()+stat_smooth(method="lm")+facet_grid(~ Sources)
	p=p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_rect(aes(fill = Sources), xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.1)+ ggtitle(title.here)
	print(p)
    }
   } 
    dev.off()
  }
  
  plotRes2(simres.sf.ed2, "resdemo")
  plotRes2(simres.kegg.plength2, "resdemoKEGG")

  
  
simres.kegg.plength2=cbind(Trial=rep(1,nrow(simres.kegg.plength[[1]][-1,])), data.frame(simres.kegg.plength[[1]][-1,]))
for(i in 2:length(simres.kegg.plength)){
  d=data.frame(Trial=rep(i,nrow(simres.kegg.plength[[1]][-1,])), simres.kegg.plength[[1]][-1,])
  simres.kegg.plength2= rbind(simres.kegg.plength2,d )
}

  
load("/home/praveen/Paurush/COSBI/interactomes/Results/Simulation/keggWithNoise.rda")
plot.analysis(simres.plength, file="simres.plength")
rm(simres.plength)



load('/home/praveen/Paurush/COSBI/interactomes/Results/EvalInp/consistencySession.rda')
load('/home/praveen/consist_pathway_cross_talk_session.rda')

# All saved Results
setwd("/home/praveen/Paurush/COSBI/interactomes/Results/EvalInp")

save( simres.rand.ed, file="random_graph_direct_edges.rda")

names=c("Biogrid","HPRD","IntAct","STRING")
> scoreCon(mat2[[1]][-1])
$Con
[1] 0.2603503

$Div
[1] 0.2175594

> scoreCon(mat2[[1]][-2])
$Con
[1] 0.357254

$Div
[1] 0.1997251

> scoreCon(mat2[[1]][-3])
$Con
[1] 0.3535538

$Div
[1] 0.215115

> scoreCon(mat2[[1]][-4])
$Con
[1] 0.2554991

$Div
[1] 0.207418

> scoreCon(mat2[[2]][-4])
$Con
[1] 0.3427854

$Div
[1] 0.2332951

> scoreCon(mat2[[2]][-3])
$Con
[1] 0.4021828

$Div
[1] 0.2436097

> scoreCon(mat2[[2]][-2])
$Con
[1] 0.4042003

$Div
[1] 0.2233487

> scoreCon(mat2[[2]][-1])
$Con
[1] 0.3239321

$Div
[1] 0.1961911

> scoreCon(mat2[[2]])
$Con
[1] 0.3583155

$Div
[1] 0.22792

> scoreCon(mat2[[3]])
$Con
[1] 0.2376344

$Div
[1] 0.2778858

> scoreCon(mat2[[3]][-1])
$Con
[1] 0.1463896

$Div
[1] 0.2066349

> scoreCon(mat2[[3]][-2])
$Con
[1] 0.2961445

$Div
[1] 0.3021665

> scoreCon(mat2[[3]][-3])
$Con
[1] 0.2116457

$Div
[1] 0.3196685

> scoreCon(mat2[[3]][-4])
$Con
[1] 0.3168459

$Div
[1] 0.2729159



bestsetMat=data.frame(All=c(0.322,0.30, 0.35,0.237),BioGrid=c(0.586,0.255,0.324,0.14), HPRD=c(0.661, 0.255, 0.3427, 0.296), IntAct=c(0.6525, 0.353,0.402, 0.212), STRING=c(0.6705, 0.357, 0.404, 0.317)) 



pathDB=list(biocarta,spike,nci,kegg, reactome)


g1=convertIdentifiers(reactome[[1]], type="symbol")
g1=pathwayGraph(pathway=g1, edge.types=NULL)
for(i in 2:length(reactome)){
  g2=convertIdentifiers(reactome[[i]], type="symbol")
  g2=pathwayGraph(pathway=g2, edge.types=NULL)
  g1=join(g1,g2)
}


compreactome=g1

mynodes=intersect(intersect(intersect(intersect(nodes(compbiocarta),nodes(compspike)),nodes(compnci)),nodes(compKegg)),nodes(compreactome))

simframe=data.frame(b1='x', b2='y',lcs='xy', score=0)
for(i in 1:length(kegg)){
  for(j in 1:length(reactome)){
  x=LCS(names(kegg)[1], names(reactome)[2])
  if(x$LCS==NULL)
  simframe=rbind(simframe, data.frame(b1=x$a, b2=x$b, lcs=x$LCS, score=x$QSI))
  }
  }
  
  
  
  
  #### Join probability conc calculation
  computeConcProb2=function(cm){
  p1=p2=0
  rr=nrow(cm)*(nrow(cm)-1)/2
  for(i in 1:(ncol(cm)-1)){
    t=i+1
    for(j in t:ncol(cm)){
      for(s in 1:(nrow(cm)-1)){
	count1=count2=0
	u=s+1
	for(v in u:nrow(cm)){
	  if((cm[s,i] > cm[s,j]) & (cm[v,i] >= cm[v,j])){
	    count1=count1+1
	  }
	  if((cm[s,i] < cm[s,j]) & (cm[v,i] <= cm[v,j])){
	    count2=count2+1
	  }
	}
	if(count1!=0){
	  temp=count1/rr
	 # cat(temp)
	  p1=p1+log(temp)
	}
	if(count2!=0){
	  temp=count2/rr
	  p2=p2+log(temp)
	}
      }
    #print(p1)
    #print(p2)
    }
  }
  #score=exp(p1/(ncol(cm)*(ncol(cm)-1)/2))+exp(p2/(ncol(cm)*(ncol(cm)-1)/2))
  score=(exp((p1)/(ncol(cm)*(ncol(cm)-1)/2))+exp((p2)/(ncol(cm)*(ncol(cm)-1)/2)))/2
  #print(score)
  return(score)
}




kappCalc=function(cm){

}



FCCS=function(cm){
  memb_prob=fuzzyCluster(cm)[[1]]
  p=0
  for(i in 1:ncol(memb_prob)){
    p=p+prod(memb_prob[,i])
  }
  return(p)
}


fuzzyCluster=function(X, k="Default"){
    if (missing(X)) 
        stop("The data set must be given")
    if (is.null(X)) 
        stop("The data set X is empty")
    X = as.matrix(X)
    if (any(is.na(X))) 
        stop("The data set X must not contain NA values")
    if (!is.numeric(X)) 
        stop("The data set X is not a numeric data.frame or matrix")
    if(!is.numeric(k)|k!="Default")
	stop("k has to be a positiv4 integer (0 < k < n/2) or 'Default'")
    n = nrow(X)
    startU = NULL
    if(k=="Default"|k>=ceil(n/2)){
      k = ceil(n/2)
    }
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
        if (length(ana) == 0){
            ana = 0
        }
        if (ana == 1) {
            cat(" ", fill = TRUE)
            cat("If you want to run the same clustering algorithm, specify 1: ", fill = TRUE)
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
memb_prob = cc[[1]]
