##############################################################################
## PROJECT: CONSISTENCY IN GRAPHS # Scoring functions
## AUTHOR: PAURUSH PRAVEEN
## DATE: 17 SEPTEMBER 2015
## APPROACH: 
## DEPENDENCIES: 
## CONTACT: praveen@cosbi.eu
##############################################################################

library(irr) # Kendall and Kappa
library(Agreement) # CCC
library(fclust) # Fuzzy clustering
library(RBGL)
library(igraph)
library(e1071)



#########################################
### Sample graph from KEGG ##############
#########################################
load('/home/praveen/Paurush/Prior/StrPrior/Data/KEGGgraphs.rda')
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
#########################################
## Calculate Concordance probability
#########################################
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


computeConcProb=function(cm){
  col=ncol(cm)
  n=nrow(cm)
  cmat=matrix(0,nrow=n,ncol=n)
  for(i in 1:nrow(cm)){
    for(j in 1:nrow(cm)){
      S1=cm[i,]
      S2=cm[j,]
      count=0
      for(k in 1:col){
	for(l in 1:col){
	  if(S1[k] > S1[l]){
	    if(S2[k]>S2[l]){
	      count=count+1
	    }
	  }
	}
      }
      count=count/(col^2)
      cmat[i,j]=count
    }
  }
  #score=exp(p1/(ncol(cm)*(ncol(cm)-1)/2))+exp(p2/(ncol(cm)*(ncol(cm)-1)/2))
  #score=(exp((p1)/(ncol(cm)*(ncol(cm)-1)/2))+exp((p2)/(ncol(cm)*(ncol(cm)-1)/2)))/2
  #print(score)
  score=mean(c(cmat))
  return(score)
}

#########################################
## Calculate Kendall's W
#########################################
kedall_w=function(ratings=cm){
  w=kendall(ratings, correct = FALSE)
  score=w$value
  return(score)
}

#########################################
## Calculate Lin score
#########################################
calc_linCCC=function(cm){
  data=t(cm)
  score=unified.agreement(dataset=data, var=colnames(data), k=ncol(data), m=1)[[1]][1,2]
  return(score)
}
# Calculate Kohen's Kappa, only for categorical values
 kappCalc=function(cm){
 w=kappam.light(ceil(data))
 score=w$value
 return(score)
}

#########################################
## Calculate FCCS
#########################################

FCCS=function(cm){
# print(head(cm))
#   memb_prob=FKM(rbind(cm), stand=FALSE, k=ceiling(nrow(cm/2)))[[1]]
  k=max(ceiling(nrow(cm)/2),2)
  memb_prob=cmeans(cm,centers=k,iter.max=100, rate.par=0.01,verbose=FALSE,  dist = "euclidean", m=3, method='cmeans')$membership
  memb_prob[memb_prob==0]<-0.01
#   print(memb_prob)
  p=0
#   print(memb_prob)
  for(i in 1:ncol(memb_prob)){
    p=p+prod(memb_prob[,i])
  }
  rel=((1/ncol(memb_prob))^nrow(memb_prob))*ncol(memb_prob)
#   print(rel)
#   print(p)
  return(p/rel)
}

#########################################
## Fuzzy clustering requires Fclust
#########################################
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
    n = nrow(X)
    startU = NULL
    k = ceiling(n/2)
    m = 2
    ent = 1
    vp = rep(1, k)
    RS = 1
    conv = 1e-09
    maxit = 1e+6
    delta = NULL
    ana = 1
    ok = 1
    param = 0
#     while (ana == 1) {
#        
            clust = fuzzy_k_means(X, k, m, RS, startU)
#         
#         ana = scan("", n = 1)
#         if (length(ana) == 0){
#             ana = 0
#         }
#         if (ana == 1) {
#             cat(" ", fill = TRUE)
#             cat("If you want to run the same clustering algorithm, specify 1: ", fill = TRUE)
#             ok = scan("", n = 1)
#             if (length(ok) == 0) {
#                 ok = 0
#             }
#             if (ok != 1) {
#                 cat(" ", fill = TRUE)
#                 cat("You chose to run a new clustering algorithm", 
#                   fill = TRUE)
#                 cat("If you want to use the same parameter values, specify 1: ", 
#                   fill = TRUE)
#                 param = scan("", n = 1)
#                 if (length(param) == 0) {
#                   param = 0
#                 }
#             }
#             else param = 0
# #              
#         }
#     }
    return(clust)
}

# cc=fuzzyCluster(Mc[,1:(ncol(Mc)-1)])
# memb_prob = cc[[1]]
#########################################
## Old simulation function
#########################################

set.seed(123)
simulateTest2=function(n=c(10,20,30,40,50), noise=c(5,10,15,20,25,50), sources=c(3,4,5,6), graph.type="random",sgen="edges"){
resmat=matrix(0, nrow=1, ncol=7)
resmat=data.frame(resmat)
colnames(resmat)=c("Sources","Size", "Noise", "Kendall", "CP", "FCCS", "Lin")
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
      if(graph.type=="random"){
	M=1:4
	p=0.2
	ag=c(as(randomGraph(V, M, 0.2),"matrix"))
      }
      if(graph.type=="power.law"){
      ag=barabasi.game(n[i])
      ag=c(as(igraph.to.graphNEL(ag),"matrix"))
      }
      if(graph.type=="KEGG"){
	ag=sampleKEGGPathway(graphs,nacc,n)
	ag=as(ag, "matrix")
      }
      for(j in 1:length(noise)){
	#cat("Noise \t")
	#print(noise[j])
	#cat("\n")
	  if(sgen=="edges"){
	    ss=generateNoisy(ag, noise[j], source)
	    rownames(ss)=src[1:nrow(ss)]
	  }
	  if(sgen=="plength"){
	    ss=generateNoisy(ag, noise[j], source)
	    ss=plength(ss)
# 	    print(ss)
	    rownames(ss)=src[1:nrow(ss)]
	  }
	
	kend=kedall_w(ss)
 	print(kend)
	concP=computeConcProb(ss)
# 	print(concP)
	lin=calc_linCCC(ss)
# 	print(lin)
# 	print(head(ss))
# 	print(head(t(ss)))
        fccs=FCCS(ss)
	resmat=rbind(resmat, data.frame(Sources=source, Size=n[i], Noise=noise[j], Kendall=kend, CP=concP, FCCS=fccs, Lin=lin))
	cat(paste("Source =", source, "| Nodes = ", n[i], "| Noise = ", noise[j], "| Kendall = ", kend, "| CP = ", concP, "| FCCS = ", fccs,"| Lin = ", lin, "\n", sep=" "))
      }
    }
  }
#   resmat=resmat[-1,]
  return(resmat)
 }

  
#########################################
## Simulation function
#########################################
 

simulate_score_test = function(vert=c(10,20,30,40,50), noise=c(5,10,15,20,25,50), num_sources=c(3,4,5,6), evidence="plength", net_type="power.law", trials=10){
  set.seed(123)
 # resmat= data.frame(Trial=double() ,Sources=double(), Size=double(), Noise=double(),  Cor=double(), Kendall=double(), CP=double(), FCCS=double(), Lin=double())
#   dimnames(resmat)[[2]]=c("Sources","Size", "Noise", "Kendall", "CP", "FCCS", "Lin")
#   dimnames(resmat)[[3]]=paste("trial",1:trials)
  Trial=Sources=Size=Noise=Cor=Kendall=CP=FCCS=Lin=c()
  for(v in 1:length(vert)){
  n=vert[v]
  genes_name=paste("g",1:n,sep="#")
#   print(genes_name)
  for(t in 1:trials){
    #### What kind of graph 
    if(net_type=="random"){
	M=1:4
	p=0.2
	#ag=as(randomGraph(V, M, 0.2), "matrix")
	ag=matrix(1,ncol=n[no],nrow=n[no])
	colnames(ag)=rownames(ag)=genes_name
	diag(ag)=0
    }
    if(net_type=="power.law"){
      ag=barabasi.game(n, power=runif(1,2,3),directed=FALSE)
      ag=as(igraph.to.graphNEL(ag), "matrix")
      colnames(ag)=rownames(ag)=genes_name
      diag(ag)=0
    }
    if(net_type=="KEGG"){
	cat("#")
# 	ag=sampleKEGGPathway(graphs,nacc,n)
	ag=kg[[t]]
	cat("$*")
	ag=as(ag, "matrix")
	colnames(ag)=rownames(ag)=genes_name
	diag(ag)=0
    }
    for(ns in 1:length(num_sources)){
	  for(no in 1:length(noise)){
    ##### Create new networks with noise
	slist=list()
	slist[[1]]=ag
	for (so in 2:num_sources[ns]){
	  if(noise[no]==0){
	    slist[[so]] = ag
	  }else
	  slist[[so]]=permute.net2(as(ag,"graphNEL"),noise[no])
	}
#      print(slist)
    #### Path length or direct edges
	ss=matrix(0,nrow=num_sources[ns], ncol=n^2)
	colnames(ss)=paste("E",1:ncol(ss),sep="_")
	if(evidence=="edges"){
	  ss[1,]=(c(slist[[1]]))
# 		print(ss)
	  for(sss in 2:length(slist)){
	    s_here=data.frame(c(slist[[sss]]))
# 	  	print(s_here)
	    colnames(s_here)=colnames(ss)
	    ss=rbind(ss,s_here)
	  }
	}
#     print(ss)
	if(evidence=="plength"){
	  G = as(slist[[1]], "graphNEL")
# 	print(class(G))
	  pl=sp.path.compute(genelist=nodes(G), genelist2=NULL, G, verbose=FALSE)[[2]]
	  pl=1/pl
	  pl[pl==Inf]<-0
	  diag(pl)=0
	  ss[1,]=c(pl)
# 	print(ss)
	  for(sss in 2:length(slist)){
# 	  print(sss)
	    G = as(slist[[sss]], "graphNEL")
	    pl_here=sp.path.compute(genelist=nodes(G), genelist2=NULL, G, verbose=FALSE)[[2]]
	    pl_here=1/pl_here
	    pl_here[pl_here==Inf]<-0
	    diag(pl_here)=0
# 	  	print(pl_here)
	    ss[sss,]=c(pl_here)
# 	  	print(ss)
	  }
	}
# 	print(ss)
	rownames(ss)=paste("S",(1:nrow(ss)))
	kend=kedall_w(ss)
	lin=calc_linCCC(ss)
# 	print(lin)
	concP=computeConcProb(ss)
# 	print(concP)
	fccs=FCCS(ss)
	cr=0
	for(cors in 1:num_sources[ns]){
	  cr=cr+sd(c(ss[cors,]))
	}
	cor=cr/num_sources[ns]
	Trial=c(Trial,t)
	Sources=c(Sources,num_sources[ns])
	Size=c(Size,n)
	Noise=c(Noise,noise[no])
	Cor=c(Cor,cor)
	Kendall=c(Kendall, kend)
	CP=c(CP,concP)
	FCCS=c(FCCS,fccs)
	Lin=c(Lin,lin)
# 	resmat=rbind(resmat, data.frame(Trial=t ,Sources=num_sources[ns], Size=n, Noise=noise[no], Cor=cor, Kendall=kend, CP=concP, FCCS=fccs, Lin=lin))
# 	print(resmat)
	print(paste(t,num_sources[ns],n,noise[no],cor,kend,concP,fccs,lin, sep="|"))
      }
    }
    }
#     cat(paste(t,num_sources[ns],n,noise[no],cor,kend,concP,fccs,lin, "\n", sep="|"))
  }
  return(data.frame(Trial,Sources,Size,Noise,Cor,Kendall,CP,FCCS,Lin))
# cat(paste(t,num_sources[ns],n,noise[no],cor,kend,concP,fccs,lin, sep="|"))
}

# simres.ppi.plength3=simulate_score_test(n=c(10,20,30,40,50,60), noise=c(10,15,20,25,50), num_sources=c(3,4,5,6), evidence="plength", net_type="power.law", trials=5)
#########################################
## Simulation test wrapper
#########################################
wrapper_fun=function(n=c(10,20,30,40,50,60), noise=c(10,15,20,25,50), num_sources=c(3,4,5,6),evidence="plength", net_type="power.law", trials=10)  {
  res =data.frame(Trial=double() ,Sources=double(), Size=double(), Noise=double(), Kendall=double(), CP=double(), FCCS=double(), Lin=double())
  for(nd in 1:length(n)){
	res=rbind(res,simulate_score_test(n=n[nd],noise=noise,num_sources=num_sources, evidence, net_type, trials))
      }
   return(res)
}
#########################################
## Generate sources
#########################################
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


#########################################
## Permute network
#########################################
permute.net2 = function(net,noise){
	genes=nodes(net)
	LV=as(net,"matrix") 
	LV[LV >= 0.5]<-1 
        newLV=LV   
	i=sample(1:length(genes), max(1,(noise/100)*dim(LV)))
        j=sample(1:length(genes), max(1,(noise/100)*dim(LV)))
        for(x in 1:length(i)){
	  for(y in 1:length(j)){
	    if(LV[i[x],j[y]]==1){
		newLV [i[x],j[y]]= 0	
	    }
	    if(LV[i[x],j[y]]==0){
		newLV [i[x],j[y]]= 1	
	    }
	  }
	}
	return (newLV)
}
#########################################
## Compute shortest path
#########################################
sp.path.compute = function(genelist, genelist2=NULL, G, verbose=TRUE){
	require(RBGL)
	if(is.null(genelist2)){
		paths = c()
		pathL = matrix(0, nrow=length(genelist), ncol=length(genelist))
		pathLweight = matrix(0, nrow=length(genelist), ncol=length(genelist))
		for(i in 1:length(genelist)){
			for(j in 1:length(genelist)){
				if(i != j){
					path = sp.between(G, genelist[i], genelist[j])
					paths = c(paths, path)
					if(!is.na(path[[1]]$length)){
						pathL[i,j] = length(path[[1]][[3]][[1]])
						pathLweight[i,j] = path[[1]][[1]]
					}	
				          
				}
				if(verbose==TRUE)
				cat(".")
			}
			if(verbose==TRUE)
			cat(i)
		}
		res=list(paths, pathL, pathLweight)
		return(res)
	}
	else{
		paths = sp.between(G, genelist, genelist2)
		return(paths)
	}
}
#########################################
## Plot function
#########################################
 plotRes=function(obj, filename){
    res.here=array(NaN, dim=c(nrow(obj)/max(obj[,1]), ncol(obj)-1, max(obj[,1])))
    #dimnames(res.here)[[2]]=c("Sources", "Size", "Noise", "cor","Kendall", "CP", "FCCS", "Lin")
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
#########################################
## Create plot object
#########################################
compRes=function(simres.kegg.plength){
  obj_comp=cbind(Trial=rep(1,nrow(simres.kegg.plength[[1]][-1,])), data.frame(simres.kegg.plength[[1]][-1,]))
for(i in 2:length(simres.kegg.plength)){
  d=data.frame(Trial=rep(i,nrow(simres.kegg.plength[[1]][-1,])), simres.kegg.plength[[1]][-1,])
  obj_comp= rbind(obj_comp,d )
}
return(obj_comp)
}


## Simulation

simulateTest=function(n=c(10,20,30,40,50), noise=c(5,10,15,20,25,50), sources=c(3,4,5,6), graph.type="random",sgen="edges"){
resmat=matrix(0, nrow=1, ncol=8)
resmat=data.frame(resmat)
colnames(resmat)=c("Sources","Size", "Noise", "Cor", "Kendall", "CP", "FCCS", "Lin")
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
      if(graph.type=="random"){
	M=1:4
	p=0.2
	ag=c(as(randomGraph(V, M, 0.2),"matrix"))
      }
      if(graph.type=="power.law"){
      ag=barabasi.game(n[i])
      ag=c(as(igraph.to.graphNEL(ag),"matrix"))
      }
      if(graph.type=="KEGG"){
	ag=sampleKEGGPathway(graphs,nacc,n)
	ag=as(ag, "matrix")
      }
      for(j in 1:length(noise)){
	#cat("Noise \t")
	#print(noise[j])
	#cat("\n")
	  if(sgen=="edges"){
	    ss=generateNoisy(ag, noise[j], source)
	    rownames(ss)=src[1:nrow(ss)]
	  }
	  if(sgen=="plength"){
	    ss=generateNoisy(ag, noise[j], source)
	    ss=plength(ss)
# 	    print(ss)
	    rownames(ss)=src[1:nrow(ss)]
	  }
	
	kend=kedall_w(ss)
# 	print(kend)
	concP=computeConcProb(ss)
# 	print(concP)
	lin=calc_linCCC(ss)
        fccs=FCCS(ss)
        cr=0
        for(cors in 1:source){
	  cr=cr+cor(ss[cors,], method="Pearson")
        }
        cor=cr/source
	resmat=rbind(resmat, data.frame(Sources=source, Size=n[i], Noise=noise[j], Cor=cor, Kendall=kend, CP=concP, FCCS=fccs, Lin=lin))
	cat(paste("Source=", source, "| Nodes=", n[i], "| Noise=", noise[j], "| Cor=", cor,"| Kendall=", kend, "| CP=", concP, "| FCCS=", fccs,"| Lin=", lin, "\n", sep=" "))
      }
    }
  }
#   resmat=resmat[-1,]
  return(resmat)
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

## Plength calculation
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


# kg=list()
# for(i in 1:5){
# kg[[i]]=sampleKEGGPathway(graphs, nacc, 40)
# }


# Computing for MAT ppi databases
#names(set_here)=c("Biogrid","HPRD","IntAct","STRING")
run_ppi_test=function(mat,names){
lins=cp=rem=I=c()
for(i in 1:length(mat2)){
  set_here=mat2[[i]]
  ss=matrix(0,nrow=length(set_here), ncol=nrow(set_here[[1]])^2)
  for(j in 1:length(set_here)){
    ss[j,]=c(set_here[[j]])
  }
  ss=1/ss
  ss[is.na(ss)]<-0
  rownames(ss)=names
#   print(ss)
  lin=calc_linCCC(ss)
  concP=computeConcProb(ss)
  lins=c(lins,lin)
  cp=c(cp,concP)
  rem=c(rem,"all")
  I=c(I,i)
  cat(paste("Lin=",lin, "ConcP=",concP, "rem=","all", "I=",i, "\n"))
  for(k in 1:length(set_here)){
    newss=ss[-k,]
    lins=c(lins,calc_linCCC(newss))
    cp=c(cp,computeConcProb(newss))
    rem=c(rem,names(k))
    I=c(I,i)
#      print(dataframe(I,rem,cp,lins))
  }
}

return(data.frame(I,rem,cp,lins))
}

# run_ppi_test(mat2,names)


crossTalker=function(dbmat, dblist, markers=NULL){
    lin=concP=rem=c()
    neib=colnames(dbmat[[1]])
    for(i in 2:length(dbmat)){
      neib=intersect(neib, colnames(dbmat[[i]]))
    }
    ss=matrix(0,nrow=length(dbmat), ncol=length(neib)^2)
    for(j in 1:length(dbmat)){
      ss[j,]=c(dbmat[[j]][neib,neib])
    }
    ss=1/ss
    ss[is.na(ss)] <- 0
    rownames(ss)=names(dblist)
#     print(ss)
    lin=c(lin,as.numeric(calc_linCCC(ss)))
#     print(lin)
    concP=c(concP,computeConcProb(ss))
#     print(concP)
    rem=c(rem,"all")
    if(!is.null(markers)){
      com_gene=intersect(neib,markers)
      print(com_gene)
      for(k in 1:length(com_gene)){
	gl=setdiff(neib,com_gene[k])
	print(gl)
	ss=matrix(0,nrow=length(dbmat), ncol=length(gl)^2)
	rem=c(rem,com_gene[k])
	for(l in 1:length(dbmat)){
	  ss[l,]=c(dbmat[[l]][gl,gl])
	}
	ss=1/ss
	ss[is.na(ss)] <- 0
# 	print(ss)
	rownames(ss)=names(dblist)
	lin=c(lin,as.numeric(calc_linCCC(ss)))
	concP=c(concP,computeConcProb(ss))
	  
      }
    }
      
      return(list(LIN=lin, CP=concP, REM=rem))
      
    }

#  test=crossTalker(dbmat=msp_gluc, dblist=db_gluc, markers=genes)
# 
# for()
# allscores=list()
# allscores[[1]]=crossTalker(dbmat=msp_gluc, dblist=db_gluc, markers=genes)
# allscores[[2]]=crossTalker(dbmat=msp_insulin, dblist=db_insulin, markers=genes)
# allscores[[3]]=crossTalker(dbmat=msp_toll, dblist=db_toll, markers=genes) 
# allscores[[4]]=crossTalker(dbmat=msp_nfkb, dblist=db_nfkb, markers=genes) 

crossTalker2=function(dbmat, dblist, markers){
    neib=colnames(dbmat[[1]])
    for(i in 2:length(dbmat)){
      neib=intersect(neib, colnames(dbmat[[i]]))
    }
    x=1
    gn=nei=trlist1=trlist2=list()
    tab=data.frame(matrix(0, nrow=1,ncol=4))
    colnames(tab)=c("db", "marker", "cs1", "cs2")
    for(j in 1:length(dbmat)){
      hm=names(dblist)[j]
      interface=intersect(markers, colnames(dbmat[[j]]))
      ss1 = ss2 = data.frame(matrix(0,ncol=length(neib), nrow=1))
      
      for(k in 1:length(interface)){
	r=c()
	for(l in 1:length(dbmat)){
# 	print(paste("####", hm, sep=" "))
# 	print(interface[k])
# 	print(names(dblist)[l])
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
	cs1=as.numeric(calc_linCCC(ss1))
	cs2=computeConcProb(ss1)
	tab=rbind(tab, data.frame(db=names(dblist)[j], marker=interface[k], cs1=cs1, cs2=cs2))
	print(tab)
      }
      
    }
    tab=tab[-1,]
    res=list(table=tab, interface=neib, nb=nei, ge=gn, ilist1=trlist1,ilist2=trlist2)
    return(res)
}




# 
# allscores_1=list()
# allscores_1[[1]]=crossTalker2(dbmat=msp_gluc, dblist=db_gluc, markers=genes)
# allscores_1[[2]]=crossTalker2(dbmat=msp_insulin, dblist=db_insulin, markers=genes)
# allscores_1[[3]]=crossTalker2(dbmat=msp_toll, dblist=db_toll, markers=genes) 
# allscores_1[[4]]=crossTalker2(dbmat=msp_nfkb, dblist=db_nfkb, markers=genes) 
# 

# restab=data.frame(matrix(0, nrow=1,ncol=13))
#     colnames(restab)=c("Pathway", "Marker", "Frequency","CS1","CS2")
#     for(i in 1:length(allscores)){
#       mytab.here=allscores[[i]]$table
#       markers=unique(mytab.here$marker)
#       print(markers)
#       for(j in 1:length(markers)){
# 	  tab.here=as.data.frame(subset(allscores[[i]]$table, marker==markers[j]))
# 	  #print(tab.here)
# 	  #print(markers[j])
# 	  mn=apply(tab.here[,c(3:7)],2,mean)
# 	  sdv=apply(tab.here[,c(3:7)],2,sd) 
# 	  mn=t(data.frame(mn))
# 	  sdv=t(data.frame(sdv))
# 	  colnames(mn)=c("CS1","CS2","RNK1", "RNK2", "RNK3")
# 	  colnames(sdv)=c("SDCS1", "SDCS2", "SDRNK1","SDRNK2","SDRNK3")
# 	  restab=rbind(restab, data.frame(Pathway=names(allscores)[i],Marker=markers[j], Frequency=nrow(tab.here),cbind(mn,sdv)))
#       }
#     }
# 
#     
  