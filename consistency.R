##############################################################################
## PROJECT: CONSISTENCY IN PPI GRAPHS
## AUTHOR: PAURUSH PRAVEEN
## DATE: 13 MARCH 2014
## APPROACH: 
## DEPENDENCIES:
## CONTACT: praveen@cosbi.eu
##############################################################################


library(igraph)
library(RBGL)
library(ggplot2)
library(reshape)



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
  return(s1)
  }

cm=compMat

##### Concordance Probability Method Score
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
	  p1=p1+log(temp)
	}
	if(count2!=0){
	  temp=count2/rr
	  p2=p2+log(temp)
	}
      }
    }
  }
  score=(exp((p1)/(ncol(cm)*(ncol(cm)-1)/2))+exp((p2)/(ncol(cm)*(ncol(cm)-1)/2)))/2
  pv=friedTest(cm)
  res=list(Score=score, p.value=pv)#, conc=pv$conc)
  return(res)
}

#### Friedman text function

friedTest=function(cm){
  x=friedman.test(t(cm))
# z=t(cm)
# group=c(1:ncol(z))
# z=kendall.post(z, group, nperm = 999, mult = "holm")
  y=x$p.value
# res=list(conc=z, pv=y)
  return(y)
 }



#### Variance function
computeVarScore=function(cm){
  s=mean(apply(cm,2,sd))
  s=1/s
  return(s)
}


#### Compute correlation matrix
computeCorScore=function(cm){
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

###### Noisy network
generateNoisy=function(g, noise, source){
  ss=data.frame(matrix(0,ncol=length(g), nrow=1))
  colnames(ss)=paste("E", 1:length(g), sep="_")
  #gg=c(as(g, "matrix"))
  for(so in 1:source){
  if(noise==0){
    gh = g
    temp = matrix(gh, nrow=1, ncol=ncol(ss))
    colnames(temp)=colnames(ss)
    ss = rbind(ss, temp)
  }
  else{
#     print(so)
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
  }
  ss=ss[-1,]
  return(ss)
}



 
                        
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


#### Generate adjacency network
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

#### Compute shortest path length

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

#### Compute shortest path matrices
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
  if(length(index)>1){
    for(i in 1:length(mylist)){
      G=join(G,pathwayGraph(pathway=convertIdentifiers(mylist[[i]], type="symbol"),edge.types=NULL))
    }
  }
return(G)
}

#### Function for cross-talk

crossTalker=function(dbmat, dblist, markers){
    neib=colnames(dbmat[[1]])
    for(i in 2:length(dbmat)){
      neib=intersect(neib, colnames(dbmat[[i]]))
    }
    x=1
    gn=nei=trlist1=trlist2=list()
    tab=data.frame(matrix(0, nrow=1,ncol=11))
    colnames(tab)=c("db", "marker", "cs1", "cs2", "conc1","conc2","pv1","pv2","rank1", "rank2", "rank3")
    for(j in 1:length(dbmat)){
      cat(j)
      hm=names(dblist)[j]
      interface=intersect(markers, colnames(dbmat[[j]]))
      ss1 = ss2 = data.frame(matrix(0,ncol=length(neib), nrow=1))
      for(k in 1:length(interface)){
	r=c()
	for(l in 1:length(dbmat)){
	#print(paste("####", hm, sep=" "))
	#print(interface[k])
	#print(names(dblist)[l])
	#cat("\n")
	if(interface[k] %in% nodes(dblist[[l]])){
	  r=c(r,betweenness(graph=igraph.from.graphNEL(dblist[[l]]), v = interface[k], directed = TRUE))#r=c(r,page.rank(graph=igraph.from.graphNEL(dblist[[l]]), algo ="prpack",vids = interface[k], directed = TRUE))$vector
	}else{
	  r=c(r,0)
	}
	if(interface[k] %in% colnames(dbmat[[l]])){
	cat(".")
	gn[[x]]=interface[k]
	nei[[x]]=neib
	tr1=as.matrix(dbmat[[l]][interface[k],neib]) # Outgoing edges
	tr2=as.matrix(dbmat[[l]][neib,interface[k]]) # Incoming edges
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
	trlist1[[x]] = ss1
	trlist2[[x]] = ss2
	x=x+1
	cs1=computeConcProb(ss1)
# 	print(cs1)
	cs2=computeConcProb(ss2)
# 	print(cs2)
	kconc1=computeKendalConc(ss1)
# 	print(kconc1)
	kconc2=computeKendalConc(ss2)
# 	print(kconc2)
	tab=rbind(tab, data.frame(db=names(dblist)[j], marker=interface[k], cs1=cs1$Score, cs2=cs2$Score, conc1=kconc1,conc2=kconc2,pv1=cs1$p.value,pv2=cs2$p.value,rank1=r[1], rank2=r[2], rank3=r[3]))
      }
      cat("\n")
    }
    tab=tab[-1,]
    res=list(table=tab, interface=neib, nb=nei, ge=gn, ilist1=trlist1,ilist2=trlist2)
    return(res)
}
 
dd=crossTalker(dbmat=msp_gluc, dblist=db_gluc, markers=genes)
dd1=crossTalker(dbmat=msp_toll, dblist=db_toll, markers=genes)
dd2=crossTalker(dbmat=msp_insulin, dblist=db_insulin, markers=genes)
dd3=crossTalker(dbmat=msp_nfkb, dblist=db_nfkb, markers=genes)
### implement kendall concordance
computeKendalConc=function(mat){
  sources=nrow(mat)
  dpoint=ncol(mat)
  R_i=mat
  for(i in 1:sources){
    R_i[i,]=rank(1/mat[i,],  ties.method = "average")
    
  }
  R_i=apply(R_i,2,sum)
  R_i_bar=sum(R_i)/dpoint
  S=0
  for(j in 1:dpoint){
    S=S+(R_i[j]-R_i_bar)^2
  }
  W=12*S/((sources^2)*((dpoint^3)-dpoint))
  return(W)
}


