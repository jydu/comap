# +---------------------------------------------------------+
# | CoMap p-value computation -- R functions                |
# | Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 28/06/06 |
# | Modified on: 15/01/2007                                 |
# | Use a sliding window and conditional p-values           |
# +---------------------------------------------------------+

# Get p-values codes:
pval<-function(x)
{
  return(as.character(symnum(x, cutpoints=c(0,0.001,0.01,0.05,0.1,1), symbols=c("***","**","*",".","NS"))))
}

round.pval<-function(p)
{
  for(i in 1:length(p))
  {
    if(!is.na(p[i])) if(p[i]<0) p[i]<-0
  }
  return(p)
}

# Test each group.
# Performs type I test (see article).
# data: The clustered groups (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# sim: The simulation result (data.frame).
# window: The size of the sliding window, in percent.
# Return a data.frame similar to data, with a new column named 'p.value'.
test<-function(data, sim, group.sizes, window)  
{ 
  cn<-"p.value"
  data[,cn]<-rep(NA,nrow(data))
  for(i in group.sizes)
  {
    data.group<-data[data$Size==i,]
    if(nrow(data.group) > 0)
    {
      sim.group<-sim[sim$Size==i,]
      ws<-(max(sim.group$Nmin) - min(sim.group$Nmin)) * window / 2

      for(j in 1:nrow(data.group))
      {
        nmin<-data.group[j,"Nmin"]
        dmax<-data.group[j,"Dmax"]
        d<-sim.group[sim.group$Nmin > nmin - ws & sim.group$Nmin < nmin + ws,]
        data[which(data$Group==data.group[j, "Group"]),cn]<-(sum(d$Dmax <= dmax) + 1) / (nrow(d) + 1)
      }
    }
  }
  data$code<-pval(data$p.value)
  return(data)
}

# Test each group.
# Performs type II test (see article).
# data: The clustered groups (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# sim: The simulation result (data.frame).
# window: The size of the sliding window for Nmin, in percent.
# window2: The size of the sliding window for delta, in percent.
# Return a data.frame similar to data, with a new column named 'p.value'.
test2<-function(data, sim, group.sizes, window, window2)
{ 
  cn<-"p.value"
  data[,cn]<-rep(NA,nrow(data))
  for(i in group.sizes)
  {
    data.group<-data[data$Size==i,]
    if(nrow(data.group) > 0)
    {
      sim.group<-sim[sim$Size==i,]
      ws<-(max(sim.group$Nmin) - min(sim.group$Nmin)) * window / 2
      ws2<-(max(sim.group$Delta) - min(sim.group$Delta)) * window2 / 2

      for(j in 1:nrow(data.group))
      {
        nmin<-data.group[j,"Nmin"]
        delta<-data.group[j,"Delta"]
        dmax<-data.group[j,"Dmax"]
        d<-sim.group[sim.group$Nmin > nmin - ws & sim.group$Nmin < nmin + ws
                   & sim.group$Delta > delta - ws2 & sim.group$Delta < delta + ws2,]
        data[which(data$Group==data.group[j, "Group"]),cn]<-(sum(d$Dmax <= dmax) + 1) / (nrow(d) + 1)
      }
    }
  }
  data$code<-pval(data$p.value)
  return(data)
}

# when nested groups are detected, keep only the level with the most significant value
belongsto<-function(group1, group2)
{
  s1<-strsplit(substr(group1, start=2, stop=nchar(group1)-1),";")[[1]]
  res<-logical(length(group2))
  for(i in 1:length(group2))
  {
    s2<-strsplit(substr(group2[i], start=2, stop=nchar(group2[i])-1),";")[[1]]
    res[i]<-all(s1 %in% s2)
  }
  return(res)
}
size<-function(groups)
{
  sizes<-numeric(length(groups))
  for(i in 1:length(groups))
  {
    s<-strsplit(substr(groups[i], start=2, stop=nchar(groups[i])-1),";")[[1]]
    sizes[i]<-length(s)
  }
  return(sizes)
}
build.cliques<-function(pred,cng,logFile="")
{
  cat("Building cliques...\n")
  cat(file=logFile,date(),"\n")
  cliques<-as.list(as.character(pred$Group))
  names(cliques)<-as.character(pred$Group)

  #Now sort cliques:

  groups<-names(cliques)
  groups<-groups[order(pred$Size,pred$p.value)]
  while(length(groups)>0)
  {
    group<-groups[1]
    candidates<-which(belongsto(group,groups))
    candidates<-candidates[-1]
    if(length(candidates) > 0)
    {
      superGroups<-groups[candidates]
      superGroup<-superGroups[which.min(size(superGroups))]
      if(cng)
      {
        #Check groups p-values:
        groupPVal     <-pred[pred$Group==group     ,"p.value"]
        superGroupPVal<-pred[pred$Group==superGroup,"p.value"]
        #This keeps the most significant group:
        if(groupPVal <= superGroupPVal)
        {
          #Small group is better, remove big group:
          cat("Removing group",superGroup,"[p-value",superGroupPVal,"] for group",group,"[p-value",groupPVal,"]\n",file=logFile,append=TRUE)
          cliques[[superGroup]]<-NULL
          groups<-groups[-which(groups==superGroup)]
          #... and test again current group!
        }
        else
        {
          #Big group is better, remove small group:
          cat("Removing group",group,"[p-value",groupPVal,"] for group",superGroup,"[p-value",superGroupPVal,"]\n",file=logFile,append=TRUE)
          cliques[[group]]<-NULL
          groups<-groups[-1]
        }
        #This keeps the smallest group:
#        cat("Removing group",superGroup,"[p-value",superGroupPVal,"] for group",group,"[p-value",groupPVal,"]\n",file=logFile,append=TRUE)
#        cliques[[superGroup]]<-NULL
#        groups<-groups[-which(groups==superGroup)]
#        #... and test again current group!
      }
      else
      {
        cat("Merging group",group,"with",superGroup,"[",length(groups),"remaining ]\n",file=logFile,append=TRUE)
        # Merging...
        cliques[[superGroup]]<-append(cliques[[superGroup]],cliques[[group]])
        # Removing old group:
        cliques[[group]]<-NULL
        groups<-groups[-1]
      }
    }
    else
    {
      #Otherwise lonely clique
      cat("Group",group,"has no super clique.\n",file=logFile,append=TRUE)
      groups<-groups[-1]
    }
  }

  # Now sort table:
  groups<-character(0)
  cliquesId<-numeric(0)
  count<-1
  for(i in names(cliques))
  {
    clique<-cliques[[i]]
    groups<-c(groups,clique)
    cliquesId<-c(cliquesId,rep(count,length(clique)))
    count<-count+1
  }

  rownames(pred)<-as.character(pred$Group)
  pred<-pred[groups,]
  pred$Clique<-cliquesId
  return(pred)
}

# Test all groups, and format results.
# data: The clustered groups (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# sim: The simulation result (data.frame).
# window: The size of the sliding window, in percent.
# method: Add a column with the name of the method used (eg Volume, Charge, etc.).
# level: The maximum p-value level for groups to output.
# cng: Tell if correction for nested groups must be applied.
# logFile: Where to write the groups removed.
get.pred<-function(data, sim, group.sizes, window, method="", level=0.05, cng, logFile)
{
  pred<-test(data, sim, group.sizes, window)
  pred<-pred[!is.na(pred$p.value),]
  pred$p.value<-round.pval(pred$p.value)
  pred<-pred[pred$p.value <= level & pred[,"Const"] == "no",]
  if(nrow(pred)==0) return(pred)
  pred<-pred[order(pred$p.value),]
  if(method != "")
  {
    pred[,"Method"]<-rep(method, nrow(pred))
  }
  n1<-sum(data$Size %in% group.sizes & data$Const == "no")
  n2<-nrow(pred)
  pred<-build.cliques(pred,cng,logFile)
  n3<-nrow(pred)
  p1<-sum(dbinom(n2:n1,n1,level))
  p2<-sum(dbinom(n3:n1,n1,level))
  cat(n1, "tested,", n2, "sign.,", n3, "indep. p-value in [", p1, ",", p2, "].\n")
  return(pred)
}

# Test all groups, and format results.
# data: The clustered groups (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# sim: The simulation result (data.frame).
# window: The size of the sliding window for Nmin, in percent.
# window2: The size of the sliding window for delta, in percent.
# method: Add a column with the name of the method used (eg Volume, Charge, etc.).
# level: The maximum p-value level for groups to output.
# cng: Tell if correction for nested groups must be applied.
# logFile: Where to write the groups removed.
get.pred2<-function(data, sim, group.sizes, window, window2, method="", level=0.05, cng, logFile)
{
  pred<-test2(data, sim, group.sizes, window, window2)
  pred<-pred[!is.na(pred$p.value),]
  pred$p.value<-round.pval(pred$p.value)
  pred<-pred[pred$p.value <= level & pred[,"Const"] == "no",]
  if(nrow(pred)==0) return(pred)
  pred<-pred[order(pred$p.value),]
  if(method != "")
  {
    pred[,"Method"]<-rep(method, nrow(pred))
  }
  n1<-sum(data$Size %in% group.sizes & data$Const == "no")
  n2<-nrow(pred)
  pred<-build.cliques(pred,cng,logFile)
  n3<-nrow(pred)
  p1<-sum(dbinom(n2:n1,n1,level))
  p2<-sum(dbinom(n3:n1,n1,level))
  cat(n1, "tested,", n2, "sign.,", n3, "indep. p-value in [", p1, ",", p2, "].\n")
  return(pred)
}

