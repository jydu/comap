# +---------------------------------------------------------+
# | CoMap p-value computation -- R functions                |
# | Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 28/06/06 |
# | Modified on: 15/01/2007                                 |
# | Use a sliding window and conditional p-values           |
# | Modified on: 10/07/2009                                 |
# | Compute threshold p-value for a given FDR               |
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

progress<-function(i=0,max=NA,size=50)
{
  if(i < 1)
  {
    cat("|",paste(rep("-",size),collapse=""),"|\n",sep="")
  }
  else
  {
    if(i>max) return()
    if(i == 1) cat("|")
    step<-ceiling(max/size)
    x<-i %% step
    if(x == 0) cat("=");
    if(i == max) cat(paste(rep("=", size-floor(max/step)), collapse=""));
    if(i == max) cat("|", fill=TRUE)
  }
}

# Test each group.
# Performs type I test (see article).
# data: The clustered groups (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# sim: The simulation result (data.frame).
# window: The size of the sliding window, in percent.
# grid.Rate: precompute a grid of size round(1/window). This will speed up the computation, but may be unappropriate if there are not enough simulations.
# grid.Stat: also precompute a grid for Statistics, for certain p-values (not fully tested yet, you should not use this option!).
# verbose: Tell if progress bars and messages should be output.
# Return a data.frame similar to data, with a new column named 'p.value'.
test<-function(data, sim, group.sizes, window, min.nobs, grid.Rate=FALSE, grid.Stat=FALSE, statName="Stat", rateName="Nmin", lower=FALSE, verbose=TRUE)
{ 
  data[,"p.value"]<-rep(NA,nrow(data))
  data[,"nobs"]<-rep(NA,nrow(data))
  for(i in group.sizes)
  {
    group.filter<-which(data$Size==i)
    tot<-length(group.filter)
    if(tot > 0)
    {
      sim.group<-sim[sim$Size == i,]
      if(grid.Rate)
      #We discretize rate for speeding up computations
      {
        if(verbose) cat("Computing grid for group size", i,".\n") 
        grid.size<-round(1/window)
        ma<-max(sim.group[,rateName])
        mi<-min(sim.group[,rateName])
        grid.bounds<-mi + (0:grid.size)*(ma - mi)/grid.size
        grid.f<-cut(sim.group[, rateName], breaks=grid.bounds, labels=FALSE)
        grid<-split(sim.group[, statName], grid.f)
        if(grid.Stat)
        #We also discretize Stat for speeding up computations
        {
          grid.stat<-c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,5e-2,1e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1,8e-1,9e-1,1)
          if(!lower) grid.stat<-sort(1-grid.stat)
          nobs<-list()
          pval<-list()
          for(j in as.character(1:length(unique(grid.f))))
          {
            d<-grid[[j]]
            nobs[[j]]<-length(d)
            if(length(d) < min.nobs)
            {
              cat("Rate category", j, "will be ignored because of unsuficient number of points.\n")
            }
            else
            {      
              d<-sort(d)
              g<-unique(round(grid.stat * (length(d) + 1) - 1))
              if(lower) pval[[j]]<-grid.stat[g > 0]
              else pval[[j]]<-1. - grid.stat[g > 0]
              g<-g[g > 0]
              grid[[j]]<-c(-Inf, d[g])
            }
          }
          if(verbose) progress()

          tmp.pval<-numeric(tot)
          tmp.nobs<-numeric(tot)
          tmp.indx<-numeric(tot)
          for(j in 1:tot)
          {
            progress(j, tot)
            nmin<-as.character(findInterval(data[group.filter[j], rateName], grid.bounds, rightmost.closed=TRUE))
            d<-grid[[nmin]]
            n<-nobs[[nmin]]
            n<-length(d)
            p<-pval[[nmin]]
            stat<-data[group.filter[j], statName]
            if(n < min.nobs)
            {
              tmp.pval[j]<-NA
            }
            else
            {
              tmp.pval[j]<-p[findInterval(stat, d, rightmost.closed=TRUE)]
            } 
            tmp.nobs[j]<-n
            tmp.indx[j]<-group.filter[j]
          }
          #Now update data:
          data[tmp.indx, c("p.value", "nobs")]<-cbind(tmp.pval, tmp.nobs)
        }
        else
        # We do not descretize Stat and compute exact p-values.
        {
          for(j in as.character(1:length(unique(grid.f))))
          {
            d<-grid[[j]]
          }
          if(verbose) progress()

          tmp.pval<-numeric(tot)
          tmp.nobs<-numeric(tot)
          tmp.indx<-numeric(tot)
          for(j in 1:tot)
          {
            if(verbose) progress(j, tot)
            nmin<-as.character(findInterval(data[group.filter[j], rateName], grid.bounds, rightmost.closed=TRUE))
            d<-grid[[nmin]]
            n<-length(d)
            stat<-data[group.filter[j], statName]
            if(n < min.nobs)
            {
              tmp.pval[j]<-NA
            }
            else
            {
              if(lower)
              {
                tmp.pval[j]<-(sum(d <= stat) + 1) / (length(d) + 1)
              }
              else
              {
                tmp.pval[j]<-(sum(d >= stat) + 1) / (length(d) + 1)
              }
            }   
            tmp.nobs[j]<-n
            tmp.indx[j]<-group.filter[j]
          }
          #Now update data:
          data[tmp.indx, c("p.value", "nobs")]<-cbind(tmp.pval, tmp.nobs)
        }
      }
      else
      #This is the exact procedure, without any discretization.
      {
        if(verbose) cat("Analysing groups of size", i, "\n")
        ws<-(max(sim.group[,rateName]) - min(sim.group[,rateName])) * window / 2
        if(verbose) progress()
        
        tmp.pval<-numeric(tot)
        tmp.nobs<-numeric(tot)
        tmp.indx<-numeric(tot)
        for(j in 1:tot)
        {
          if(verbose) progress(j, tot)
          nmin<-data[group.filter, rateName][j]
          d<-sim.group[sim.group[,rateName] > nmin - ws & sim.group[,rateName] < nmin + ws,]
          if(nmin < 0.01)
          {
            #data[group.filter[j], "p.value"]<-1. #Threshold for conserved sites.
            tmp.pval[j]<-1.
          }
          else
          {
            stat<-data[group.filter, statName][j]
            if(nrow(d) < min.nobs)
            {
              #data[group.filter[j], "p.value"]<-NA
              tmp.pval[j]<-NA
            }
            else
            {
              if(lower)
              {
                #data[group.filter[j], "p.value"]<-(sum(d[, statName] <= stat) + 1) / (nrow(d) + 1)
                tmp.pval[j]<-(sum(d[, statName] <= stat) + 1) / (nrow(d) + 1)
              }
              else
              {
                #data[group.filter[j], "p.value"]<-(sum(d[, statName] >= stat) + 1) / (nrow(d) + 1)
                tmp.pval[j]<-(sum(d[, statName] >= stat) + 1) / (nrow(d) + 1)
              }
            } 
          }
          #data[group.filter[j], "nobs"]<-nrow(d)
          tmp.nobs[j]<-nrow(d)
          tmp.indx[j]<-group.filter[j]
        }
        data[tmp.indx, c("p.value", "nobs")]<-cbind(tmp.pval, tmp.nobs)
      }
    }
  }
  data$code<-pval(data$p.value)
  return(data)
}

# when nested groups are detected, keep only the level with the most significant value

# group1 belongs to group2? group2 may be a vector of groups.
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
# cng: correction for nested groups (deprecated).
# If one group contains a smaller, most significant group, we remove it from the list.
build.cliques<-function(pred,cng,logFile="")
{
  cat("Building cliques...\n")
  if(!is.null(logFile))
    cat(file=logFile,date(),"\n")
  cliques<-as.list(as.character(pred$Group))
  names(cliques)<-as.character(pred$Group)

  #Now sort cliques:

  groups<-names(cliques)
  groups<-groups[order(pred$Size,pred$p.value)]
  while(length(groups)>0)
  {
    group<-groups[1]
    # 'candidates' contains all supergroups of this one:
    candidates<-which(belongsto(group,groups))
    # The groups belongs to himself... ignore this.
    candidates<-candidates[-1]
    if(length(candidates) > 0)
    {
      superGroups<-groups[candidates]
      # Test the smaller supergroup:
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
          if(!is.null(logFile))
          {
            cat("Removing group",superGroup,"[p-value",superGroupPVal,"] for group",group,"[p-value",groupPVal,"]\n",file=logFile,append=TRUE)
          }
          cliques[[superGroup]]<-NULL
          groups<-groups[-which(groups==superGroup)]
          #... and test again current group!
        }
        else
        {
          #Big group is better, remove small group:
          if(!is.null(logFile))
          {
            cat("Removing group",group,"[p-value",groupPVal,"] for group",superGroup,"[p-value",superGroupPVal,"]\n",file=logFile,append=TRUE)
          }
          cliques[[group]]<-NULL
          groups<-groups[-1]
        }
#        #This keeps the smallest group:
#        cat("Removing group",superGroup,"[p-value",superGroupPVal,"] for group",group,"[p-value",groupPVal,"]\n",file=logFile,append=TRUE)
#        cliques[[superGroup]]<-NULL
#        groups<-groups[-which(groups==superGroup)]
#        #... and test again current group!
      }
      else
      {
        if(!is.null(logFile))
        {
          cat("Merging group",group,"with",superGroup,"[",length(groups),"remaining ]\n",file=logFile,append=TRUE)
        }
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
      if(!is.null(logFile))
      {
        cat("Group",group,"has no super clique.\n",file=logFile,append=TRUE)
      }
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

# New correction for nested groups:
# erase nested (!)
ernest<-function(pred, logFile="", verbose=TRUE)
{
  tmp<-unique(pred$Size)
  if(length(tmp)==1 && tmp==2) return(pred) #Nothing nested here (pairwise analysis).

  if(verbose) cat("Building cliques...\n")
  if(!is.null(logFile))
    cat(file=logFile,date(),"\n")
  groups<-as.character(pred$Group)
  row.names(pred)<-groups
  groups<-groups[order(pred$Size)]
  
  # First remove all groups that contain a smaller and most significant group.
  # We rely on while loops since the size of the 'groups' vector may change during iterations, and
  # hence need to be reevaluated each time.
  n<-length(groups)
  i<-1
  while(i < n)
  {
    group<-groups[i]
    groupPVal<-pred[group,"p.value"]
    j<-i+1
    while(j <= n)
    {
      #cat("1!",i,j,n,"\n")
      superGroup<-groups[j]
      superGroupPVal<-pred[superGroup,"p.value"]
      if(belongsto(group,superGroup) & groupPVal < superGroupPVal)
      {
        # Remove supergroup:
        groups<-groups[-j]
        n<-n-1
        if(!is.null(logFile))
        {
          cat("Removing group",superGroup,"[p-value",superGroupPVal,"] for group",group,"[p-value",groupPVal,"]\n",file=logFile,append=TRUE)
        }
      }
      else
      {
        j<-j+1
      }
    }
    i<-i+1
  }

  # Then keep only the most significant size:
  i<-length(groups)
  while(i > 1)
  {
    superGroup<-groups[i]
    superGroupPVal<-pred[superGroup,"p.value"]
    j<-i-1
    while(j >= 1)
    {
      #cat("2!",i,j,length(groups),"\n")
      group<-groups[j]
      groupPVal<-pred[group,"p.value"]
      if(belongsto(group,superGroup))
      {
        # Remove subgroup:
        groups<-groups[-j]
        i<-i-1
        if(!is.null(logFile))
        {
          cat("Removing group",group,"[p-value",groupPVal,"] for group",superGroup,"[p-value",superGroupPVal,"]\n",file=logFile,append=TRUE)
        }
      }
      j<-j-1
    }
    i<-i-1
  }
  return(pred[groups,])
}

# Test all groups.
# data: The clustered groups (data.frame).
# sim: The simulation result (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# window: The size of the sliding window, in percent.
# verbose: Tell if progress bars and messages should be output.
get.pred<-function(data, sim, group.sizes, window, min.nobs, verbose=TRUE, ...)
{
  pred<-test(data, sim, group.sizes, window, min.nobs, grid.Rate=FALSE, verbose=verbose, ...)
  pred<-pred[!is.na(pred$p.value),]
  if (!is.null(pred$Const)) pred<-pred[as.character(pred$Const) == "no",]
  if (nrow(pred) > 0)
    pred$p.value<-round.pval(pred$p.value)
  return(pred)
}

# Correction for multiple testing:
fdrcalc<-function(sim, fdr, simindex, group.sizes, window, min.nobs, cng, ...)
{
  sim<-sim[sim$Size %in% group.sizes,]
  p<-numeric(0)
  progress()
  count<-1
  for(i in simindex)
  {
    progress(count, length(simindex))
    count<-count+1
    sima<-sim[sim$Rep + 1 == i,]
    simr<-sim[sim$Rep + 1 != i,]
    pred<-get.pred(sima, simr, group.sizes, window, min.nobs, verbose=FALSE, ...)
    if(cng) pred<-ernest(pred, NULL, verbose=FALSE)
    p<-c(p,pred$p.value)
  }
  p<-sort(p)
  return(list(threshold=p[round(length(p)*fdr)], p.values=p))
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
format.pred<-function(data, sim, group.sizes, window, min.nobs, method="", level=0.05, cng, logFile, fdr=0.05, nfdr=10, ...)
{
  cat("Computing p-values for data set...\n");
  pairs<-FALSE
  if (is.null(data$Size)) #For pairwise analysis. 
  {
    if(!is.null(sim$Size) && (length(unique(sim$Size)) > 1 || unique(sim$Size)  != 2))
    {
       error("Data are from a pairwise analysis, but simulations result from a clustering analysis.")
    }
    data$Size<-2
    sim$Size<-2
    group.sizes<-2
    pairs<-TRUE
  }
  pred<-get.pred(data, sim, group.sizes, window, min.nobs, ...)
  nbtests<-sum(pred$nobs >= min.nobs)
  if (!is.na(level))
    pred<-pred[!is.na(pred$p.value) & pred$p.value <= level,]
  if(nrow(pred)==0) return(pred)
  pred<-pred[order(pred$p.value),]
  if(method != "")
  {
    pred[,"Method"]<-rep(method, nrow(pred))
  }
  n1<-sum(data$Size %in% group.sizes)
  n2<-nrow(pred)
  #cat(n2,n1,level,"\n")
  p1<-sum(dbinom(n2:n1,n1,level))
  cat(n1, "tested,", n2, "sign., global p-value >", p1, "].\n")
  if(cng)
  {
    pred<-ernest(pred,logFile)
    n3<-nrow(pred)
    #cat(n3,n1,level,"\n")
    p2<-sum(dbinom(n3:n1,n1,level))
    cat(n1, "tested,", n3, "indep., global p-value <", p2, "].\n")
  }
  if(!is.na(fdr))
  {
    cat("Computing FDR...\n");
    if(pairs)
    {
      cat("Only pairs are tested, perform Benjamini and Hochberg FDR procedure.\n")
      x<-sort(pred[,"p.value"])
      test <- which(x <= ((1:length(x)) * fdr/nbtests))
      t<-ifelse(length(test) == 0, 0, x[max(test)])
      cat("Significance threshold at level", fdr, "is", t, "(", nbtests, "performed )\n")
      pred$FDR<-ifelse(pred$p.value <= t, "yes", "no")
    }
    else
    {
      cat("Groups result from a clustering approach, perform Dutheil and Galtier FDR procedure.\n")
      
      f<-fdrcalc(sim, fdr, 1:nfdr, group.sizes, window, min.nobs, cng, ...)
      pred$FDR<-ifelse(pred$p.value <= f$threshold, "yes", "no")
    }
    n4<-sum(pred$FDR == "yes")
    cat(n4, "groups remain significant after correction for multiple testing.\n")
  }
  return(pred)
}

merge2<-function(d1, d2, ...)
{
  if (is.null(d1)) return(d2)
  if (is.null(d2)) return(d1)
  return(merge(d1, d2, ...))
}

