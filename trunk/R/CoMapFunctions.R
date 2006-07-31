# +---------------------------------------------------------+
# | CoMap p-value computation -- R functions                |
# | Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 28/06/06 |
# +---------------------------------------------------------+


# Return the position of the maximum value of a multi-dimensional array
# t: the array to parse
which.max2 <- function(t)
{
  return(which(t==max(t), arr.ind=TRUE))
}

# Get p-values codes:
pval<-function(x)
{
  return(as.character(symnum(x, cutpoints=c(0,0.001,0.01,0.05,0.1,1), symbols=c("***","**","*",".","NS"))))
}

# Distribution discretization. Perform type I test (see article).
#
# sim: The simulation result (data.frame).
# group.size: Group size for which the distribution must be computed.
# n, n.nmin, n.dmax: number of classes to use for discretization for each axe.
# Returns a n.nmin times n.dmax array, with -1 values in non-significant cells.
disc.dist<-function(sim, group.size, n=10, n.nmin=n, n.dmax=n)
{
  sim<-sim[sim$Size==group.size,]
  d.lim<-seq(0,max(sim$Dmax),len=n.dmax+1)
  n.lim<-seq(0,max(sim$Nmin),len=n.nmin+1)
  dmax.f<-cut(sim$Dmax, breaks=d.lim)
  nmin.f<-cut(sim$Nmin, breaks=n.lim)
  t<-table(dmax.f, nmin.f)
  t<-t/sum(t)
  dimnames(t)[[1]]<-d.lim[-1]
  dimnames(t)[[2]]<-n.lim[-1]
  pvalues<-t*0 - 1 #A table with the same dimensions but filled with -1.
  a<-0;
  while(max(t) >= 0) {
    #cat(a,"\n")
    d<-which.max2(t)
    i<-d[1,1]
    j<-d[1,2]
    a<-a+t[i,j]
    t[i,j]<- -1
    for(ii in i:n.dmax) {
      if(t[ii,j] > -1) {
        a<-a+t[ii,j]
        t[ii,j]<- -1
      }
    }
    pvalues[i,j]<-1-a
    for(ii in i:n.dmax) {
      if(pvalues[ii,j] == -1) {
        pvalues[ii,j]<-1-a
      }
    }
  }
  return(pvalues)
}

# Test each group.
# Performs type I test (see article).
# data: The clustered groups (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# sim: The simulation result (data.frame).
# n, n.nmin, n.dmax: number of classes to use for discretization for each axe.
# Return a data.frame similar to data, with a new column named Pred 'level', telling if a group is significant at the given level.
test<-function(data, sim, group.sizes, n=10, n.nmin=n, n.dmax=n)  
{ 
  cn<-"p.value"
  data[,cn] <- numeric(nrow(data))
  # Compute dist for each group size:
  dist<-list();
  for(i in group.sizes) {
    dist[[as.character(i)]]<-disc.dist(sim,group.size=i,n.nmin=n.nmin,n.dmax=n.dmax)
  }
  
  for(i in 1:nrow(data)) {
    group<-data[i,"Size"]
    nmin <-data[i,"Nmin"]
    dmax <-data[i,"Dmax"]
    nmin.i<-which(nmin<=as.numeric(colnames(dist[[as.character(group)]])))[1]
    dmax.i<-which(dmax<=as.numeric(rownames(dist[[as.character(group)]])))[1]
    if(is.null(dist[[as.character(group)]])) {
      data[i,cn]<-NA
    } else {
      if(is.na(nmin.i)) { # Value is out of simulation range.
        warning(paste("Nmin value out of range for group",i,"."))
        data[i,cn]<-NA
      } else {
        data[i,cn]<-dist[[as.character(group)]][dmax.i,nmin.i]
      }
    }
  }
  data$code<-pval(data$p.value)
  return(data)
}

# Distribution discretization. Perform type II test (see article).
#
# sim: The simulation result (data.frame).
# group.size: Group size for which the distribution must be computed.
# n, n.nmin, n.dmax, n.delta: number of classes to use for discretization for each axe.
# Returns a n.nmin times n.dmax array, with -1 values in non-significant cells.
disc.dist2<-function(sim, group.size, n=10, n.nmin=n, n.dmax=n, n.delta=n)
{
  sim<-sim[sim$Size==group.size,]
  d.lim<-seq(0,max(sim$Dmax),len=n.dmax+1)
  n.lim<-seq(0,max(sim$Nmin),len=n.nmin+1)
  m.lim<-seq(0,max(sim$Delta),len=n.delta+1)
  dmax.f<-cut(sim$Dmax, breaks=d.lim)
  nmin.f<-cut(sim$Nmin, breaks=n.lim)
  delta.f<-cut(sim$Delta, breaks=m.lim)
  t<-table(dmax.f, nmin.f, delta.f)
  t<-t/sum(t)
  r<-array(rank(t), dim=dim(t))
  dimnames(t)[[1]]<-d.lim[-1]
  dimnames(t)[[2]]<-n.lim[-1]
  dimnames(t)[[3]]<-m.lim[-1]
  pvalues<-t*0 - 1 #A table with the same dimensions but filled with -1.
  a<-0;
  while(max(t) >= 0) {
    #cat(a,"\n")
    d<-which.max2(t)
    i<-d[1,1]
    j<-d[1,2]
    k<-d[1,3]
    a<-a+t[i,j,k]
    t[i,j,k]<- -1
    for(ii in i:n.dmax) {
      for(kk in 1:k) {
        if(t[ii,j,kk] > -1) {
          a<-a+t[ii,j,kk]
          t[ii,j,kk]<- -1
        }
      }
    }
    pvalues[i,j,k]<-1-a
    for(ii in i:n.dmax) {
      for(kk in 1:k) {
        if(pvalues[ii,j,kk] == -1) {
          pvalues[ii,j,kk]<-1-a
        }
      }
    }
  }
  return(pvalues)
}

# Test each group.
# Performs type II test (see article).
# data: The clustered groups (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# sim: The simulation result (data.frame).
# n, n.nmin, n.dmax, n.delta: number of classes to use for discretization for each axe.
# Return a data.frame similar to data, with a new column named Pred 'level', telling if a group is significant at the given level.
test2<-function(data, sim, group.sizes, n=10, n.nmin=n, n.dmax=n, n.delta=n)
{ 
  cn<-"p.value"
  data[,cn] <- numeric(nrow(data))
  # Compute dist for each group size:
  dist<-list();
  for(i in group.sizes) {
    dist[[as.character(i)]]<-disc.dist2(sim,group.size=i,n.nmin=n.nmin,n.dmax=n.dmax,n.delta=n.delta)
  }
  
  for(i in 1:nrow(data)) {
    group<-data[i,"Size"]
    nmin <-data[i,"Nmin"]
    dmax <-data[i,"Dmax"]
    delta<-data[i,"Delta"]
    nmin.i<-which(nmin<=as.numeric(colnames(dist[[as.character(group)]])))[1]
    dmax.i<-which(dmax<=as.numeric(rownames(dist[[as.character(group)]])))[1]
    delta.i<-which(delta<=as.numeric(dimnames(dist[[as.character(group)]])[[3]]))[1]
    if(is.null(dist[[as.character(group)]])) {
      data[i,cn]<-NA
    } else {
      if(is.na(nmin.i)) { # Value is out of simulation range.
        warning(paste("Nmin value out of range for group",i,"."))
        data[i,cn]<-NA
      } else {
        data[i,cn]<-dist[[as.character(group)]][dmax.i,nmin.i,delta.i]
      }
    }
  }
  data$code<-pval(data$p.value)
  return(data)
}

# Compute how many simulated group are observed for each group size and each nmin.
# sim: The simulation result (data.frame).
# group.sizes Group sizes to test.
# n.nmin Number of classes to use for nmin discretization.
compute.eff<-function(sim, group.sizes, n.nmin=10)
{
  l<-list()
  for(i in group.sizes)
  {
    sim.i<-sim[sim$Size==i,]
    n.lim<-seq(0,max(sim.i$Nmin),len=n.nmin+1)
    nmin.f<-cut(sim.i$Nmin, breaks=n.lim)
    t<-table(nmin.f)
    dimnames(t)[[1]]<-n.lim[-1]
    l[[as.character(i)]]<-t
  }
  return(l)
}

# Display method.
# Test all groups at the 5%, 1% and 0.1% level with a type I test (see article).
# Add a 'Conf' column telling how many simulated groups have a given nmin and group size.
# Optionaly add a column specifying the type of vector used.
# data: The clustered groups (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# sim: The simulation result (data.frame).
# level: The maximum p-value level for groups to output.
# conf.level: remove groups with a conf value lower than this parameter.
# Return a data.frame similar to data, with 3 new column named Pred 0.05, Pred 0.01 and Pred 0.001.
# A 10*10 discretization grid is applied, with one tail tests.
get.pred<-function(data, sim, group.sizes, method="", level=0.05, conf.level=50)
{
  pred<-test(data, sim, group.sizes)
  conf<-compute.eff(sim, group.sizes);

  pred<-pred[!is.na(pred$p.value),]
  pred<-pred[pred$p.value <= level & pred[,"Const"] == "no",]
  if(nrow(pred)==0) return(pred)
  pred<-pred[order(pred$p.value),]
  if(method != "")
  {
    pred[,"Method"]<-rep(method, nrow(pred))
  }
  # Confidence:
  pred[,"Conf"]<-numeric(nrow(pred))
  for(i in 1:nrow(pred))
  {
    conf.i<-conf[[as.character(pred[i,"Size"])]]
    nmin.i<-which(pred[i,"Nmin"] < as.numeric(rownames(conf.i)))[1]
    pred[i,"Conf"]<-conf.i[nmin.i]
  }
  return(pred[pred$Conf>=conf.level,])
}

# Display method.
# Test all groups at the 5%, 1% and 0.1% level with a type II test (see article).
# Add a 'Conf' column telling how many simulated groups have a given nmin and group size.
# Optionaly add a column specifying the type of vector used.
# data: The clustered groups (data.frame).
# sim: The simulation result (data.frame).
# group.sizes: A vector of sizes to test. 2:20 will test all groups with a size <= 20 for instance.
# level: The maximum p-value level for groups to output.
# conf.level: remove groups with a conf value lower than this parameter.
# Return a data.frame similar to data, with 3 new column named Pred 0.05, Pred 0.01 and Pred 0.001.
# A 10*10 discretization grid is applied, with one tail tests.
get.pred2<-function(data, sim, group.sizes, method="", level=0.05, conf.level=50)
{
  pred<-test2(data, sim, group.sizes)
  conf<-compute.eff(sim, group.sizes);

  pred<-pred[!is.na(pred$p.value),]
  pred<-pred[pred$p.value <= level & pred[,"Const"] == "no",]
  if(nrow(pred)==0) return(pred)
  pred<-pred[order(pred$p.value),]
  if(method != "")
  {
    pred[,"Method"]<-rep(method, nrow(pred))
  }
  # Confidence:
  pred[,"Conf"]<-numeric(nrow(pred))
  for(i in 1:nrow(pred))
  {
    conf.i<-conf[[as.character(pred[i,"Size"])]]
    nmin.i<-which(pred[i,"Nmin"] < as.numeric(rownames(conf.i)))[1]
    pred[i,"Conf"]<-conf.i[nmin.i]
  }
  return(pred[pred$Conf>=conf.level,])
}
