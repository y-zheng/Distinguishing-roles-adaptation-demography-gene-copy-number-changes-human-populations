mhc.evolve.yonly = function (pop, fitlist, rec, ne.n) {
 ne = length(pop$cpn)/2
 par1 = sample(1:ne,ne.n,replace=TRUE,prob=pop$fit) #Choose parent individuals (diploid)
 par2 = sample(1:ne,ne.n,replace=TRUE,prob=pop$fit) 
 
 rec.help = c((1-rec)/2,(1-rec)/2,rec/2,rec/2)
 rec.type = sample(1:4,2*ne.n,replace=TRUE,prob=rec.help) # 1-> Gamete from paternal copy, 2-> Gamete from maternal copy, 3-> Paternal first half, maternal second half, 4-> vv
 
 cpn.p = pop$cpn[c(par1,par2)]
 cpn.m = pop$cpn[c(par1+ne,par2+ne)]
 
 cpn.o = rep(0,2*ne.n)
 cpn.o[rec.type==1] = cpn.p[rec.type==1]
 cpn.o[rec.type==2] = cpn.m[rec.type==2]
 for (i in which(rec.type==3)) {
  rec.bp = c(ceiling(runif(1,0,1)*cpn.p[i]),ceiling(runif(1,0,1)*cpn.m[i]))
  cpn.o[i] = rec.bp[1]+(cpn.m[i]-rec.bp[2])
 }
 for (i in which(rec.type==4)) {
  rec.bp = c(ceiling(runif(1,0,1)*cpn.m[i]),ceiling(runif(1,0,1)*cpn.p[i]))
  cpn.o[i] = rec.bp[1]+(cpn.p[i]-rec.bp[2])
 }
  
 y = cpn.o[1:ne.n]+cpn.o[(ne.n+1):(2*ne.n)]
 fit = fitlist[y]
 if (sum(is.na(fit))>0) {fit[is.na(fit)] = 0}
 
 return(list(cpn = cpn.o, fit = fit, y = y))

}





arg = commandArgs(TRUE)
r = as.numeric(arg[1])
sx = as.numeric(arg[2])
sy = as.numeric(arg[3])
gname = arg[4]

k1 = 0.95
k2 = 1.05
fitlist = rep(0,300)
fitlist[1] = 1
fitlist[2] = 1
for (y in 3:300) {
 t1 = (1+sx)^(sum(k1^(0:(y-1))))
 t2 = (1-sy)^(sum(k2^(0:(y-3))))
 fitlist[y] = t1*t2
}


for (stp in 1:100) {
for (ne.bn in c(100,500,1000)) {

 infile = paste("Bottleneck_Initial/Initial_state_hapcpn_",gname,"_rep_",stp,".txt",sep="")
 inpop = read.table(infile,header=FALSE)
 newcpn = sample(rep(inpop[,1],inpop[,2]))
 newy = newcpn[1:10000]+newcpn[10001:20000]
 pop.start = list(cpn = newcpn, fit = fitlist[newy], y = newy)
 
 if (sum(is.na(pop.start$fit))>0) {pop.start$fit[is.na(pop.start$fit)] = 0}
 
 for (num in 1:100) {
  pop = pop.start
  ofile = paste("Bottleneck/",gname,"/bn",ne.bn,"/",gname,"_Bottleneck_",ne.bn,"_rep_",(num+100*(stp-1)),".txt",sep="")
  for (gen in 1:5000) {
   pop = mhc.evolve.yonly(pop,fitlist,r,ne.bn)
  }
  for (gen in 1:1000) {
   pop = mhc.evolve.yonly(pop,fitlist,r,10000)
   if (gen %% 50 == 0) {
    write.table(t(table(pop$y)),ofile,row.names=FALSE,quote=FALSE,append=TRUE)
   }
  }
 
 
 }


}


}
