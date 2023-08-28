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


mhc.migration = function (pop1, pop2, nmig) {#pop1 is source, pop2 is recipient. Returns updated pop2.
 ne.1 = length(pop1$y)
 ne.2 = length(pop2$y)

 mig.ind = sample(ne.1, nmig)

 return(list(cpn = c(pop2$cpn[1:ne.2],pop1$cpn[mig.ind],pop2$cpn[ne.2+(1:ne.2)],pop1$cpn[ne.2+mig.ind]), fit = c(pop2$fit,pop1$fit[mig.ind]), y = c(pop2$y,pop1$y[mig.ind])))
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


grr = log(1450/280)/5095
ne.t = round(280*exp(grr*(0:5095)))

grr = log(19600/890)/895
ne.ceu = round(890*exp(grr*(0:895)))

grr = log(42200/560)/895
ne.chb = round(560*exp(grr*(0:895)))

ofile.ceu = paste("GADMA_Output/",gname,"_GADMA_CEU.txt",sep="")
ofile.chb = paste("GADMA_Output/",gname,"_GADMA_CHB.txt",sep="")

for (stp in 1:100) {

 infile = paste("Bottleneck_Initial/Initial_state_hapcpn_",gname,"_rep_",stp,".txt",sep="")
 inpop = read.table(infile,header=FALSE)
 newcpn = sample(rep(inpop[,1],inpop[,2]))
 newy = newcpn[1:10000]+newcpn[10001:20000]
 pop.yri = list(cpn = newcpn, fit = fitlist[newy], y = newy)
 
 if (sum(is.na(pop.yri$fit))>0) {pop.yri$fit[is.na(pop.yri$fit)] = 0}
 
 for (num in 1:100) {
  
  pop.euas = pop.yri
  for (gen in 1:5096) {
   pop.euas = mhc.evolve.yonly(pop.euas,fitlist,r,ne.t[gen])
   nmig = rpois(1,ne.t[gen]*5.74e-4)
   if (nmig > 0) {pop.euas = mhc.migration(pop.yri,pop.euas,nmig)}
   
  }
  
  pop.ceu = pop.euas
  pop.chb = pop.euas
  for (gen in 1:896) {
   pop.ceu = mhc.evolve.yonly(pop.ceu,fitlist,r,ne.ceu[gen])
   pop.chb = mhc.evolve.yonly(pop.chb,fitlist,r,ne.chb[gen])

   nmig = rpois(1,ne.ceu[gen]*2.8e-5)
   if (nmig > 0) {pop.ceu = mhc.migration(pop.yri,pop.ceu,nmig)}
   nmig = rpois(1,ne.chb[gen]*2.1e-5)
   if (nmig > 0) {pop.chb = mhc.migration(pop.yri,pop.chb,nmig)}

   nmig = rpois(1,ne.ceu[gen]*1.88e-4)
   if (nmig > 0) {pop.ceu = mhc.migration(pop.chb,pop.ceu,nmig)}
   nmig = rpois(1,ne.chb[gen]*6.8e-5)
   if (nmig > 0) {pop.chb = mhc.migration(pop.ceu,pop.chb,nmig)}


  }
 
  write.table(t(table(pop.ceu$y)),ofile.ceu,row.names=FALSE,quote=FALSE,append=TRUE)
  write.table(t(table(pop.chb$y)),ofile.chb,row.names=FALSE,quote=FALSE,append=TRUE)
 
 }




}
