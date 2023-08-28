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
ne = as.numeric(arg[1])
r = as.numeric(arg[2])
ratio_rs = as.numeric(arg[3])
ratio_ss = as.numeric(arg[4])

s1 = r/ratio_rs
s2 = s1/ratio_ss

ofile = paste("Y_gridtest_Ne_",ne,"_r_",r,"_ratio1_",ratio_rs,"_ratio2_",ratio_ss,".txt",sep="")

k1 = 0.95
k2 = 1.05
fitlist = rep(0,100)
fitlist[1] = 1
fitlist[2] = 1
for (y in 3:100) {
 t1 = (1+s1)^(sum(k1^(0:(y-1))))
 t2 = (1-s2)^(sum(k2^(0:(y-3))))
 fitlist[y] = t1*t2
}

pop = list(cpn = rep(5,2*ne), fit = rep(fitlist[5],ne), y = rep(10,ne))


for (gen in 1:5000000) {
 pop = mhc.evolve.yonly(pop,fitlist,r,ne)
 if (gen %% ne == 0) {
  outvec = c(gen,mean(pop$y),sd(pop$y),median(pop$y))
  write(outvec,ofile,ncolumns=4,append=TRUE)
 }
}
