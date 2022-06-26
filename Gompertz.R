library(deSolve)
a = .000002
b= 7


Gomp <- function (t, x, parms) {
  with(as.list(c(x, parms)), {
    dS = (p*log(K/(S+R))*(1-8*d)*S)
    dR = p*log(K/(S+R))*R
    #dP = -0.002*((S+R)^5)*P
    dP = -a*exp(b*(S+R))*P
    return(list(c(dS, dR, dP)))
  })
}

Time = seq(0,7*365,by=1)
avec_G= c(p = .002335398, K=1, d = 1, X = 0)
sans_G <- c(p = .002335398, K=1, d = 0, X = 0)
State_0_G=c(S=.75*.73, R = .75*.27 ,P=1)

out_c_G = data.frame(ode(y = State_0_G,times = Time,func = Gomp,parms = avec_G))
out_nt_G = data.frame(ode(y = State_0_G,times = Time,func = Gomp,parms = sans_G))


out1 = out_c_G
m = .75*.7
M=.75
t1
#adaptative = function(out1, m, M, func, pav, psa){
  t1=which((out1$S+out1$R)<m)[1]
  out_ad=out1[out1$time<out1$time[t1],]
  State=c(S=out1$S[t1], R = out1$R[t1], P = out1$P[t1])
  out = data.frame(ode(y = State,times = Time, func = Gomp, parms = sans_G))
  
  t1=which((out$S+out$R)>M)[1]
  f1=0
  f2=0
  
  while (f1==0){
    if (is.na(t1)){
      out$time=out$time+max(out_ad$time)
      out_ad=rbind(out_ad,out)
      f1=1}
    
    else{
      apartde=out[out$time<out$time[t1],]
      apartde$time=apartde$time+max(out_ad$time)
      out_ad=rbind(out_ad,apartde)
      if (f2==0){
        State =c(S=out$S[t1], R = out$R[t1], P = out$P[t1])
        out = data.frame(ode(y = State,times = Time,func = Gomp,parms = avec_G))
        f2=1
        t1=which((out$S+out$R)<m)[1]
      }
      
      else {
        State = c(S=out$S[t1], R = out$R[t1], P = out$P[t1])
        out = data.frame(ode(y = State,times = Time,func = Gomp,parms = sans_G))
        f2=0
        t1=which((out$S+out$R)>M)[1]
      }
    }
  }
  
  #return (out_ad)
#}
out_ad_G=out_ad


#out_ad_G = adaptative(out1 = out_c_G, m= 75*.7, M= .75, func = Gomp, pav = avec_G, psa = sans_G)


plot(out_nt_G$time/365,out_nt_G$S+out_nt_G$R,ylim =c(0,1),
     xlab="Time (Years)",ylab="Total Tumour Burden",type="l",col="red")
lines(out_c_G$time/365,out_c_G$S+out_c_G$R,col="blue")
lines(out_ad_G$time/365,out_ad_G$S+out_ad_G$R,col="green")
mtext("No therapy",col="red",at=5,padj=2.2)
mtext("continuous",col="blue",at=3,padj=10)
mtext("Adaptive",col="green",at=6,padj=6)
abline(h = c(.31, .22), v = 1)


surv1 = function(t,l) {
  exp(-0.003*t*l)
}





par(mfrow=c(1,2))
#plot tumor growth
plot(out_nt_G$time/365,out_nt_G$S+out_nt_G$R,
     xlab="Time (Years)",ylab="Total Tumour Burden",
     ylim=c(0,1),type="l",col="red")
lines(out_c_G$time/365,out_c_G$S+out_c_G$R,col="blue")
lines(out_ad_G$time/365,out_ad_G$S+out_ad_G$R,col="green")
mtext("No therapy",col="red",at=5.5,padj=8)
mtext("Continuous",col="blue",at=5.5,padj=9)
mtext("Adaptive",col="green",at=5.5,padj=10)
#abline(h = .75*.8, v=56/365)

plot(out_nt_G$time/365,out_nt_G$P,ylim =c(0,1),
     xlab="Time (Years)",ylab="Survival Probability",type="l",col="red")
lines(out_c_G$time/365, out_c_G$P,col="blue")
#lines(out_ad_G$time/365,surv1(out_ad_G$time,(out_ad_G$R+out_ad_G$S)) ,col="green")
lines(out_ad_G$time/365,out_ad_G$P,col="green")
#lines(Time/365,surv1(Time,(out_gp$t2)) ,col="brown")
mtext("No therapy",col="red",at=5,padj=2)
mtext("Continuous",col="blue",at=5,padj=3.5)
#mtext("Adaptive",col="green",at=5,padj=5)
mtext("Adaptative",col="green",at=5,padj=5)
abline(h = c(.26), v =5)