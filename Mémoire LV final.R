library(deSolve)

a1 = .007
alpha= 5
a = .000009456
b = 7

LV = function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dS = gS*S*(1-(S+R)/K) - d*S*.0133
    dR = gR*R*(1-(R+S)/K)
    dP = -a*exp(b*(S+R))*P
    #dP = -0.003*((S+R)^2)*P
    list(c(dS, dR, dP))
  })
}
avec=c(gS = 0.0033, gR = 0.0017,K=1, a2 = 0.003, d=1)
sans=c(gS = 0.0033, gR = 0.0017,K=1, a2 = 0.003, d=0)

State_0=c(S=.55, R = .2, P = 1)
State_appres = c(S=.74, R =.01 , P=1)

Time = seq(0,7*365,by=1)

out_nt = data.frame(ode(y = State_0,times = Time,func = LV,
                        parms = sans))
out_c = data.frame(ode(y = State_0,times = Time,func = LV,
                       parms = avec))



out_nt_ar = data.frame(ode(y = State_appres,times = Time,func = LV,
                           parms = sans))
out_c_ar = data.frame(ode(y = State_appres,times = Time,func = LV,
                          parms = avec))

adaptative = function(out_c, m, M){
  t1=which((out_c$S+out_c$R)<m)[1]
  out_ad=out_c[out_c$time<out_c$time[t1],]
  State=c(S=out_c$S[t1], R = out_c$R[t1], P = out_c$P[t1])
  out = data.frame(ode(y = State,times = Time,func = LV,
                       parms = sans))
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
        out = data.frame(ode(y = State,times = Time,func = LV,
                             parms = avec))
        f2=1
        t1=which((out$S+out$R)<m)[1]
      }
      else {
        State = c(S=out$S[t1], R = out$R[t1], P = out$P[t1])
        out = data.frame(ode(y = State,times = Time,func = LV,
                             parms = sans))
        f2=0
        t1=which((out$S+out$R)>M)[1]
      }
    }
  }
  return (out_ad)
}

out_ad=adaptative(out_c,  M= .75, m = .75*.7)
#out_ad_ar=adaptative(out_c_ar,  M= 0.75, m = .75*.7)



par(mfrow=c(1,2))
#plot tumor growth
plot(out_nt$time/365,out_nt$S+out_nt$R,
     xlab="Time (Years)",ylab="Total Tumour Burden",
     ylim=c(0,1),type="l",col="red")
lines(out_c$time/365,out_c$S+out_c$R,col="blue")
lines(out_ad$time/365,out_ad$S+out_ad$R,col="green")
mtext("No therapy",col="red",at=6,padj=8)
mtext("Continuous",col="blue",at=6,padj=9)
mtext("Adaptive",col="green",at=6,padj=10)
#abline(h=.6, v= .1538)

# plot survival
plot(out_nt$time/365, out_nt$P, type="l",col="red",
     xlab="Time (Years)", ylab="Survival Probability")
lines(out_c$time/365,out_c$P,col="blue")
lines(out_ad$time/365,out_ad$P,col="green")
mtext("No therapy",col="red",at=5.5,padj=2)
mtext("Continuous",col="blue",at=5.5,padj=3.5)
mtext("Adaptive",col="green",at=5.5,padj=5)
abline(h = c(.26), v =c( 5))









#plot growth avec apparition de la résistnace

#plot(out_nt_ar$time/365,out_nt_ar$S+out_nt_ar$R,
#     xlab="Time (Years)",ylab="Total Tumour Burden",
#     ylim=c(0,1),type="l",col="red")
#lines(out_c_ar$time/365,out_c_ar$S+out_c_ar$R,col="blue")
#lines(out_ad_ar$time/365,out_ad_ar$S+out_ad_ar$R,col="green")
#mtext("No therapy",col="red",at=5,padj=2.2)
#mtext("continuous",col="blue",at=4,padj=5)
#mtext("Adaptive",col="green",at=6,padj=10)

#surv1 = function(t,l) {
#  exp(-0.003*t*l)
#}


#lines(out_ad$time/365, surv1(out_ad$time, (out_ad$R+out_ad$S)),col="purple")
#mtext("Adaptive calcul direct",col="purple",at=5,padj=6)

#lines(out_c$time/365, surv1(out_c$time, (out_c$R+out_c$S)),col="pink")
#lines(out_nt$time/365,surv1(out_nt$time,out_gp$t1) ,col="orange")


#lines(out_ad$time/365, surv2((out_ad$R+out_ad$S)),col="orange")

#mtext("Adaptive nouveau",col="orange",at=5,padj=7)

