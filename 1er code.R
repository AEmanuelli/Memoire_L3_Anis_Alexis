library(deSolve)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dS = gS*S*(1-(S+R)/K)-d*0.0133*S
    dR = gR*R*(1-(S+R)/K)
    return(list(c(dS, dR)))
  })
}
Pars1 <- c(gS = 0.0033, gR = 0.0017, d = 1, K=1)
Pars2 = c(gS = 0.0033, gR = 0.0017, d = 0, K=1)

State <- c(S = .55, R = .20)
Time = seq(0,2555, by = 1)
Time_1 <- seq(0, 7, by = 1/365)

out1 = as.data.frame(ode(func = LotVmod, y = State, parms = Pars1, times = Time))
out2 = as.data.frame(ode(func = LotVmod, y = State, parms = Pars2, times = Time))


#On a regardé "à la main" chaque fois qu'un seuil était atteint
State1 = c(S = 0.304, R = 0.206) #0.0.511
State2 = c(S = 0.487, R = 0.262) #0.75
State3 = c(S = 0.239, R = 0.272) #0.0.511
State4 = c(S = 0.397, R = 0.353) #0.75
State5 = c(S = 0.136, R = 0.374) #0.0.511
State6 = c(S = 0.215, R = 0.5) #0.75

Time <- seq(0, 2555, by = 1)
Time2 <- seq(0, 50, by = 1)
Time3 <- seq(51, 444, by = 1)
Time4 <- seq(445, 505, by = 1)
Time6 <- seq(506, 929, by = 1)
Time7 <- seq(930, 1020, by = 1)
Time8 <- seq(1021, 1508, by = 1)
TimeInfini = seq(1509, 2555, by = 1)

out1 = as.data.frame(ode(func = LotVmod, y = State, parms = Pars1, times = Time))
out2 = as.data.frame(ode(func = LotVmod, y = State, parms = Pars2, times = Time))
out8 = as.data.frame(ode(func = LotVmod, y = State, parms = Pars1, times = Time2))
out9 = as.data.frame(ode(func = LotVmod, y = State1, parms = Pars2, times = Time3))
out10 = as.data.frame(ode(func = LotVmod, y = State2, parms = Pars1, times = Time4))
out11 = as.data.frame(ode(func = LotVmod, y = State3, parms = Pars2, times = Time6))
out12 = as.data.frame(ode(func = LotVmod, y = State4, parms = Pars1, times = Time7))
out13 = as.data.frame(ode(func = LotVmod, y = State5, parms = Pars2, times = Time8))
out14 = as.data.frame(ode(func = LotVmod, y = State6, parms = Pars1, times = TimeInfini))


outinfini = c(rowSums(out8[,-1]), rowSums(out9[,-1],), 
              rowSums(out10[,-1]), rowSums(out11[,-1]), 
              rowSums(out12[,-1]), rowSums(out13[,-1]),
              rowSums(out14[,-1]))



out = data.frame(a =rowSums(out1[,-1]),B = rowSums(out2[,-1]), AD = outinfini)


#taille de la tumeur 
plot(Time/365,out$a,
     xlab="Time (Years)",ylab="Total Tumour Burden",
     ylim=c(0,1),type="l",col="blue")
lines(Time/365,out$B,col="red")
lines(Time/365,out$AD,col="green")
mtext("No therapy",col="red",at=5,padj=2.2)
mtext("Continuous",col="blue",at=4,padj=5)
mtext("Adaptive",col="green",at=6,padj=10)

Surv = function(t,l){
  exp(-.003*t*l)
}

#fonction de survie 
plot(Time/365,Surv(Time,out$a),
     xlab="Time (Years)",ylab="Survival Probability",
     ylim=c(0,1),type="l",col="blue")
lines(Time/365,Surv(Time,out$B),col="red")
lines(Time/365,Surv(Time,out$AD),col="green")
mtext("No therapy",col="red",at=5,padj=2)
mtext("Continuous",col="blue",at=5,padj=3.5)
mtext("Adaptive",col="green",at=5,padj=5)
