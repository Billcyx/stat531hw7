library(foreach)
library(doParallel)
library(doRNG)
library(pomp)
library(tidyverse)

registerDoParallel(36)  # change this number locally, you don't have 36 cores.
registerDoRNG(123456)

junk <- foreach(
  i=1:10,.combine=c,.packages=c("pomp")
) %dopar% {
  Sys.sleep(2)  # This doesn't do anything, just "sleeps" for 2 seconds.
}

read_csv(paste0("https://kingaa.github.io/sbied/stochsim/",
                "Measles_Consett_1948.csv")) %>%
  select(week,reports=cases) -> meas

seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*E/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")


seir_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 1;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

seir_dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
  ")

seir_rmeas <- Csnippet("
  reports = rbinom(H,rho);
  ")

meas |>
  pomp(times="week",t0=0,
       rprocess=euler(seir_step,delta.t=1/7),
       dmeasure=seir_dmeas,
       rinit=seir_init,
       paramnames=c("N","Beta","mu_EI","mu_IR","rho","eta"),
       statenames=c("S","E","I","R","H"),
       params=c(Beta=15,mu_IR=0.5,rho=0.5,mu_EI=0.8,k=10,eta=0.06,N=38000)
  ) -> measSIR


# check that the basic particle filter is working
measSIR |>
  pfilter(Np=1000) -> pf


fixed_params <- c(N=38000, mu_IR=2, k=10, mu_EI = 0.8)



tic <- Sys.time()
foreach(i=1:10,.combine=c,
        .packages=c("pomp")
) %dopar% {
  measSIR |> pfilter(Np=5000)
} -> pf

pf |> logLik() |> logmeanexp(se=TRUE) -> L_pf
L_pf
toc <- Sys.time()

pf[[1]] |> coef() |> bind_rows() |>
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) |>
  write_csv("measles_params.csv")


#############################################################################local search
bake(file="local_search.rds",{
  foreach(i=1:20,.combine=c,
          .packages=c("pomp")
  ) %dopar% {
    measSIR |>
      mif2(
        Np=2000, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw_sd(Beta=0.02, rho=0.02, eta=ivp(0.02)),
        partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
        paramnames=c("Beta","rho","eta")
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- nbrOfWorkers()
  mifs_local
}) -> mifs_local
t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")

mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")



bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.combine=rbind,
          .packages=c("pomp")
  ) %dopar% {
    evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- nbrOfWorkers()
  results
}) -> results

t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")

pairs(~loglik+Beta+eta+rho,data=results,pch=16)

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


#############################################################################global search 
if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
}

set.seed(2062379496)

runif_design(
  lower=c(Beta=5,rho=0.2,eta=0),
  upper=c(Beta=80,rho=0.9,eta=1),
  nseq=400
) -> guesses

mf1 <- mifs_local[[1]]


bake(file="global_search.rds",
     dependson=guesses,{
       foreach(guess=iter(guesses,"row"), .combine=rbind,
               .packages=c("pomp")
       ) %dopar% {
         mf1 |>
           mif2(params=c(guess,fixed_params)) |>
           mif2(Nmif=100) -> mf
         replicate(
           10,
           mf |> pfilter(Np=1) |> logLik()
         ) |>
           logmeanexp(se=TRUE) -> ll
         mf |> coef() |> bind_rows() |>
           bind_cols(loglik=ll[1],loglik.se=ll[2])
       } -> results
       attr(results,"ncpu") <- nbrOfWorkers()
       results
     }) |>
  filter(is.finite(loglik)) -> results

t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-50) |>
  bind_rows(guesses) |>
  mutate(type=if_else(is.na(loglik),"guess","result")) |>
  arrange(type) -> all

pairs(~loglik+Beta+eta+rho, data=all, pch=16, cex=0.3,
      col=ifelse(all$type=="guess",grey(0.5),"red"))


################################################################global search ends
