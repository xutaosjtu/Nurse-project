library(RHmm)

levels(metabo_aggre$SW_Nr)

obs = list(); i = 1
for(p in levels(metabo_aggre$SW_Nr)){
  tmp = subset(metabo_aggre, subset=SW_Nr==p&shift=="Tagschicht")
  obs[[i]] = tmp[,valid_measures[20:32]]
  i = i+1
}
names(obs)=levels(metabo_aggre$SW_Nr)


obs = sapply(obs, function(x) {sapply(x, function(y) {y[which(is.na(y))] = mean(y, na.rm = T); return(y)})})

Res_obs = HMMFit(obs, dis="NORMAL", nStates=2)
summary(Res_obs)


state_seq = viterbi(Res_obs, obs)


