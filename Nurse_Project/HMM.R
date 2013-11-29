library(RHmm)

levels(metabo_aggre$SW_Nr)

obs = list(); i = 1
for(p in levels(metabo_aggre$SW_Nr)){
  tmp = subset(metabo_aggre, subset=SW_Nr==p&shift=="Tagschicht")
  obs[[i]] = tmp$Gln
  i = i+1
}
names(obs)=levels(metabo_aggre$SW_Nr)


obs = sapply(obs, function(x) {sapply(x, function(y) {y[which(is.na(y))] = mean(y, na.rm = T); return(y)})})
participants = levels(metabo_aggre$SW_Nr)[sapply(obs, function(x)  !is.null(dim(x)))]
obs = obs[sapply(obs, function(x)  !is.null(dim(x)))]

## separate the cases and controls

Res_obs_control = HMMFit(obs[control], dis="NORMAL", nStates=3)
Res_obs_case = HMMFit(obs[case], dis="NORMAL", nStates=3)

state_seq = viterbi(Res_obs, obs[control])

require(depmixS4)