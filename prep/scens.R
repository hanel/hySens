library(hySens)
data(amerika)
data(delty)
data(stat)

bilan = amerika

#scen = sce[EXP=='obs', .(DTM, P = sP, T = sT)]

dta = data.table(bil.get.data(bilan))[, .(DTM, month = month(DTM), P, T)]
sce = delty[dta, on = 'month', allow = TRUE]
dta = dta[, .(month = month(DTM), SID = 'CTRL', PER = per_label, EXP = 'obs', dPR = 1, dTAS = 0, DTM, P, T)]
sce = rbind(dta, sce)

sce[, sP := P * dPR]
sce[, sT := T + dTAS]



run = function(scen, bilan){

  bil.set.values(bilan, scen)
  bil.pet(bilan)
  bil.run(bilan)[, -c("H", "WEI", "R", "B"), with = FALSE]

  }

calc_cordex = function(bilan, per_label = 'CTRL'){

  data("delty")
  res = sce[, run(data.table(DTM, P = sP, T = sT), bilan), by = .(SID, PER, EXP) ]
  copy(res)
}

calc_sens_mean = function(bilan, f = 0.05, samples = 20, sdP = stat[, mean(sdP, na.rm = TRUE)], sdT = stat[, mean(sdT, na.rm = TRUE)], per_label = "CTRL"){

  data(stat)
  data(cyc)
  rng = stat[, .(var = c('mn', 'mx'), P = extendrange(meanP, f = f), T = extendrange(meanT, f = f))]
  dp = rng[, seq(P[1], P[2], length = samples)]
  dt = rng[, seq(T[1], T[2], length = samples)]
  gr = data.table(expand.grid(dP = dp, dT = dt))
  gr[, iT := rep(1:samples, each = samples)]
  gr[, iP := rep(1:samples)]
  d = copy(dta[EXP=='obs'])
  d = cyc[d, on = 'month']

  i = 1
  SCE = list()
  for (i in 1:nrow(gr)){
    d[, c('dPR', 'dTAS', 'iP', 'iT') := .(gr[i, dP], gr[i, dT], gr[i, iP], gr[i, iT])]
    d[, c('dP', 'dT') := .((mP * sdP) + dPR, (mT * sdT) + dTAS)]
    d[, c('sP', 'sT') := .(P * dP, T + dT)]
    SCE[[i]] = copy(d)
  }

  SCE = rbindlist(SCE)
  res = SCE[, run(data.table(DTM, P = sP, T = sT), bilan), by = .(SID, PER, EXP, iP, iT) ]
  copy(res)

}

calc_sens_sd = function(bilan, f = 0.05, samples = 20, meanP = stat[, mean(meanP, na.rm = TRUE)], meanT = stat[, mean(meanT, na.rm = TRUE)], per_label = "CTRL"){

  data(stat)
  data(cyc)
  rng = stat[, .(var = c('mn', 'mx'), P = extendrange(sdP, f = f), T = extendrange(sdT, f = f))]
  dp = rng[, seq(P[1], P[2], length = samples)]
  dt = rng[, seq(T[1], T[2], length = samples)]
  gr = data.table(expand.grid(dP = dp, dT = dt))
  gr[, iT := rep(1:samples, each = samples)]
  gr[, iP := rep(1:samples)]
  d = copy(dta[EXP=='obs'])
  d = cyc[d, on = 'month']

  i = 1
  SCE = list()
  for (i in 1:nrow(gr)){
    d[, c('sdP', 'sdT', 'iP', 'iT') := .(gr[i, dP], gr[i, dT], gr[i, iP], gr[i, iT])]
    d[, c('dP', 'dT') := .((mP * sdP) + meanP, (mT * sdT) + meanT)]
    d[, c('sP', 'sT') := .(P * dP, T + dT)]
    SCE[[i]] = copy(d)
  }

  SCE = rbindlist(SCE)
  res = SCE[, run(data.table(DTM, P = sP, T = sT), bilan), by = .(SID, PER, EXP, iP, iT) ]
  copy(res)

}



require(ggplot2)

delty[, nP := (dPR-mean(dPR))/sd(dPR), by = .(SID, PER, EXP)]
delty[, nT := (dTAS-mean(dTAS))/sd(dTAS), by = .(SID, PER, EXP)]

ggplot(delty) + geom_line(aes(x=month, y = nP, group = paste(SID, PER, EXP)), alpha = .3, size = 2) #+ facet_grid(EXP~PER)
ggplot(delty) + geom_boxplot(aes(x=factor(month), y = nP))
ggplot(delty) + geom_line(aes(x=month, y = dTAS, group = SID)) + facet_grid(EXP~PER)
ggplot(delty) + geom_boxplot(aes(x=factor(month), y = nT))


cyc = delty[, .(mP = median(nP, na.rm  = TRUE), mT = median(nT)), by = month]

stat = delty[, .(meanP = mean(dPR), sdP = sd(dPR), meanT = mean(dTAS), sdT = sd(dTAS)), by = .(SID, PER, EXP)]

mstat = melt(stat, id = 1:3, measure.vars = patterns("mean", "sd"), value.name = c('mean', 'sd'))
mstat[, variable := factor(variable, labels = c('P', 'T'))]

vstat = melt(stat, id = 1:3, measure.vars = list(c(4:5), c(6:7)), value.name = c('P', 'T'))
vstat[, variable := factor(variable, labels = c('mean', 'sd'))]

ggplot(mstat) + geom_point(aes(x = mean, y = sd, col = EXP)) + facet_wrap(~variable, scale = 'free', ncol = 2)

ggplot(vstat) + geom_point(aes(x = P, y = T, shape = EXP, col = PER), size = 3, alpha = .5) + facet_wrap(~variable, scale = 'free', ncol = 2)


setwd('/home/hanel/r-packages/hySens/data/')
save(stat, file = 'stat.rda')
save(cyc, file = 'cyc.rda')

require(hySens)
data("amerika")
bilan = bil.clone(amerika)

MEAN = calc_sens_mean(bilan)
SD = calc_sens_sd(bilan)
CORDEX = calc_cordex(bilan)

annual = function(DTM){
  return(rep('ANN', length(DTM)))
}

season = function(DTM){
  mon = month(DTM)
  c('DJF', 'MAM', 'MAM', 'MAM', 'JJA', 'JJA', 'JJA', 'SON', 'SON', 'SON', 'DJF', 'DJF')[mon]

}

stats = function(bilan, MEAN, SD, CORDEX, fun = mean, var = "RM", type = annual){

  dta = data.table(bil.get.data(bilan))[, .(DTM, month = month(DTM), P, T)]
  dta = dta[, .(month = month(DTM), SID = 'CTRL', PER = "CTRL", EXP = 'obs', dPR = 1, dTAS = 0, DTM, P, T)]
  d = copy(dta[EXP=='obs'])

  mr = MEAN[, fun(eval(parse(text = var))), by = .(dX, dY, SEASON = type(DTM)) ]
  sr = SD[, fun(eval(parse(text = var))), by = .(dX, dY, SEASON = type(DTM)) ]
  cres = CORDEX[, .(fun(eval(parse(text = var)))), by = .(SID, PER, EXP, SEASON = type(DTM))]
  cres = stat[cres, on = c('SID', 'PER', 'EXP')]

  mr[, V1:= V1/cres$V1[1]]
  sr[, V1:= V1/cres$V1[1]]
  cres[, V1:= V1/V1[1]]

  rbindlist(list(MEAN = mr, SD = sr, CORDEX = cres), fill = TRUE, idcol = "SENS")

}
