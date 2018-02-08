#prep_data = function(bilan){


#   assign('dta', data.table(copy(dta)), .GlobalEnv)
#   assign('sce', data.table(copy(sce)), .GlobalEnv)
# }


run = function(scen, bilan){

  b = bil.clone(bilan)
  bil.set.values(b, scen)
  bil.pet(b)
  bil.run(b)[, -c("H", "WEI", "R", "B"), with = FALSE]

}

calc_cordex = function(bilan, per_label = 'CTRL'){
  dta = data.table(bil.get.data(bilan))[, .(DTM, month = month(DTM), P, T)]
  sce = delty[dta, on = 'month', allow = TRUE]
  dta = dta[, .(month = month(DTM), SID = 'CTRL', PER = "CTRL", EXP = 'obs', dPR = 1, dTAS = 0, DTM, P, T)]
  sce = rbind(dta, sce)

  sce[, sP := P * dPR]
  sce[, sT := T + dTAS]

  data("delty")
  res = sce[, run(data.table(DTM, P = sP, T = sT), bilan), by = .(SID, PER, EXP) ]
  copy(res)
}

calc_sens_mean = function(bilan, f = 0.05, samples = 20, sdP = stat[, mean(sdP, na.rm = TRUE)], sdT = stat[, mean(sdT, na.rm = TRUE)], per_label = "CTRL"){

  data(stat)
  data(cyc)
  data("delty")

  dta = data.table(bil.get.data(bilan))[, .(DTM, month = month(DTM), P, T)]
  sce = delty[dta, on = 'month', allow = TRUE]
  dta = dta[, .(month = month(DTM), SID = 'CTRL', PER = "CTRL", EXP = 'obs', dPR = 1, dTAS = 0, DTM, P, T)]
  sce = rbind(dta, sce)

  sce[, sP := P * dPR]
  sce[, sT := T + dTAS]

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
  res = SCE[, run(data.table(DTM, P = sP, T = sT), bilan), by = .(SID, PER, EXP, iP, iT, dX = dPR, dY = dTAS) ]
  copy(res)

}

calc_sens_sd = function(bilan, f = 0.05, samples = 20, meanP = stat[, mean(meanP, na.rm = TRUE)], meanT = stat[, mean(meanT, na.rm = TRUE)], per_label = "CTRL"){


  data(stat)
  data(cyc)
  data("delty")

  dta = data.table(bil.get.data(bilan))[, .(DTM, month = month(DTM), P, T)]
  sce = delty[dta, on = 'month', allow = TRUE]
  dta = dta[, .(month = month(DTM), SID = 'CTRL', PER = "CTRL", EXP = 'obs', dPR = 1, dTAS = 0, DTM, P, T)]
  sce = rbind(dta, sce)

  sce[, sP := P * dPR]
  sce[, sT := T + dTAS]

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
  res = SCE[, run(data.table(DTM, P = sP, T = sT), bilan), by = .(SID, PER, EXP, iP, iT, dX = sdP, dY = sdT) ]
  copy(res)

}


annual = function(DTM){
  return(rep('ANN', length(DTM)))
}

season = function(DTM){
  mon = month(DTM)
  c('DJF', 'MAM', 'MAM', 'MAM', 'JJA', 'JJA', 'JJA', 'SON', 'SON', 'SON', 'DJF', 'DJF')[mon]

}

stats = function(bilan, MEAN, SD, CORDEX, fun = mean, var = "RM", type = annual, diff_type = "multiplicative"){

  dta = data.table(bil.get.data(bilan))[, .(DTM, month = month(DTM), P, T)]
  dta = dta[, .(month = month(DTM), SID = 'CTRL', PER = "CTRL", EXP = 'obs', dPR = 1, dTAS = 0, DTM, P, T)]
  d = copy(dta[EXP=='obs'])

  mr = MEAN[, fun(eval(parse(text = var))), by = .(dX, dY, SEASON = type(DTM)) ]
  sr = SD[, fun(eval(parse(text = var))), by = .(dX, dY, SEASON = type(DTM)) ]
  cres = CORDEX[, .(fun(eval(parse(text = var)))), by = .(SID, PER, EXP, SEASON = type(DTM))]
  cres = stat[cres, on = c('SID', 'PER', 'EXP')]

  obs = cres[EXP=='obs', .(SEASON, CTRL = V1)]

  mr = obs[mr, on = 'SEASON']
  sr = obs[sr, on = 'SEASON']
  cres = obs[cres, on = 'SEASON']

  diffun = if (diff_type=="multiplicative") {function(a, b)a/b} else {function(a, b)a-b}

  mr[, V1:= diffun(V1,CTRL)]
  sr[, V1:= diffun(V1,CTRL)]
  cres[, V1:= diffun(V1,CTRL)]

  rbindlist(list(MEAN = mr, SD = sr, CORDEX = cres), fill = TRUE, idcol = "SENS")

}
