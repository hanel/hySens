run = function(scen, bilan){

  b = bil.clone(bilan)
  bil.set.values(b, scen)
  bil.pet(b)
  bil.run(b)[, -c("H", "WEI", "R", "B"), with = FALSE]

}

wat_res = function(DTM, Q, S = NULL, SA = 0, R = NULL, Y = NULL, P = NULL, EV = NULL, ALT = NULL, EAS = NULL, WU = NULL){
  #S - zasoba [m3], SA- plocha zatopy[m2], R-zabezpecenost, Y-nadlepseny odtok, P-srazky, EV-rada evaporace, ALT-nadmorska vyska pro vypocet evaporace z nomogramu z normy
  #EAS je elvation-area-storage relationship, ale to asi vetsinou nebude k dispozici, WU - odber pripadne pritok do nadrze, dle znamenka
  #pro optimalizaci objemu je zapotrebi dodat nejaky pocatecni objem a zatopenou plochu
  if (is.null(S)) reser = as.wateres(data.frame(DTM, Q), storage = Y*365*24*3600, area = SA, eas = EAS) else reser = as.wateres(data.frame(DTM, Q), storage = S, area = SA, eas = EAS)
  if (!is.null(P)) reser = set_precipitation(reser, P)
  if (!is.null(EV)) reser = set_evaporation(reser, EV)
  if (!is.null(ALT)) reser = set_evaporation(reser, altitude = ALT)
  if (!is.null(WU)) reser = set_wateruse(reser, values = WU)
  #vlastni vypocet storage-reliability-yield relationship
  if (is.null(S)) res = summary(reser, reliability = R, yield = Y, get_series = F, upper_limit = 100)
  if (is.null(Y)) res = summary(reser, storage = S, reliability = R, get_series = F, upper_limit = 100)
  if (is.null(R)) res = summary(reser, storage = S, yield = Y, get_series = F)
  res
}


calc_cordex = function(bilan, per_label = 'CTRL', CA=NULL, S = NULL, SA = 0, R = NULL, Y = NULL, EV = NULL, ALT = NULL, EAS = NULL, WU = NULL){ #doplnene promenne pro wateres
  #CA - plocha povodi [km2]
  dta = data.table(bil.get.data(bilan))[, .(DTM, month = month(DTM), P, T)]
  sce = delty[dta, on = 'month', allow = TRUE]
  dta = dta[, .(month = month(DTM), SID = 'CTRL', PER = "CTRL", EXP = 'obs', dPR = 1, dTAS = 0, DTM, P, T)]
  sce = rbind(dta, sce)

  sce[, sP := P * dPR]
  sce[, sT := T + dTAS]

  data("delty")
  res = sce[, run(data.table(DTM, P = sP, T = sT), bilan), by = .(SID, PER, EXP) ]

  if (!is.null(CA)) {
    res[, Q := (RM*CA)/(3.6*24*30.5)] #pridani prutoku (prepocet z mm na m3/s)
    res = res[, wat_res(DTM, Q, S = S, SA = SA, R = R, Y = Y, P, EAS = EAS, EV = EV, ALT = ALT, WU = WU),  by = .(SID, PER, EXP)]
  }

  copy(res)
}


calc_sens_mean = function(bilan, f = 0.05, samples = 20, sdP = stat[, mean(sdP, na.rm = TRUE)], sdT = stat[, mean(sdT, na.rm = TRUE)], CA=NULL, S = NULL, SA = 0, R = NULL, Y = NULL, EV = NULL, ALT = NULL, EAS = NULL, WU = NULL, per_label = "CTRL"){

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

  if (!is.null(CA)) {
    res[, Q := (RM*CA)/(3.6*24*30.5)] #pridani prutoku (prepocet z mm na m3/s)
    res = res[, wat_res(DTM, Q, S = S, SA = SA, R = R, Y = Y, P, EAS = EAS, EV = EV, ALT = ALT, WU = WU),  by = .(SID, PER, EXP, iP, iT, dX, dY)]
  }

  copy(res)

}


calc_sens_sd = function(bilan, f = 0.05, samples = 20, meanP = stat[, mean(meanP, na.rm = TRUE)], meanT = stat[, mean(meanT, na.rm = TRUE)], CA=NULL, S = NULL, SA = 0, R = NULL, Y = NULL, EV = NULL, ALT = NULL, EAS = NULL, WU = NULL, per_label = "CTRL"){


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

  if (!is.null(CA)) {
    res[, Q := (RM*CA)/(3.6*24*30.5)] #pridani prutoku (prepocet z mm na m3/s)
    res = res[, wat_res(DTM, Q, S = S, SA = SA, R = R, Y = Y, P, EAS = EAS, EV = EV, ALT = ALT, WU = WU),  by = .(SID, PER, EXP, iP, iT, dX, dY)]
  }

  copy(res)

}


annual = function(DTM){
  return(rep('ANN', length(DTM)))
}

season = function(DTM){
  mon = month(DTM)
  c('DJF', 'MAM', 'MAM', 'MAM', 'JJA', 'JJA', 'JJA', 'SON', 'SON', 'SON', 'DJF', 'DJF')[mon]

}

stats = function(MEAN, SD, CORDEX, fun = mean, var = 'RM', type = annual, diff_type = 'multiplicative'){
  #myslim, ze tohle tu neni potreba ?
  #dta = data.table(bil.get.data(bilan))[, .(DTM, month = month(DTM), P, T)]
  #dta = dta[, .(month = month(DTM), SID = 'CTRL', PER = "CTRL", EXP = 'obs', dPR = 1, dTAS = 0, DTM, P, T)]
  #d = copy(dta[EXP=='obs'])

  if('DTM' %in% names(MEAN)){

    mr = MEAN[, fun(eval(parse(text = var))), by = .(dX, dY, SEASON = type(DTM)) ]
    sr = SD[, fun(eval(parse(text = var))), by = .(dX, dY, SEASON = type(DTM)) ]
    cres = CORDEX[, .(fun(eval(parse(text = var)))), by = .(SID, PER, EXP, SEASON = type(DTM))]
    cres = stat[cres, on = c('SID', 'PER', 'EXP')]
    obs = cres[EXP=='obs', .(SEASON, CTRL = V1)]
    mr = obs[mr, on = 'SEASON']
    sr = obs[sr, on = 'SEASON']
    cres = obs[cres, on = 'SEASON']

  } else {

    mr = copy(MEAN)
    setnames(mr, eval(var), 'V1')
    sr = copy(SD)
    setnames(sr, eval(var), 'V1')
    cres = copy(CORDEX)
    cres = stat[cres, on = c('SID', 'PER', 'EXP')]
    setnames(cres, eval(var), 'V1')
    obs = cres[EXP=='obs', .(CTRL = V1)]
    mr[, c('CTRL','SEASON'):=.(as.numeric(obs), 'ANN')]
    sr[, c('CTRL','SEASON'):=.(as.numeric(obs), 'ANN')]
    cres[, c('CTRL','SEASON'):=.(as.numeric(obs), 'ANN')]

  }

  diffun = if (diff_type=="multiplicative") {function(a, b)a/b} else {function(a, b)a-b}

  mr[, V1:= diffun(V1, CTRL)]
  sr[, V1:= diffun(V1, CTRL)]
  cres[, V1:= diffun(V1, CTRL)]

  resul = rbindlist(list(MEAN = mr, SD = sr, CORDEX = cres), fill = TRUE, idcol = "SENS")
  resul
}
