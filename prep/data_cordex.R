require(mha)

setwd('/home/hanel/DATA/geo')
cr = readOGR('cr0.shp')

setwd('/media/mha/ADATA NH13/CORDEX_MONTHLY/')

getX = function(x, s, i){sapply(strsplit(x, s), function(y)y[i])}
lst = data.table(FNAM = dir())
lst[, SID:=gsub('pr_|tas_|\\.nc', '', FNAM)]
lst[, N:= .N, by = SID]
lst = lst[N==2, ]
lst[, RES:=getX(SID, '_', 1)]
lst[, GCM:=getX(SID, '_', 2)]
lst[, EXP:=getX(SID, '_', 3)]
lst[, EID:=getX(SID, '_', 4)]
lst[, RCM:=getX(SID, '_', 5)]
lst[, RID := paste(RES, GCM, EID, RCM, sep = '_')]
lst[, NN:=.N, by = RID]

lst = lst[!duplicated(SID)]
lst = lst[EXP!='historical']
lst = lst[, .(EXP = c('historical', EXP)), by = .(SID, GCM, RCM, EID)]

#lst = lst[, .(FN = paste0(c('pr_', 'tas_'), SID, '.nc')), by = .(SID, GCM, RCM, EID, EXP)]


i = 1
DTA = list()
for (i in 1:nrow(lst)){
  cat(i, '\n')
  n = lst[i, gsub(getX(SID, '_', 3), EXP, SID  )]
  p = brick(paste0('pr_', n, '.nc'))
  pr = mask(p, cr)
  t = brick(paste0('tas_', n, '.nc'))
  ta = mask(t, cr)
  
  apr = as.array(pr)
  ata = as.array(ta)
  apr = apply(apr, 3, mean, na.rm = TRUE) * 60 * 60 * 24
  ata = apply(ata, 3, mean, na.rm = TRUE) - 272.15
  
  AP = data.table(DTM = p@z$Date, PR = apr)
  AT = data.table(DTM = t@z$Date, TAS = ata)
  dta = AP[AT, on = 'DTM']  
  dta[, EXP := lst[i, EXP]]
  dta[, SID := lst[i, SID]]
  DTA[[length(DTA)+1]] = copy(dta)
}

DTA = rbindlist(DTA)

DTA[, PER := NA_character_]
DTA[year(DTM) %in% c(1970:1999), PER := 'CTRL']
DTA[year(DTM) %in% c(2021:2050), PER := '2035']
DTA[year(DTM) %in% c(2070:2099), PER := '2085']
DTA[year(DTM) %in% c(2041:2070), PER := '2055']
DTA = DTA[!is.na(PER)]

DTA[TAS< -200, TAS := TAS + 272.15]

sta = DTA[, .(PR = mean(PR), TAS = mean(TAS)), by = .(EXP, SID, PER, month(DTM))]
ctrl = sta[PER=='CTRL']
setnames(ctrl, c('PR', 'TAS'), c('cPR', 'cTAS'))

sta = ctrl[, .(month, SID, cPR, cTAS)][sta, on = c('SID', 'month')]

sta = sta[, .(dPR = PR/cPR, dTAS = TAS - cTAS), by = .(month, SID, PER, EXP)][EXP!='historical']

setwd('/home/hanel/LAPV/soft/')
saveRDS(sta, 'delty.rds')
