amerika = bilan::bil.new(type = 'm', file = 'VD_Amerika.bil')
dta = bilan::bil.get.data(amerika)
dta = dta[data.table::year(dta$DTM) %in% c(1970:1999),]
bilan::bil.set.values(amerika, dta)
bilan::bil.pet(amerika)
amerika
