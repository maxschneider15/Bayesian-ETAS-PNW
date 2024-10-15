# Bayesian-model-runs-paper.R
# Author: Max Schneider
# Date Created: Apr 2024

setwd("~/Downloads/Takeout 3/Drive/M9/")
# setwd("/Users/maxs/Downloads/M9")
library(etasFLP)
library(Rcpp)
library(scales)
library(ggplot2)
library(bayesianETAS)
library(cowplot)
library(gridExtra)
library(shotGroups)
library(SpatialEpi)
library(xtable)
data(world2HiresMapEnv)
w2hr <- map_data("world2Hires")
can.df <- data.frame("long"=(w2hr$long-360), "lat"=w2hr$lat, "group"=w2hr$group)
lat <- c(47.6062, 47.0379, 48.4284, 46.8523, 47.6588, 44.0582, 45.5051, 46.2087, 42.2249)
lon <- c(-122.3321, -122.9007, -123.3656, -121.7603, -117.4260, -121.3133, -122.6750, -119.1199, -121.7817)
city <- c("Seattle", "Olypmia", "Victoria", "Mt. Rainier", "Spokane", "Bend", "Portland", "Kennewick", "Klamath Falls")
cities_map <- data.frame(city, lon, lat, stringsAsFactors = FALSE)
states <- map_data("state")
usa <- map_data("usa")
canada <- map_data("worldHires", "Canada")
# setwd("C:/Users/Hank/Dropbox/PNW-Data-Processing/D. identify-swarms/Swarm Analysis/")
source("src/functions/functions-bayesian-ETAS.R")
# source("src/functions/functions-bETAS-summaries.R")
# source("src/functions/functions-to-calculate-lambda.R")
# source("src/functions/functions-to-simulate-ETAS.R")
source("src/functions/functions-spatiotemporal-ETAS.R")
source("src/functions/functions-MLE-ETAS.R")
# source("src/modeling-files/bayesianETAS/maxLikelihoodETAS.R")
# source("src/modeling-files/bayesianETAS/my-funcs/sampleETASPosterior-MS.R")
# sourceCpp("src/modeling-files/b-ETAS/bayesianETAS/my-funcs/bETAS-MS.cpp")
# sourceCpp("src/modeling-files/b-ETAS/bayesianETAS/my-funcs/st-bETAS-MS-v2.cpp")
sourceCpp("src/modeling-files/b-ETAS/bayesianETAS/my-funcs/st-bETAS-MS-v2-20210815.cpp")
sourceCpp("src/modeling-files/b-ETAS/bayesianETAS/my-funcs/st-bETAS-MS-priors.cpp")

intl <- read.csv("data/processed-data/pnw-intl-ETAS-slab-lowM.csv")
intl <- intl[order(intl$time),]
intl.aux <- intl %>% filter(lat >= 41, lat <= 50, lon >= -126, lon <= -115.5, mag >= 2, is.na(remove),
                            potential_duplicates_20210412 == F)
intl.t <- intl %>% filter(lat >= 42, lat <= 49, lon >= -125, lon <= -116.5, mag >= 2, is.na(remove),
                          potential_duplicates_20210412 == F)
intl.t.area <- (range(intl.t$x)[2] - range(intl.t$x)[1]) * (range(intl.t$y)[2] - range(intl.t$y)[1])

intl.Nt.box <- latlong2grid(rbind(c(-125, 49), c(-116.5, 45)))
intl.Nt.box$x <- intl.Nt.box$x - min(intl.Nt.box$x)
intl.Nt.box$y <- intl.Nt.box$y - min(intl.Nt.box$y)
intl.Nt.area <- intl.Nt.box$x[2]*intl.Nt.box$y[1]

maxT.pnw <- round(as.numeric(max(intl$time)-min(intl$time) + 1), 0)
maxT.pnw.target <- round(as.numeric(max(intl.t$time)-min(intl.t$time) + 1), 0)
cat2D.true <- c(0.1, 0.02, 1.7, 0.05, 1.08, 0.1, 1.5)
# cat2D.Gevals <- read.csv("data/bayesianETAS-output/2020-experiments/st-model/cat2/Gevals-tightlong-cat2D.csv")
# cat2D.Gevals.mat <- as.matrix(cat2D.Gevals)
# cat2D.Gevals.mat <- unname(cat2D.Gevals.mat)

cat3D.true <- c(0.1, 0.006, 2.2999, 0.05, 1.08, 0.1, 1.5)
numMCMCSamples <- 500 # 500=What Gordon uses
M0 <- 2.0; approx=FALSE; this.M0 <- 2.0
burnin <- 500
sims <- sims+burnin
cat2inits <- c(0.1, 0.02, 2.2999, 0.05, 1.08, 0.32, 1.5)
cat3inits <- c(0.1, 0.006, 2.2999, 0.05, 1.08, 1, 2)
bad.inits3 <- c(0.1, 0.00002, 2.2999, 0.005, 1.5, 0.8, 2.01)
crazy.inits1 <- c(10, 0.00002, 2.2999, 0.005, 1.5, 0.32, 1.5)
cat2inits.ss <- cat2inits
cat2inits.ss[3] <- 0.935*(log(10))
cat3inits[3] <- 0.935*(log(10))
bad.inits3[3] <- 0.935*(log(10))
crazy.inits1[3] <- 0.935*(log(10))

pnw.Gevals.long <- read.csv("data/bayesianETAS-output/2020-experiments/st-model/pnw/Gevals-long-pnw-target.csv")
pnw.Gevals.mat <- as.matrix(pnw.Gevals.long)
pnw.Gevals.mat <- unname(pnw.Gevals.mat)

pnw.Nt.cat.noMCswarms <- intl %>% filter(lat >= 45, lat <= 49, lon >= -125, lon <= -116.5, 
                                         mag >= 2, date >= as.Date("1985-01-01"),
                                         is.na(remove), potential_duplicates_20210412 == F,
                                         !scheme2020 %in% c(1, 2))
maxT.pnw.Nt <- as.numeric(as.Date("2019-01-01")-as.Date("1985-01-01") + 1)

crustal.Nt.cat.noMCswarms <- pnw.Nt.cat.noMCswarms %>% filter(Depth2ModS > 10| is.na(Depth2ModS) )
crustal.Nt.cat.noMCswarms <- crustal.Nt.cat.noMCswarms[order(crustal.Nt.cat.noMCswarms$time),]
crustal.Nt.cat.noMCswarms$time <- crustal.Nt.cat.noMCswarms$time - min(crustal.Nt.cat.noMCswarms$time)
which.crustal.Nt.cat.noMCswarms <- which(intl.t$id %in% crustal.Nt.cat.noMCswarms$id)
crustal.Nt.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustal.Nt.cat.noMCswarms+2)]

# Fix alpha
runSTBETASPlotsPNW(cat=crustal.Nt.cat.noMCswarms, 
                   cat.label = "crustal-Nt-noMCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=5000, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/complete",
                   plot.title="PNW North Target (45-49N) Crustal Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 


###### Baseline runs, with crustal, no MCswarms, Nt region, & alpha fixed to SS
# Fix alpha.
runSTBETASPlotsPNW(cat=crustal.Nt.cat.noMCswarms, 
                   cat.label = "crustal-Nt-noMCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2, these.init.vals=cat2inits, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA/baseline",
                   plot.title="PNW North Target (45-49N) Crustal Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 

# Free alpha
runSTBETASPlotsPNW(cat=crustal.Nt.cat.noMCswarms, 
                   cat.label = "crustal-Nt-noMCswarms-Nt-cat2initsSS-freealpha", 
                   nsims=10500, this.M0=2.0, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=F,
                   cat.Gevals.mat = crustal.Nt.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA/baseline",
                   plot.title="PNW North Target (45-49N) Crustal Events M2+, 1985-2018, Low Proposal Sds (=0.05), Free Alpha ") 

##### Initial values
runSTBETASPlotsPNW(cat=crustal.Nt.cat.noMCswarms, 
                   cat.label = "crustal-Nt-noMCswarms-Nt-badinits3-fixalpha", 
                   nsims=10500, this.M0=2.0,these.init.vals=bad.inits3, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA/baseline",
                   plot.title="PNW North Target (45-49N) Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha SS") 

runSTBETASPlotsPNW(cat=crustal.Nt.cat.noMCswarms, 
                   cat.label = "crustal-noMCswarms-Nt-crazyinits1-fixalpha", 
                   nsims=10500, this.M0=2.0,these.init.vals=crazy.inits1, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA/baseline",
                   plot.title="crazyinits1, PNW North Target (45-49N) Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha SS") 

runSTBETASPlotsPNW(cat=crustal.Nt.cat.noMCswarms, 
                   cat.label = "crustal-Nt-noMCswarms-Nt-cat3inits-fixalpha", 
                   nsims=10500, this.M0=2.0,these.init.vals=cat3inits, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA/baseline",
                   plot.title="PNW North Target (45-49N) Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha SS") 



##### Swarms
pnw.Nt.cat.all <- intl %>% filter(lat >= 45, lat <= 49, lon >= -125, lon <= -116.5, 
                                  mag >= 2, date >= as.Date("1985-01-01"),
                                  is.na(remove), potential_duplicates_20210412 == F)

pnw.Nt.cat.all <- pnw.Nt.cat.all[order(pnw.Nt.cat.all$time),]
pnw.Nt.cat.all$time <- pnw.Nt.cat.all$time - min(pnw.Nt.cat.all$time)
crustal.Nt.cat.all <- pnw.Nt.cat.all %>% filter(Depth2ModS > 10| is.na(Depth2ModS) )
which.crustal.Nt.cat.all <- which(intl.t$id %in% crustal.Nt.cat.all$id)
crustal.Nt.all.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustal.Nt.cat.all+2)]
runSTBETASPlotsPNW(cat=crustal.Nt.cat.all, 
                   cat.label = "crustal-Nt-all-Nt-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.all.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA",
                   plot.title="PNW North Target (45-49N) Crustal Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 

pnw.Nt.cat.noCswarms <- intl %>% filter(lat >= 45, lat <= 49, lon >= -125, lon <= -116.5, 
                                        mag >= 2, date >= as.Date("1985-01-01"),
                                        is.na(remove), potential_duplicates_20210412 == F,
                                        !scheme2020 %in% c(1))
pnw.Nt.cat.noCswarms <- pnw.Nt.cat.noCswarms[order(pnw.Nt.cat.noCswarms$time),]
pnw.Nt.cat.noCswarms$time <- pnw.Nt.cat.noCswarms$time - min(pnw.Nt.cat.noCswarms$time)
crustal.Nt.cat.noCswarms <- pnw.Nt.cat.noCswarms %>% filter(Depth2ModS > 10| is.na(Depth2ModS) )
which.crustal.Nt.cat.noCswarms <- which(intl.t$id %in% crustal.Nt.cat.noCswarms$id)
crustal.Nt.noCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustal.Nt.cat.noCswarms+2)]
runSTBETASPlotsPNW(cat=crustal.Nt.cat.noCswarms, 
                   cat.label = "crustal-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2.0, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA",
                   plot.title="PNW North Target (45-49N) Crustal Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 

#### Tectonics
# TODO: remove this catalog? I don't think I ever model the crustal W only
crustal.W.Nt.cat.noMCswarms <- pnw.Nt.cat.noMCswarms %>% filter(Depth2ModS > 10 )
which.crustal.W.Nt.cat.noMCswarms <- which(intl.t$id %in% crustal.W.Nt.cat.noMCswarms$id)
crustal.W.Nt.cat.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustal.W.Nt.cat.noMCswarms+2)]
intl.Nt.W.box <- latlong2grid(rbind(c(-125, 49), c(-121, 45)))
intl.Nt.W.box$x <- intl.Nt.W.box$x - min(intl.Nt.W.box$x)
intl.Nt.W.box$y <- intl.Nt.W.box$y - min(intl.Nt.W.box$y)
intl.Nt.W.area <- intl.Nt.W.box$x[2]*intl.Nt.W.box$y[1]


deep.Nt.cat.noMCswarms <- pnw.Nt.cat.noMCswarms %>% filter(ABO_ModS == "Below" )
which.deep.Nt.cat.noMCswarms <- which(intl.t$id %in% deep.Nt.cat.noMCswarms$id)
deep.Nt.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.deep.Nt.cat.noMCswarms+2)]


# Fix alpha
runSTBETASPlotsPNW(cat=deep.Nt.cat.noMCswarms, 
                   cat.label = "deep-Nt-noMCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.W.area, this.isAlphaFixed=T, # Note West area only
                   cat.Gevals.mat = deep.Nt.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA",
                   plot.title="PNW North Target (45-49N) deep Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 

deep.Nt.cat.noCswarms <- pnw.Nt.cat.noCswarms %>% filter(ABO_ModS == "Below" )
which.deep.Nt.cat.noCswarms <- which(intl.t$id %in% deep.Nt.cat.noCswarms$id)
deep.Nt.noCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.deep.Nt.cat.noCswarms+2)]
runSTBETASPlotsPNW(cat=deep.Nt.cat.noCswarms, 
                   cat.label = "deep-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.W.area, this.isAlphaFixed=T, # Note West area only
                   cat.Gevals.mat = deep.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA",
                   plot.title="PNW North Target (45-49N) deep Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 

pnw.Nt.cat.all <- intl %>% filter(lat >= 45, lat <= 49, lon >= -125, lon <= -116.5, 
                                  mag >= 2, date >= as.Date("1985-01-01"),
                                  is.na(remove), potential_duplicates_20210412 == F)

deep.Nt.cat.all <- pnw.Nt.cat.all %>% filter(ABO_ModS == "Below")
which.deep.Nt.cat.all <- which(intl.t$id %in% deep.Nt.cat.all$id)
deep.Nt.all.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.deep.Nt.cat.all+2)]
runSTBETASPlotsPNW(cat=deep.Nt.cat.all, 
                   cat.label = "deep-Nt-all-Nt-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = deep.Nt.all.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA",
                   plot.title="PNW North Target (45-49N) Deep Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 
deep.Nt.cat.noCswarms <- pnw.Nt.cat.noCswarms %>% filter(ABO_ModS == "Below" )
which.deep.Nt.cat.noCswarms <- which(intl.t$id %in% deep.Nt.cat.noCswarms$id)
deep.Nt.noCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.deep.Nt.cat.noCswarms+2)]
runSTBETASPlotsPNW(cat=deep.Nt.cat.noCswarms, 
                   cat.label = "deep-Nt-noCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T, # Note West area only
                   cat.Gevals.mat = deep.Nt.noCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA",
                   plot.title="PNW North Target (45-49N) deep Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 


# Crustal+ Deep together
crustaldeep.Nt.cat.noMCswarms <- pnw.Nt.cat.noMCswarms %>% filter(Depth2ModS > 10| is.na(Depth2ModS) | ABO_ModS == "Below" )
which.crustaldeep.Nt.cat.noMCswarms <- which(intl.t$id %in% crustaldeep.Nt.cat.noMCswarms$id)
crustaldeep.Nt.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustaldeep.Nt.cat.noMCswarms+2)]

# Fix alpha
runSTBETASPlotsPNW(cat=crustaldeep.Nt.cat.noMCswarms, 
                   cat.label = "crustaldeep-Nt-noMCswarms-Nt-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2, these.init.vals=cat2inits.ss, 
                   this.spat.area = intl.Nt.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustaldeep.Nt.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.Nt, this.plot.loc="pnw/BSSA",
                   plot.title="PNW North Target (45-49N) crustal+deep Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 

########
#### PNW and SZ PRIORS
########
# Crustal and deep with priors
sz.lit <- read.csv("data/etas-estimates-subduction-zones.csv", 
                   colClasses = c(rep("character", 8), rep("numeric", 16), "character"))


mu.SZ.Unif.hypers <- c(min(sz.lit$mulow, na.rm=T), max(sz.lit$muhigh, na.rm=T))
K.SZ.Unif.hypers <- c(min(sz.lit$Klow, na.rm=T), max(sz.lit$Khigh, na.rm=T))
alpha.SZ.Unif.hypers <- c(min(sz.lit$alphalow, na.rm=T), max(sz.lit$alphahigh, na.rm=T))
c.SZ.Unif.hypers <- c(min(sz.lit$clow, na.rm=T), max(sz.lit$chigh, na.rm=T))
p.SZ.Unif.hypers <- c(min(sz.lit$plow, na.rm=T), max(sz.lit$phigh, na.rm=T))
d.SZ.Unif.hypers <- c(min(sz.lit$dlow, na.rm=T), max(sz.lit$dhigh, na.rm=T))
q.SZ.Unif.hypers <- c(min(sz.lit$qlow, na.rm=T), max(sz.lit$qhigh, na.rm=T))
gamma.SZ.Unif.hypers <- c(min(sz.lit$gammalow, na.rm=T), max(sz.lit$gammahigh, na.rm=T))
SZ.Unif.pgammas.vec <- rep(FALSE, 7)

# SZ.inits <- c(0.1, 0.02, 2.3, 0.005, 1.08, 0.32, 1.6)
SZ.inits <- c(0.1, 0.02, 2, 0.005, 1.08, 0.32, 1.6)

runSTBETASPlotsPNWPriors(cat=crustal.Nt.cat.noMCswarms, 
                         cat.label = "crustal-noMCswarms-Nt-SZinits-fixalpha-JGpriors", 
                         these.init.vals=SZ.inits, 
                         cat.Gevals.mat=crustal.Nt.noMCswarms.Gevals.mat, 
                         this.plot.loc="pnw/BSSA/priors",
                         this.spat.area=intl.Nt.area, nsims=10500, this.M0=2, 
                         this.maxT=maxT.pnw.Nt, this.isAlphaFixed=T, 
                         this.pgammas_vec = SZ.Unif.pgammas.vec,
                         these.mu.hypers=mu.SZ.Unif.hypers, 
                         these.K.hypers=K.SZ.Unif.hypers, 
                         these.alpha.hypers=alpha.SZ.Unif.hypers,
                         these.c.hypers=c.SZ.Unif.hypers, 
                         these.p.hypers=p.SZ.Unif.hypers,
                         these.d.hypers=d.SZ.Unif.hypers, 
                         these.q.hypers=q.SZ.Unif.hypers,
                         this.JG.region = "pnw",
                         plot.title="Cat2 Inits, PNW North (45-49N) Events M2+, 2004-2018, Low Proposal Sds (=0.05), \nAlpha = 2.3, Ginterp, Using JG's PNW Priors") 
runSTBETASPlotsPNWPriors(cat=deep.Nt.cat.noMCswarms, 
                         cat.label = "deep-noMCswarms-Nt-SZinits-fixalpha-JGpriors", 
                         these.init.vals=SZ.inits, 
                         cat.Gevals.mat=deep.Nt.noMCswarms.Gevals.mat, 
                         this.plot.loc="pnw/BSSA/priors",
                         this.spat.area=intl.Nt.area, nsims=10500, this.M0=2, 
                         this.maxT=maxT.pnw.Nt, this.isAlphaFixed=T, 
                         this.pgammas_vec = SZ.Unif.pgammas.vec,
                         these.mu.hypers=mu.SZ.Unif.hypers, 
                         these.K.hypers=K.SZ.Unif.hypers, 
                         these.alpha.hypers=alpha.SZ.Unif.hypers,
                         these.c.hypers=c.SZ.Unif.hypers, 
                         these.p.hypers=p.SZ.Unif.hypers,
                         these.d.hypers=d.SZ.Unif.hypers, 
                         these.q.hypers=q.SZ.Unif.hypers,
                         this.JG.region = "pnw",
                         plot.title="Cat2 Inits, PNW North (45-49N) Events M2+, 2004-2018, Low Proposal Sds (=0.05), \nAlpha = 2.3, Ginterp, Using JG's PNW Priors") 

# mu, K, alpha, c, p, d, q
# Set PNW priors from Paul Bodin
mu.a.crustal <- 9.22
mu.b.crustal <- 74.7

# These are for alpha fixed to 2.3
K.a.crustal <- 1.45
K.b.crustal <- 285

# # These are for alpha fixed to 1.7
# K.a.crustal <- 1.49
# K.b.crustal <- 26.4

alpha.min.crustal <- 0.3
alpha.max.crustal <- 2.5

c.a.crustal <- 0.717
c.b.crustal <- 35.8

p.min.crustal <- 0.407 
p.max.crustal <- 2.643

d.a.crustal <- 2.477884
d.b.crustal <- 0.05948049

q.min.crustal <- 1.012 
q.max.crustal <- 2.073
PB.pgammas.vec <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
runSTBETASPlotsPNWPriors(cat=crustal.Nt.cat.noMCswarms, 
                         cat.label = "crustal-noMCswarms-Nt-cat2initsSS-fixalpha-PBpriors", 
                         these.init.vals=cat2inits.ss, 
                         cat.Gevals.mat=crustal.Nt.noMCswarms.Gevals.mat, 
                         this.plot.loc="pnw/BSSA/priors",
                         this.spat.area=intl.Nt.area, nsims=10500, this.M0=2, 
                         this.maxT=maxT.pnw.Nt, this.isAlphaFixed=T, 
                         this.pgammas_vec = PB.pgammas.vec,
                         these.mu.hypers=c(mu.a.crustal, mu.b.crustal), 
                         these.K.hypers=c(K.a.crustal, K.b.crustal), 
                         these.alpha.hypers=c(alpha.min.crustal, alpha.max.crustal),
                         these.c.hypers=c(c.a.crustal, c.b.crustal), 
                         these.p.hypers=c(p.min.crustal, p.max.crustal),
                         these.d.hypers=c(d.a.crustal, d.b.crustal), 
                         these.q.hypers=c(q.min.crustal, q.max.crustal),
                         this.PB.region = "crustal",
                         plot.title="Cat2 Inits, PNW South (42-45N) Events M2+, 2004-2018, Low Proposal Sds (=0.05), \nAlpha = 2.3, Ginterp, Using JG's PNW Priors") 
mu.a.deep <- 160
mu.b.deep <- 5350

# These are for alpha fixed to 2.3
K.a.deep <- round(0.698, 2)
K.b.deep <- 209

# # These are for alpha fixed to 1.7
# K.a.deep <- 0.746
# K.b.deep <- 19.8

alpha.min.deep <- 0
alpha.max.deep <- 1.4

# Same for deep and crustal
c.a.deep <- round(0.717, 2)
c.b.deep <- 35.8

# Same for deep and crustal
p.min.deep <- round(0.407 , 2)
p.max.deep <- round(2.643, 2)

d.a.deep <- round(0.8877908, 2)
d.b.deep <- round(0.01382139, 2)

q.min.deep <- round(1.024, 2)
q.max.deep <- round(2.775, 2)
runSTBETASPlotsPNWPriors(cat=deep.Nt.cat.noMCswarms, 
                         cat.label = "deep-noMCswarms-Nt-cat2initsSS-freealpha-PBpriors", 
                         these.init.vals=cat2inits.ss, 
                         cat.Gevals.mat=deep.Nt.noMCswarms.Gevals.mat, 
                         this.plot.loc="pnw/BSSA/priors",
                         this.spat.area=intl.Nt.area, nsims=10500, this.M0=2, 
                         this.maxT=maxT.pnw.Nt, this.isAlphaFixed=F, 
                         this.pgammas_vec = PB.pgammas.vec,
                         these.mu.hypers=c(mu.a.deep, mu.b.deep), 
                         these.K.hypers=c(K.a.deep, K.b.deep), 
                         these.alpha.hypers=c(alpha.min.deep, alpha.max.deep),
                         these.c.hypers=c(c.a.deep, c.b.deep), 
                         these.p.hypers=c(p.min.deep, p.max.deep),
                         these.d.hypers=c(d.a.deep, d.b.deep), 
                         these.q.hypers=c(q.min.deep, q.max.deep),
                         this.PB.region = "deep",
                         plot.title="Cat2 Inits, PNW South (42-45N) Events M2+, 2004-2018, Low Proposal Sds (=0.05), \nFree Alpha, Ginterp, Using JG's PNW Priors") 



##########
# South
cat2inits.ss.M2 <- cat2inits.ss
cat2inits.ss.M2[3] <- 0.882*(log(10))
pnw.St.cat.noCswarms <- intl %>% filter(lat >= 42, lat < 45, lon >= -125, lon <= -116.5, 
                                        mag >= 2, date >= as.Date("2004-01-01"),
                                        is.na(remove), potential_duplicates_20210412 == F,
                                        !scheme2020 %in% c(1))
pnw.St.cat.noCswarms$time <- pnw.St.cat.noCswarms$time - min(pnw.St.cat.noCswarms$time)
time.zero.pnw <- as.POSIXct("1970-01-01 01:01:01", format = "%Y-%m-%d %H:%M:%S",  tz = "UTC")
maxT.pnw.St <- as.numeric(as.Date("2019-01-01")-as.Date("2004-01-01") + 1)

which.pnw.St.cat.noCswarms <- which(intl.t$id %in% pnw.St.cat.noCswarms$id)
pnw.St.cat.noCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.pnw.St.cat.noCswarms+2)]

intl.St.box <- latlong2grid(rbind(c(-125, 45), c(-116.5, 42)))
intl.St.box$x <- intl.St.box$x - min(intl.St.box$x)
intl.St.box$y <- intl.St.box$y - min(intl.St.box$y)
intl.St.area <- intl.St.box$x[2]*intl.St.box$y[1]

pnw.St.cat.noMCswarms <- intl %>% filter(lat >= 42, lat < 45, lon >= -125, lon <= -116.5,
                                         mag >= 2, date >= as.Date("2004-01-01"),
                                         is.na(remove), potential_duplicates_20210412 == F,
                                         !scheme2020 %in% c(1, 2))
pnw.St.cat.noMCswarms$time <- pnw.St.cat.noMCswarms$time - min(pnw.St.cat.noMCswarms$time)
which.pnw.St.cat.noMCswarms <- which(intl.t$id %in% pnw.St.cat.noMCswarms$id)
pnw.St.cat.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.pnw.St.cat.noMCswarms+2)]

crustal.St.cat.noMCswarms <- pnw.St.cat.noMCswarms %>% filter(Depth2ModS > 10| is.na(Depth2ModS) )
which.crustal.St.cat.noMCswarms <- which(intl.t$id %in% crustal.St.cat.noMCswarms$id)
crustal.St.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustal.St.cat.noMCswarms+2)]

# Fix alpha
runSTBETASPlotsPNW(cat=crustal.St.cat.noMCswarms, 
                   cat.label = "crustal-St-noMCswarms-St-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2, these.init.vals=cat2inits.ss.M2, 
                   this.spat.area = intl.St.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.St.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.St, this.plot.loc="pnw/BSSA/south",
                   plot.title="PNW North Target (45-49N) Crustal Events M2+, 1985-2018, Low Proposal Sds (=0.05), Alpha = 2.3") 

runSTBETASPlotsPNWPriors(cat=crustal.St.cat.noMCswarms, 
                         cat.label = "crustal-noMCswarms-St-SZinits-fixalpha-JGpriors", 
                         these.init.vals=SZ.inits, 
                         cat.Gevals.mat=crustal.St.noMCswarms.Gevals.mat, 
                         this.plot.loc="pnw/BSSA/priors",
                         this.spat.area=intl.St.area, nsims=10500, this.M0=2, 
                         this.maxT=maxT.pnw.Nt, this.isAlphaFixed=T, 
                         this.pgammas_vec = SZ.Unif.pgammas.vec,
                         these.mu.hypers=mu.SZ.Unif.hypers, 
                         these.K.hypers=K.SZ.Unif.hypers, 
                         these.alpha.hypers=alpha.SZ.Unif.hypers,
                         these.c.hypers=c.SZ.Unif.hypers, 
                         these.p.hypers=p.SZ.Unif.hypers,
                         these.d.hypers=d.SZ.Unif.hypers, 
                         these.q.hypers=q.SZ.Unif.hypers,
                         this.JG.region = "pnw",
                         plot.title="Cat2 Inits, PNW North (45-49N) Events M2+, 2004-2018, Low Proposal Sds (=0.05), \nAlpha = 2.3, Ginterp, Using JG's PNW Priors") 

runSTBETASPlotsPNWPriors(cat=crustal.St.cat.noMCswarms, 
                         cat.label = "crustal-noMCswarms-St-cat2initsSS-fixalpha-PBpriors", 
                         these.init.vals=cat2inits.ss, 
                         cat.Gevals.mat=crustal.St.noMCswarms.Gevals.mat, 
                         this.plot.loc="pnw/BSSA/priors",
                         this.spat.area=intl.St.area, nsims=10500, this.M0=2, 
                         this.maxT=maxT.pnw.Nt, this.isAlphaFixed=T, 
                         this.pgammas_vec = PB.pgammas.vec,
                         these.mu.hypers=c(mu.a.crustal, mu.b.crustal), 
                         these.K.hypers=c(K.a.crustal, K.b.crustal), 
                         these.alpha.hypers=c(alpha.min.crustal, alpha.max.crustal),
                         these.c.hypers=c(c.a.crustal, c.b.crustal), 
                         these.p.hypers=c(p.min.crustal, p.max.crustal),
                         these.d.hypers=c(d.a.crustal, d.b.crustal), 
                         these.q.hypers=c(q.min.crustal, q.max.crustal),
                         this.PB.region = "crustal",
                         plot.title="Cat2 Inits, PNW South (42-45N) Events M2+, 2004-2018, Low Proposal Sds (=0.05), \nAlpha = 2.3, Ginterp, Using JG's PNW Priors") 


###
# Mc=2.3
cat2inits.ss.M2.3 <- cat2inits
cat2inits.ss.M2.3[3] <- 0.981*(log(10))
pnw.St.comp1.cat.noMCswarms <- intl %>% filter(lat >= 42, lat < 45, lon >= -125, lon <= -116.5, 
                                               mag >= 2.3, date >= as.Date("2004-01-01"),
                                               is.na(remove), potential_duplicates_20210412 == F,
                                               !scheme2020 %in% c(1, 2))
pnw.St.comp1.cat.noMCswarms$time <- pnw.St.comp1.cat.noMCswarms$time - min(pnw.St.comp1.cat.noMCswarms$time)
time.zero.pnw <- as.POSIXct("1970-01-01 01:01:01", format = "%Y-%m-%d %H:%M:%S",  tz = "UTC")
maxT.pnw.St <- as.numeric(as.Date("2019-01-01")-as.Date("2004-01-01") + 1)

which.pnw.St.comp1.cat.noMCswarms <- which(intl.t$id %in% pnw.St.comp1.cat.noMCswarms$id)
pnw.St.comp1.cat.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.pnw.St.comp1.cat.noMCswarms+2)]

intl.St.comp1.box <- latlong2grid(rbind(c(-125, 45), c(-116.5, 42)))
intl.St.comp1.box$x <- intl.St.comp1.box$x - min(intl.St.comp1.box$x)
intl.St.comp1.box$y <- intl.St.comp1.box$y - min(intl.St.comp1.box$y)
intl.St.comp1.area <- intl.St.comp1.box$x[2]*intl.St.comp1.box$y[1]

crustal.St.comp1.cat.noMCswarms <- pnw.St.comp1.cat.noMCswarms %>% filter(Depth2ModS > 10| is.na(Depth2ModS) )
which.crustal.St.comp1.cat.noMCswarms <- which(intl.t$id %in% crustal.St.comp1.cat.noMCswarms$id)
crustal.St.comp1.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustal.St.comp1.cat.noMCswarms+2)]

# Fix alpha
runSTBETASPlotsPNW(cat=crustal.St.comp1.cat.noMCswarms, 
                   cat.label = "crustal-St-comp1-noMCswarms-St-comp1-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2.3, these.init.vals=cat2inits.ss.M2.3, 
                   this.spat.area = intl.St.comp1.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.St.comp1.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.St, this.plot.loc="pnw/BSSA/south",
                   plot.title="PNW South Target (42-45N) Crustal Events M2.3+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar") 

###
# Mc=2.5
cat2inits.ss.M2.5 <- cat2inits
cat2inits.ss.M2.5[3] <- 1.07*(log(10))

pnw.St.comp2.cat.noMCswarms <- intl %>% filter(lat >= 42, lat < 45, lon >= -125, lon <= -116.5, 
                                               mag >= 2.5, date >= as.Date("2004-01-01"),
                                               is.na(remove), potential_duplicates_20210412 == F,
                                               !scheme2020 %in% c(1,2))
pnw.St.comp2.cat.noMCswarms$time <- pnw.St.comp2.cat.noMCswarms$time - min(pnw.St.comp2.cat.noMCswarms$time)
time.zero.pnw <- as.POSIXct("1970-01-01 01:01:01", format = "%Y-%m-%d %H:%M:%S",  tz = "UTC")
maxT.pnw.St <- as.numeric(as.Date("2019-01-01")-as.Date("2004-01-01") + 1)

which.pnw.St.comp2.cat.noMCswarms <- which(intl.t$id %in% pnw.St.comp2.cat.noMCswarms$id)
pnw.St.comp2.cat.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.pnw.St.comp2.cat.noMCswarms+2)]

intl.St.comp2.box <- latlong2grid(rbind(c(-125, 45), c(-116.5, 42)))
intl.St.comp2.box$x <- intl.St.comp2.box$x - min(intl.St.comp2.box$x)
intl.St.comp2.box$y <- intl.St.comp2.box$y - min(intl.St.comp2.box$y)
intl.St.comp2.area <- intl.St.comp2.box$x[2]*intl.St.comp2.box$y[1]

crustal.St.comp2.cat.noMCswarms <- pnw.St.comp2.cat.noMCswarms %>% filter(Depth2ModS > 10| is.na(Depth2ModS) )
which.crustal.St.comp2.cat.noMCswarms <- which(intl.t$id %in% crustal.St.comp2.cat.noMCswarms$id)
crustal.St.comp2.noMCswarms.Gevals.mat <- pnw.Gevals.mat[, c(1:2, which.crustal.St.comp2.cat.noMCswarms+2)]

# Fix alpha
runSTBETASPlotsPNW(cat=crustal.St.comp2.cat.noMCswarms, 
                   cat.label = "crustal-St-comp2-noMCswarms-St-comp2-cat2initsSS-fixalpha", 
                   nsims=10500, this.M0=2.5, these.init.vals=cat2inits.ss.M2.5, 
                   this.spat.area = intl.St.comp2.area, this.isAlphaFixed=T,
                   cat.Gevals.mat = crustal.St.comp2.noMCswarms.Gevals.mat, this.maxT=maxT.pnw.St, this.plot.loc="pnw/BSSA/south",
                   plot.title="PNW South Target (42-45N) Crustal Events M2.5+, 1985-2018, Low Proposal Sds (=0.05), Self-Similar") 

