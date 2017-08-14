setwd("D:/Box Sync/Data/")
library(tidyverse)
library(ggthemes)
library(GGally)
library(Hmisc)

dat <- read_csv("compiledLabeled.csv")
dat <- select(dat, -Shift.1, -Shift.2, - Group)
dat <- filter(dat, !grepl("p53|Abl", Target))
dat$Id <- dat$Ring %% 4

head(dat)

dat.dcast <<- dcast(dat,
                    Treatment + `Time Point` + `Cell Line` + Id ~ Target,
                    value.var = "Net Shift",
                    fun.aggregate = mean)
head(dat.dcast)
dat.dcast$`Cell Line` <- factor(dat.dcast$`Cell Line`)

plot <- ggplot(dat.dcast, 
               aes(x = `pRb Ser807/11`, 
                   y = `pS6 Ser235/6`, 
                   color = Treatment,
                   shape = interaction(`Cell Line`, `Time Point`)))

plot + geom_point(size = 6) + theme_few()

dat.scale <- scale(dat.dcast[,5:17])
dat.scaled <- cbind(dat.dcast[,1:4], dat.scale)
colnames(dat.scaled)[3] <- "CellLine"
colnames(dat.scaled)[2] <- "TimePoint"
dat.scaled$TimePoint <- factor(dat.scaled$TimePoint)
colnames(dat.scaled)[8] <- "pGSK3B Ser9"

pairs(dat.dcast[5:17], col = factor(dat.dcast$`Cell Line`))
pairs(dat.scaled[5:17], col = factor(dat.scaled$`Cell Line`))

ggpairs(dat.scaled, aes(color = TimePoint))


dat.dcast$Txt <- dat.dcast$Treatment
dat.dcast$Treatment <- NULL

plot(dat.dcast$`hydroxy-HIF Pro564` ~ dat.dcast$`pAkt Ser473`)
