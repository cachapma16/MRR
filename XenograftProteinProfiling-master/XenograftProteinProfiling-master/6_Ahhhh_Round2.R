setwd("D:/Google Drive/Research/Notebook/2017/07/GBM-6_10 Compiled/")
library(tidyverse)
library(ggthemes)

dat <- read_csv("compiledLabeled.csv")
dat <- select(dat, -Shift.1, -Shift.2, - Group)
dat <- filter(dat, !grepl("4E", Target))
dat$Id <- dat$Ring %% 4

head(dat)

dat.dcast <<- dcast(dat,
                    Treatment + `Time Point` + Cell.Line + Id ~ Target,
                    value.var = "Net Shift",
                    fun.aggregate = mean)
head(dat.dcast)
dat.dcast$Cell.Line <- factor(dat.dcast$Cell.Line)

plot.1 <- ggplot(dat.dcast, 
               aes(x = `pAkt Ser473`, 
                   y = `pMAPK Thr202/Tyr204`, 
                   color = Treatment,
                   shape = Cell.Line)) +
        geom_point(size = 6) + theme_few()

plot.1

ggsave(plot.1, filename = "Akt vs MAPK Flow cytometry esque.png", width = 8, height = 6)

head(dat.dcast)

plot.2 <- ggplot(dat.dcast, 
               aes(x = `pRb Ser807/11`, 
                   y = `pS6 Ser235/6`, 
                   color = Treatment,
                   shape = Cell.Line)) +
        geom_point(size = 6) + theme_few()

ggsave(plot.2, filename = "Rb vs S6 Flow cytometry esque.png", width = 8, height = 6)

plot(dat.dcast$`hydroxy-HIF Pro564` ~ dat.dcast$`pAkt Ser473`)
