setwd("D:/Box Sync/Data/")
library(tidyverse)
library(RColorBrewer)
library(reshape2)

dat <- read_csv("compiledLabeled.csv")
dat <- select(dat, -Experiment, -Shift.1, -Shift.2, -Experiment, -Replicate)

dat.1 <- filter(dat, !grepl("p53|Abl", Target))
dat.2 <- filter(dat.1, !grepl("Abl|p53|HIF", Target))

dat.6 <- filter(dat.2, `Cell Line` == 6)
dat.26 <- filter(dat.2, `Cell Line` == 26)

casting.1 <- dcast(dat.1, Treatment + `Cell Line` + `Time Point` ~ 
                        Target, mean, value.var = "Net Shift")

casting.2 <- dcast(dat.2, Treatment + `Cell Line` + `Time Point` ~ 
                        Target, mean, value.var = "Net Shift")

casting.6 <- dcast(dat.6, Treatment + `Cell Line` + `Time Point` ~ 
                           Target, mean, value.var = "Net Shift")

casting.26 <- dcast(dat.26, Treatment + `Cell Line` + `Time Point` ~ 
                           Target, mean, value.var = "Net Shift")

rownames(casting.1) = paste(casting.1$Treatment, casting.1$`Cell Line`, 
                        casting.1$`Time Point`, sep = " ")

rownames(casting.2) = paste(casting.2$Treatment, casting.2$`Cell Line`, 
                            casting.2$`Time Point`, sep = " ")

rownames(casting.6) = paste(casting.6$Treatment, casting.6$`Cell Line`, 
                            casting.6$`Time Point`, sep = " ")

rownames(casting.26) = paste(casting.26$Treatment, casting.26$`Cell Line`, 
                            casting.26$`Time Point`, sep = " ")

fig.1 <- ggplot(dat, aes(interaction(Treatment, factor(`Cell Line`), 
                        factor(`Time Point`)), Target, fill = `Net Shift`)) + 
                geom_tile() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))
fig.1


casting.1[,c(1,2,3)] <- NULL
casting.2[,c(1,2,3)] <- NULL
casting.6[,c(1,2,3)] <- NULL
casting.26[,c(1,2,3)] <- NULL

casting.1.m <- data.matrix(casting.1)
casting.2.m <- data.matrix(casting.2)
casting.6.m <- data.matrix(casting.6)
casting.26.m <- data.matrix(casting.26)

hmcol <- colorRampPalette(c("darkblue", "white", "darkred"))(256)

heatmap(casting.1.m, scale = "col", col = hmcol)
dev.print(png,"heatmap_all.png",width=8,height=6,units="in",res=300)
dev.off()

heatmap(casting.1.m, scale = "col", col = hmcol)
dev.print(pdf,"heatmap_all.pdf",width=8,height=6)
dev.off()

heatmap(casting.2.m, scale = "col", col = hmcol)
dev.print(png,"heatmap_all_2.png",width=8,height=6,units="in",res=300)
dev.off()

heatmap(casting.2.m, scale = "col", col = hmcol)
dev.print(pdf,"heatmap_all_2.pdf",width=8,height=6)
dev.off()

heatmap(casting.6.m, scale = "col", col = hmcol)
dev.print(png,"heatmap_GBM6.png",width=8,height=6,units="in",res=300)
dev.off()

heatmap(casting.6.m, scale = "col", col = hmcol)
dev.print(pdf,"heatmap_GBM6.pdf",width=8,height=6)
dev.off()

heatmap(casting.26.m, scale = "col", col = hmcol)
dev.print(png,"heatmap_GBM26.png",width=8,height=6,units="in",res=300)
dev.off()

heatmap(casting.26.m, scale = "col", col = hmcol)
dev.print(pdf,"heatmap_GBM26.pdf",width=8,height=6)
dev.off()


install.packages("igraph")
install.packages("network")
install.packages("sna")
install.packages("ndtv")


install.packages("corrplot")
library(corrplot)
casting.1.m
holder <- cor(casting.1.m)
corrplot(holder, type ="upper")
View(casting.1.m)


# igraph package ----------------------------------------------------------

library(igraph)
library(corrplot)
library(Hmisc)

dat <- read_csv("compiledLabeled.csv")
dat <- select(dat, -Experiment, -Shift.1, -Shift.2, -Experiment, -Replicate, 
              -Step, -Ring, -Group)
dat <- filter(dat, !grepl("p53|Abl", Target))

PlotNetwork <- function(cellLine, treatment){
        dat.cntl <- filter(dat, Treatment == "(-)-Serum" & `Cell Line` == cellLine
                           & `Time Point` == 1)
        
        dat.1 <- filter(dat, Treatment == treatment & `Cell Line` == cellLine)
        
        dat.1 <- rbind(dat.1, dat.cntl)

        m.1 <- dcast(dat.1, Treatment + `Cell Line` + `Time Point` ~ Target, 
                                mean, value.var = "Net Shift")
        rownames(m.1) = paste(m.1$Treatment, m.1$`Cell Line`,
                              m.1$`Time Point`, sep = " ")
        
        m.1[,c(1,2,3)] <- NULL
        m.1 <- as.matrix(m.1)

        holder.t <- cor(scale(t(m.1)))
        holder <- cor(scale(m.1))
        holder.r <- rcorr(scale(m.1))
        
        custom.layout <- function(g, ...) {
                # layout.drl(g, weights = E(g)$weight, ...) # For bigger graphs
                layout.fruchterman.reingold(g, weights = E(g)$weight, ...)
        }
        
        g <- graph.adjacency(holder, mode = "undirected", weighted = TRUE, diag = FALSE)
        
        # Format edges
        E(g)$cor <- E(g)$weight
        E(g)$weight <- abs(E(g)$cor)
        E(g)$color <- ifelse(E(g)$cor < 0, "blue", "red")
        #E(g)$width <- 3*atanh(E(g)$weight)
        E(g)$width <- atanh(E(g)$weight)
        
        # Format vertices
        V(g)$size <- 3*abs(rowSums(holder))
        V(g)$color <- "grey"
        V(g)$label.color <- "black"
        
        plot(g, layout = custom.layout)
}


PlotNetwork(cellLine = c(26, 6), treatment = "Palbociclib")

# Chose the layout function


# Do the plot
plot(g, layout = custom.layout)
plot(g)
corrplot(holder, type = "upper")
dat.1
View(casting.1)

# Okay, so I think what I need to do is to go through each treatment and get
# correlationgs for each of them and then generate this network. I also should
# read the documentation to kinda know what the hell I'm doing. What about
# integrating into a single plot coloring by treatment


# Single Value Decomposition ----------------------------------------------

svd1 <- svd(scale(casting.1.m))

par(mfrow = c(1, 3))
heatmap(casting.6.m, scale = "col", col = hmcol)
plot(svd1$u[, 1], nrow(casting.1.m):1, xlab = "Row", ylab = "First left singular vector",
     pch = 19)
plot(svd1$v[, 1], xlab = "Column", ylab = "First right singular vector", pch = 19)


plot(svd1$d^2/sum(svd1$d^2), xlab = "Column", ylab = "Singular value", pch = 19)

corrplot(holder)


# Reshape -----------------------------------------------------------------

dat.test <- filter(dat, Treatment == "Palbociclib", `Time Point` == 1, 
                   `Cell Line` == 26)
View(dat.test)
