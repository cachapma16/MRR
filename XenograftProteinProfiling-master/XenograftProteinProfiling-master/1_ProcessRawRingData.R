GetName <- function(){
        # get the filename from the current working directory
        directory <- basename(getwd())
        
        # directory naming from MRR: "CHIPNAME_gaskGASK_DATE"
        # extracts and returns GASK from directory name
        name <- unlist(strsplit(directory, split = "_"))[2]
        
        # define name as global variable for use in other functions
        name <<- gsub('gask','',name) # removes "gask" from name
}

AggData <- function(loc = 'plots', filename = 'groupNames_XPP.csv') {
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # get information of chip layout from github repository
        if (!file.exists(filename)){
                git <- "https://raw.githubusercontent.com/JamesHWade/XenograftProteinProfiling/master/"
                url <- paste0(git, filename)
                filename <- basename(url)
                download.file(url, filename)
        }
        
        # read in recipe/chip layout
        recipe <- read_csv(filename)
        colnames(recipe)[1] <- "Target" # rename column to remove byte order mark
        targets <- recipe$Target
        
        # generate list of rings to analyze (gets all *.csv files)
        rings <- list.files(directory, pattern = ".csv", recursive = FALSE)
        idfile <- grepl("group", rings)
        removeFiles <- c("comments.csv", rings[idfile])
        rings <- rings[!rings %in% removeFiles]
        
        # add data to data frame corresponding for each ring in rings
        df <- lapply(rings, function(i){
                dat <- read_csv(i, col_names = FALSE)
                ringNum <- as.numeric(strsplit(i, "\\.")[[1]][1])
                recipeCol <- which(recipe$Ring == ringNum)
                tmp <- dat[,c(1,2)] # time and shift from raw data
                tmp$ring <- ringNum
                tmp$group <- recipe$Group[recipeCol]
                tmp$groupName <- as.character(recipe$Target[[recipeCol]])
                tmp$channel <- recipe$Channel[[recipeCol]]
                tmp$run <- name
                tmp$timePoint <- seq(1:nrow(dat))
                tmp
        })
        
        df <- bind_rows(df)
        
        # renames columns in df
        names(df) <- c("Time", "Shift", "Ring", "Group", "Target", "Channel",
                       "Experiment", "Time Point")
        
        # creates "plots" directory
        dir.create(loc, showWarnings = FALSE)
        
        # saves aggregated data with name_allRings.csv
        write_csv(df, paste(loc, '/', name, "_allRings.csv", sep=""))
}

SubtractControl <- function(loc = 'plots', ch, cntl){
        #load relevant libraries
        library(tidyverse)
        
        # get ring data and filter by channel
        dat <- read_csv(paste0(loc, "/", name, "_", "allRings.csv"))
        if (ch != "U"){dat <- filter(dat, Channel == ch)}
        dat <- filter(dat, Target != "Ignore")
        
        # get thermal control averages
        controls <- filter(dat, Target == cntl) %>% group_by(`Time Point`) %>%
                summarise_at("Shift", mean) %>% select(Shift) %>% unlist()
        dat$Cntl <- rep(controls, length(unique(dat$Ring)))
        
        # subtracts thermal controls from each ring
        dat.cntl <- mutate(dat, Shift = Shift - Cntl)
        
        # remove control column and control rings
        dat.cntl <- filter(dat, Target != cntl)
        dat.cntl$Cntl <- NULL
        
        write_csv(dat, paste(loc,"/", name, "_", cntl, "Control", "_ch", ch, 
                             ".csv", sep = ''))   
}

PlotRingData <- function(cntl, ch, loc = 'plots', splitPlot = FALSE){
        # loads relevant libraries and plot theme
        library(tidyverse)
        library(ggthemes)
        theme_set(theme_few(base_size = 16))
        
        # use thermally controlled data if desired
        if (cntl != "raw"){
                dat <- read_csv(paste(loc, "/", name, "_", cntl, "Control", 
                                      "_ch", ch,".csv", sep=''))
        } else if (cntl == "raw") {
                dat <- read_csv(paste(loc, "/", name, "_allRings.csv", sep=''))
                if (ch != "U") {dat <- filter(dat, Channel == ch)}
        }
        
        # configure plot and legend
        plots <- ggplot(dat, aes(x = Time, y = Shift,
                                 color = Target, group = Ring)) + 
                labs(x = "Time (min)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)")),
                     color = "Target") +
                geom_line() + 
                ggtitle(paste(name, "Ch:", ch, "Control:", cntl, sep = " "))
        
        # alternative plots with averaged clusters
        
        dat.2 <- dat %>% group_by(Target, `Time Point`) %>% 
                summarise_at(vars(Time, Shift), funs(mean, sd))
        
        plot2 <- ggplot(dat.2, aes(x = Time_mean, y = Shift_mean, 
                                   color = Target)) +
                geom_line() +
                labs(x = "Time (min)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)"))) +
                ggtitle(paste(name, "Ch:", ch, "Control:", cntl, sep = " "))
        
        plot3 <- plot2 + 
                geom_ribbon(aes(ymin = Shift_mean - Shift_sd,
                                ymax = Shift_mean + Shift_sd, 
                                linetype = NA),
                            fill = "slategrey", alpha = 1/8)
        
        if (splitPlot){
                plots <- plots + facet_grid(. ~ Channel)
        }
        
        # save plots
        ggsave(plot = plots, 
               file = paste0(loc, "/", name, "_", cntl, 
                             "Control_ch", ch, ".png"),
               width = 10, height = 6)
        ggsave(plot = plot2, 
               file = paste0(loc, "/", name, "_", cntl,
                             "Control", "_ch", ch, "_avg.png"),
               width = 10, height = 6)
        ggsave(plot = plot3, 
               file = paste0(loc, "/", name, "_", cntl,
                             "Control", "_ch", ch, "_avg.png"),
               width = 10, height = 6)
}

GetNetShifts <- function(cntl, ch, loc = 'plots', time1, time2, step = 1){
        # load relevant libraries
        library(tidyverse)
        
        # use thermally controlled data if desired
        if (cntl != "raw"){
                dat <- read_csv(paste0(loc, "/", name, "_", cntl, "Control", 
                                       "_ch", ch, ".csv"))
        } else {
                dat <- read_csv(paste0(loc, "/", name, "_", "allRings.csv"))
        }
        
        # generate list of rings and empty dataframe to store net shift data
        ringList <- unique(dat$Ring)

        # locations for each time is determined using which, min, and abs func
        dat.rings <- lapply(ringList, function(i){
                dat.ring <- filter(dat, Ring == i)
                time1.loc <- which.min(abs(dat.ring$Time - time1))
                time1.val <- dat.ring$Shift[time1.loc]
                time2.loc <- which.min(abs(dat.ring$Time - time2))
                time2.val <- dat.ring$Shift[time2.loc]
                ring <- i
                group <- unique(dat.ring$Group)
                target <- unique(dat.ring$Target)
                experiment <- unique(dat.ring$Experiment)
                channel <- unique(dat.ring$Channel)
                data.frame(i, group, target, time1.val,
                           time2.val, experiment, channel, step)
        })
        
        # renames dat.rings columns
        dat.rings <- bind_rows(dat.rings)
        names(dat.rings) <- c("Ring", "Group", "Target", "Shift.1", "Shift.2", 
                              "Experiment", "Channel", "Step")
        
        # calculate nat shift and create new column in dataframe
        dat.rings <- dat.rings %>% 
                mutate(`Net Shift` = Shift.1 - Shift.2)
        
        # save net shift data
        write_csv(dat.rings, paste0(loc, "/", name, "_netShifts_", cntl, "cntl_",
                                    "ch", ch, "_step", step, ".csv"))
}

PlotNetShifts <- function(cntl, ch, loc = 'plots', step = 1){
        # load relevant libraries and plot theme
        library(tidyverse)
        library(ggthemes)
        theme_set(theme_few(base_size = 16))

        dat <- read_csv(paste0(loc, "/", name, "_netShifts_", cntl, "cntl_",
                                       "ch", ch, "_step", step, ".csv"))
        
        # configure plot and legend
        dat.nothermal <- filter(dat, Target != "thermal")
        
        plots <- ggplot(dat.nothermal, 
                        aes(x = Target, y = `Net Shift`, fill = Target)) +
                geom_boxplot() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position="none") +
                ylab(expression(paste("Net Shift (",Delta,"pm)"))) +
                ggtitle(paste0(name, " Ch: ", ch, " Control: ", cntl))
        
        allRings <- ggplot(dat.nothermal, 
                           aes(x = factor(Ring), y= `Net Shift`, 
                               fill = Target)) +
                geom_bar(stat = "identity") +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                ylab(expression(paste("Net Shift (",Delta,"pm)"))) + 
                xlab("Ring") +
                ggtitle(paste0(name, " Ch: ", ch, " Control: ", cntl))
        
        # save plot, uncomment to save
        ggsave(plot = plots, 
               file = paste0(loc, "/",  name, "_NetShift_", cntl, "cntl_",
                             "ch", ch, "_step", step, ".png"), 
               width = 10, height = 6)
        ggsave(plot = allRings, 
               file = paste0(loc, "/", name, "_IndyRings","_NetShift_", cntl, 
                             "cntl_", "ch", ch, "_step", step, ".png"), 
               width = 12, height = 6)
}

CheckRingQuality <- function(loc = 'plots', time1, time2, varLevel = 100) {
        library(tidyverse)
        
        dat <- read_csv(paste0(loc,"/", name, "_allRings.csv"))
        dat <- subset(dat, Time > time1 & Time < time2)
        
        dat.avg <- dat %>% group_by(Ring) %>%
                summarise_at(vars(Shift), funs(var))
        
        ringWinners <- filter(dat.avg, Shift < varLevel) %>% select(Ring)
        ringLosers <- filter(dat.avg, Shift > varLevel) %>% select(Ring)
        
        write_csv(ringWinners, paste0(loc, '/', name, "_ringWinners.csv"))
        write_csv(ringLosers, paste0(loc, '/', name, "_ringLosers.csv"))
}

AnalyzeData <- function() {
        GetName()
        AggData(filename = "groupNames_XPP.csv")
        SubtractControl(ch = 1, cntl = "thermal")
        SubtractControl(ch = 2, cntl = "thermal")
        SubtractControl(ch = "U", cntl = "thermal")
        PlotRingData(cntl = "raw", ch = "U", splitPlot = TRUE)
        PlotRingData(cntl = "thermal", ch = "U", splitPlot = TRUE)
        PlotRingData(cntl = "thermal", ch = 1, splitPlot = FALSE)
        PlotRingData(cntl = "raw", ch = 1, splitPlot = FALSE)
        PlotRingData(cntl = "thermal", ch = 2, splitPlot = FALSE)
        PlotRingData(cntl = "raw", ch = 2, splitPlot = FALSE)
        GetNetShifts(cntl = "thermal", ch = 1, time1 = 51, time2 = 39, step = 1)
        GetNetShifts(cntl = "thermal", ch = 2, time1 = 51, time2 = 39, step = 1)
        PlotNetShifts(cntl = "thermal", ch = 1, step = 1)
        PlotNetShifts(cntl = "thermal", ch = 2, step = 1)
        CheckRingQuality(time1 = 10, time2 = 20)
}

AnalyzeAllData <- function() {
        foldersList <- list.dirs(recursive = FALSE)
        lapply(foldersList, function(i){
                directory <- getwd()
                setwd(i)
                AnalyzeData()
                setwd(directory)
        })
}
