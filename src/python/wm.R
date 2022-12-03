ver_mjr = R.Version()$major
ver_mnr = unlist(strsplit(R.Version()$minor,"\\."))[1]
lpath = sprintf("~/AppData/Local/R/win-library/%s.%s",ver_mjr,ver_mnr)
rversion = R.Version()
platform = Sys.info()['sysname']
print(sprintf("R version %s.%s running on %s",ver_mjr,ver_mnr,platform))
if(platform == "Windows") {
   if(!file.exists(lpath)) {
     dir.create(lpath,recursive=TRUE)
   }
   libpath <- c(lpath)
  .libPaths(libpath)
}

source('custom.R')

load.libraries(c('emmeans','sciplot','ez','psych','reshape','plyr','ggplot2','afex','dplyr','pastecs'),libpath)

source("tmcustom.R") 
source("lwheatmap.R") 
source("WMSP.R") 

df <- read.csv("fxtn-aois.csv")

# ------------------------------------------------------------------------
# main analyses
#smin <- min(df$aoi_span)
#smax <- max(df$aoi_span)
smin <- -6
smax <- 6

# -- picking out individual conditions -----------------------------------
M <- zeroWM(smin,smax)

ddf <- df
#ddf
M_20 <- TransWMatrix(M,data=ddf,
                    StimulusVar="exp_id",
                    SubjectsVar="subj",
                    SpanVar="aoi_span")
M_20 <- M
M_20

en_20 <- TransEntropy(M,data=ddf,
                    StimulusVar="exp_id",
                    SubjectsVar="subj",
                    SpanVar="aoi_span")
en_20 <- TMentrop
en_20

TransWPlot2(transMatrix=M_20,
           plotName="./figs/WM.pdf",
           plotColors=brewer.pal(9,"Oranges"),
           title="WM",
           annColor='black')

M <- zeroWM(smin,smax)
