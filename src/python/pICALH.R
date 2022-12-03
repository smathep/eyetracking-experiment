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

source('customFunc.R')
load.libraries(c('tidyverse','ggplot2','emmeans','afex','corrplot'))
df <- read_csv(file = "pICALH.csv")

st <- read_csv("fxtn-aois.csv") %>% 
	filter(key == TRUE) %>%  
	select(subj, block, shot, frequency, length, duration, 
		prev_sacc_amplitude) %>% 
	group_by(subj, block, shot, frequency) %>% 
	summarize(length = mean(length),
		dur = mean(duration),
		sac = mean(prev_sacc_amplitude)) %>% 
	ungroup() %>% 
	mutate(len = recode(length, `4` = "short", `3` = "short", 
		`10` = "long", `11` = "long"))

d <- inner_join(df, st)


# ==============================================================================
(fit <- aov_ez(d, id = "subj", dv = "pICALH", within = c("frequency", "len","block")))

emmeans(fit, ~block)
emmeans(fit, ~len)
emmeans(fit, ~frequency)

emmip(fit, ~block)
ggsave("graphs/pICALH_block.pdf")
dev.off()

emmip(fit, ~len)
ggsave("graphs/pICALH_length.pdf")
dev.off()

emmip(fit, ~frequency)
ggsave("graphs/pICALH_freq.pdf")
dev.off()

#  ============================================================================
(fit <- aov_ez(d, id = "subj", dv = "dur", within = c("frequency", "len","block")))

emmip(fit, frequency~len)

#  ============================================================================
(fit <- aov_ez(d, id = "subj", dv = "sac", within = c("frequency", "len","block")))

emmip(fit, frequency~len)

# correlations
cd <- d %>% 
  select(pICALH, sac, dur)

M <- cor(cd)
corrplot(M, type = "upper", addCoef.col = "black", diag = FALSE)
ggsave("graphs/corr.pdf")
dev.off()

