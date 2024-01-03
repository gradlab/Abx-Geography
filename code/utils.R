library(tidyverse) 
library(data.table) 

# Save a list of ggplot figures:
save_figlist <- function(obj,dir,name,driver,width=8,height=5){
	mapply(function(x,y) ggsave(x,
		file=paste0(dir,name,"_",y,".",driver), width=width, height=height),
	obj,1:length(obj))
}

tempsave <- function(fig, width=5.6, height=3.5){
	ggsave(fig,file=paste0("/Users/sk792/DropboxHarvard/LabNotebook/figures/",deparse(substitute(fig)),"_",format(Sys.time(), "%Y-%m-%d_%H_%M_%S"),".png"), width=width, height=height)
}

# Prints the length of a data frame: 
printlength <- function(x, msg){
	print(paste(msg,nrow(x),'rows'))
}

# Prints the current time with a message:
printtime <- function(msg){
	print(paste(msg, Sys.time()))
}

# Prints the difference in two timestamps:
printtimediff <- function(msg, startq, endq){
	print(
	paste(msg,
		round(difftime(endq, startq, units='mins'), digits=2),
		' mins')
	)
}

shuffledf <- function(df){
	out <- df[sample(nrow(df)),]
}

# Generate a YMD column from DT_YEAR and DT_MONTH (good for plotting):
makeYMD <- function(x){
	x %>%
	mutate(dummyday=1) %>%
	mutate(YMD=ymd(paste(DT_YEAR,"/",DT_MONTH,"/",dummyday,sep=""))) %>%
	select(-dummyday)
}

trim_dates <- function(x){
	x %>% 
	mutate(keep=case_when(
		DT_YEAR >= 2009 & DT_YEAR <= 2013 ~ 1,
		DT_YEAR == 2008 & DT_MONTH >= 7 ~ 1,
		DT_YEAR == 2014 & DT_MONTH <= 6 ~ 1, 
		TRUE ~ 0
		)) %>%
	filter(keep==1) %>%
	select(-keep)
}

makeHHS <- function(x){
	x %>% 
		mutate(HHS=case_when(
			STATE %in% c("Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont") ~ 1,
			STATE %in% c("New Jersey","New York") ~ 2,
			STATE %in% c("Delaware","Washington DC","Maryland","Pennsylvania","Virginia","West Virginia") ~ 3,
			STATE %in% c("Alabama","Florida","Georgia","Kentucky","Mississippi","North Carolina","South Carolina","Tennessee") ~ 4,
			STATE %in% c("Illinois","Indiana","Michigan","Minnesota","Ohio","Wisconsin") ~ 5,
			STATE %in% c("Arkansas","Louisiana","New Mexico","Oklahoma","Texas") ~ 6,
			STATE %in% c("Iowa","Kansas","Missouri","Nebraska") ~ 7,
			STATE %in% c("Colorado","Montana","North Dakota","South Dakota","Utah","Wyoming") ~ 8,
			STATE %in% c("Arizona","California","Hawaii","Nevada") ~ 9,
			STATE %in% c("Alaska","Idaho","Oregon","Washington") ~ 10,
			TRUE ~ NA_real_
			))
}



calc_gini <- function(r,p){
	gini_df <- tibble(rank=r, prop=p)
	out <- gini_df %>% 
		arrange(rank) %>% 
		mutate(rank=rank/max(rank)) %>% 
		mutate(rankdiff=rank-lag(rank)) %>% 
		mutate(proptrap=(prop+lag(prop))/2) %>% 
		filter(!is.na(rankdiff)) %>% 
		mutate(area_rect=prop*rankdiff, area_trap=proptrap*rankdiff) %>% 
		summarise(
			gini_coef_rect=1-((1-sum(area_rect))/0.5),
			gini_coef_trap=1-((1-sum(area_trap))/0.5)) 
	return(out)
}


HHSnames <- tibble(HHS=1:10, HHSname=c("Boston","New York","Philadelphia","Atlanta","Chicago","Dallas","Kansas City","Denver","San Francisco","Seattle"))

# Reduce to the lower 48 states: 
make_lwr48 <- function(x){
	x %>% 
		filter(!(STATE %in% c("Hawaii","Alaska")))
}

summ_overall <- function(x){
	x %>% 
		group_by(DT_YEAR, DT_MONTH) %>% 
		summarise(NRX=sum(NRX), NMEMB=sum(NMEMB), RXPKP=NRX/NMEMB*1000)
}

summ_by_state <- function(x){
	x %>% 
		group_by(DT_YEAR, DT_MONTH, STATE) %>% 
		summarise(NRX=sum(NRX), NMEMB=sum(NMEMB), RXPKP=NRX/NMEMB*1000)
}

summ_by_age <- function(x){
	x %>% 
		group_by(DT_YEAR, DT_MONTH, AGEGRP) %>% 
		summarise(NRX=sum(NRX), NMEMB=sum(NMEMB), RXPKP=NRX/NMEMB*1000)
}

summ_by_sex <- function(x){
	x %>% 
		group_by(DT_YEAR, DT_MONTH, SEX) %>% 
		summarise(NRX=sum(NRX), NMEMB=sum(NMEMB), RXPKP=NRX/NMEMB*1000)
}

# Turn DT_YEAR and DT_MONTH into a running integer index (good for regression):
DTtoindex_0 <- data.frame(
  DT_YEAR=unlist(map(2008:2018,rep,12)), 
  DT_MONTH=rep(1:12,11), 
  INDEX=0:(12*11-1))

# Store a vector with the number of days per month: 
dayspermonth <- read.csv('data/dayspermonth.csv', colClasses=c('double','double')) %>% setDT()

# Conversion data frame to find states: 
egeoloclist <- read.csv('data/EGEOLOClist_char.csv', colClasses=c('character','character'), stringsAsFactors=FALSE) %>% setDT()

# Store the pointers to the 't' data frames (membership detail)
filenames_t <- c(
	"/data/markscan_authorized/data/commercial/2008/ccaet082.sas7bdat",
	"/data/markscan_authorized/data/commercial/2009/ccaet091.sas7bdat",
	"/data/markscan_authorized/data/commercial/2010/ccaet101.sas7bdat",
	"/data/markscan_authorized/data/commercial/2011/ccaet111.sas7bdat",
	"/data/markscan_authorized/data/commercial/2012/ccaet121.sas7bdat",
	"/data/markscan_authorized/data/commercial/2013/ccaet131.sas7bdat",
	"/data/markscan_authorized/data/commercial/2014/ccaet141.sas7bdat",
	"/data/markscan_authorized/data/commercial/2015/ccaet151.sas7bdat",
	"/data/markscan_authorized/data/commercial/2016/ccaet161.sas7bdat",
	"/data/markscan_authorized/data/commercial/2017/ccaet171.sas7bdat",
	"/data/markscan_authorized/data/commercial/2018/ccaet181.sas7bdat"
	)

# Store the pointers to the 'd' data frames (prescriptions) 
filenames_d <- c(
	"/data/markscan_authorized/data/commercial/2008/ccaed082.sas7bdat",
	"/data/markscan_authorized/data/commercial/2009/ccaed091.sas7bdat",
	"/data/markscan_authorized/data/commercial/2010/ccaed101.sas7bdat",
	"/data/markscan_authorized/data/commercial/2011/ccaed111.sas7bdat",
	"/data/markscan_authorized/data/commercial/2012/ccaed121.sas7bdat",
	"/data/markscan_authorized/data/commercial/2013/ccaed131.sas7bdat",
	"/data/markscan_authorized/data/commercial/2014/ccaed141.sas7bdat",
	"/data/markscan_authorized/data/commercial/2015/ccaed151.sas7bdat",
	"/data/markscan_authorized/data/commercial/2016/ccaed161.sas7bdat",
	"/data/markscan_authorized/data/commercial/2017/ccaed171.sas7bdat",
	"/data/markscan_authorized/data/commercial/2018/ccaed181.sas7bdat"
	)

# Store the pointers to the 'o' data frames (outpatient services)
filenames_o <- c(
	"/data/markscan_authorized/data/commercial/2008/ccaeo082.sas7bdat",
	"/data/markscan_authorized/data/commercial/2009/ccaeo091.sas7bdat",
	"/data/markscan_authorized/data/commercial/2010/ccaeo101.sas7bdat",
	"/data/markscan_authorized/data/commercial/2011/ccaeo111.sas7bdat",
	"/data/markscan_authorized/data/commercial/2012/ccaeo121.sas7bdat",
	"/data/markscan_authorized/data/commercial/2013/ccaeo131.sas7bdat",
	"/data/markscan_authorized/data/commercial/2014/ccaeo141.sas7bdat",
	"/data/markscan_authorized/data/commercial/2015/ccaeo151.sas7bdat",
	"/data/markscan_authorized/data/commercial/2016/ccaeo161.sas7bdat",
	"/data/markscan_authorized/data/commercial/2017/ccaeo171.sas7bdat",
	"/data/markscan_authorized/data/commercial/2018/ccaeo181.sas7bdat"
	)

# Store the pointers to the 's' data frames (inpatient services)
filenames_s <- c(
	"/data/markscan_authorized/data/commercial/2008/ccaes082.sas7bdat",
	"/data/markscan_authorized/data/commercial/2009/ccaes091.sas7bdat",
	"/data/markscan_authorized/data/commercial/2010/ccaes101.sas7bdat",
	"/data/markscan_authorized/data/commercial/2011/ccaes111.sas7bdat",
	"/data/markscan_authorized/data/commercial/2012/ccaes121.sas7bdat",
	"/data/markscan_authorized/data/commercial/2013/ccaes131.sas7bdat",
	"/data/markscan_authorized/data/commercial/2014/ccaes141.sas7bdat",
	"/data/markscan_authorized/data/commercial/2015/ccaes151.sas7bdat",
	"/data/markscan_authorized/data/commercial/2016/ccaes161.sas7bdat",
	"/data/markscan_authorized/data/commercial/2017/ccaes171.sas7bdat",
	"/data/markscan_authorized/data/commercial/2018/ccaes181.sas7bdat"
	)

# Store the pointers to the sampled 't' data frames (membership detail)
filenames_t_sam <- c(
	"/data/markscan_authorized/data/commercial/2008/ccaet082sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2009/ccaet091sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2010/ccaet101sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2011/ccaet111sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2012/ccaet121sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2013/ccaet131sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2014/ccaet141sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2015/ccaet151sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2016/ccaet161sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2017/ccaet171sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2018/ccaet181sam.sas7bdat"
	)

# Store the pointers to the sampled 'd' data frames (prescriptions) 
filenames_d_sam <- c(
	"/data/markscan_authorized/data/commercial/2008/ccaed082sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2009/ccaed091sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2010/ccaed101sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2011/ccaed111sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2012/ccaed121sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2013/ccaed131sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2014/ccaed141sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2015/ccaed151sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2016/ccaed161sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2017/ccaed171sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2018/ccaed181sam.sas7bdat"
	)

# Store the pointers to the sampled 'o' data frames (outpatient services)
filenames_o_sam <- c(
	"/data/markscan_authorized/data/commercial/2008/ccaeo082sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2009/ccaeo091sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2010/ccaeo101sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2011/ccaeo111sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2012/ccaeo121sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2013/ccaeo131sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2014/ccaeo141sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2015/ccaeo151sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2016/ccaeo161sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2017/ccaeo171sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2018/ccaeo181sam.sas7bdat"
	)

# Store the pointers to the sampled 's' data frames (inpatient services)
filenames_s_sam <- c(
	"/data/markscan_authorized/data/commercial/2008/ccaes082sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2009/ccaes091sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2010/ccaes101sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2011/ccaes111sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2012/ccaes121sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2013/ccaes131sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2014/ccaes141sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2015/ccaes151sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2016/ccaes161sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2017/ccaes171sam.sas7bdat",
	"/data/markscan_authorized/data/commercial/2018/ccaes181sam.sas7bdat"
	)

# Store a list of the two-digit years included in the analysis: 
yearlist <- c(
	"08",
	"09",
	"10",
	"11",
	"12",
	"13",
	"14",
	"15",
	"16",
	"17",
	"18")

# Store a list of the years included in the analysis, in numeric form: 
yearlist_num_long <- c(
	2008,
	2009,
	2010,
	2011,
	2012,
	2013,
	2014,
	2015,
	2016,
	2017,
	2018)

# Define a data frame for the month indices (1:60 for five years/60 months):
index_df <- setDT(data.frame(INDEX=1:60))

birthcodes_icd9 <- c(
	"V3000",
	"V3001",
	"V301",
	"V302",
	"V3100",
	"V3101",
	"V311",
	"V312",
	"V3200",
	"V3201",
	"V321",
	"V322",
	"V3300",
	"V3301",
	"V331",
	"V332",
	"V3400",
	"V3401",
	"V341",
	"V342",
	"V3500",
	"V3501",
	"V351",
	"V352",
	"V3600",
	"V3601",
	"V361",
	"V362",
	"V3700",
	"V3701",
	"V371",
	"V372",
	"V3900",
	"V3901",
	"V391",
	"V392")

birthcodes_prefix_icd9 <- "V30|V31|V32|V33|V34|V35|V36|V37|V39"
birthcodes_prefix_icd10 <- "Z38"

birthcodes_prefix_list_icd9 <- c("V30","V31","V32","V33","V34","V35","V36","V37","V39")
birthcodes_prefix_list_icd10 <- c("Z38")

# Specify CPT codes for vaccines: 
vaxcodes <- c(
	"90669", #PCV7
	"90670", #PCV13
	"90470", #influenza
	"90630", #influenza
	"90653", #influenza
	"90654", #influenza
	"90655", #influenza
	"90656", #influenza
	"90657", #influenza
	"90658", #influenza
	"90659", #influenza
	"90660", #influenza
	"90661", #influenza
	"90662", #influenza
	"90663", #influenza
	"90664", #influenza
	"90666", #influenza
	"90668", #influenza
	"90672", #influenza
	"90673", #influenza
	"90674", #influenza
	"90682", #influenza
	"90685", #influenza
	"90686", #influenza
	"90687", #influenza
	"90688", #influenza
	"90694", #influenza
	"90724", #influenza
	"90756", #influenza
	"90723", #polio
	"90696", #polio
	"90697", #polio
	"90698", #polio
	"90712", #polio
	"90713") #polio

# Make DT_YEAR (int), DT_MONTH(int), and YM (date) from DATE column:
make_ym <- function(x){
	x %>% 
	mutate(DT_YEAR=year(DATE)) %>%
	mutate(DT_MONTH=month(DATE)) %>% 
	mutate(DT_DUMMY_DAY=1) %>%
	mutate(YM=ymd(paste0(DT_YEAR,"-",DT_MONTH,"-",DT_DUMMY_DAY)))
}


# Some functions for pulling specific drug claasses: 

pull_abx <- function(x){
	x %>% 
	filter(EVENT=="prescription") %>% 
	inner_join(select(antibiotics, NDCNUM), by=c("CODE"="NDCNUM"))
}

pull_adr <- function(x){
	x %>% 
	filter(EVENT=="prescription") %>% 
	inner_join(select(adrenals, NDCNUM), by=c("CODE"="NDCNUM"))
}

make_cumrx_fig <- function(drug_class){

	rx_df_drug <- rx_df[
		CLASS==drug_class][
		,AGE_DAYS:=as.numeric(difftime(DATE,BIRTH_DATE,units="days"))]

	cumrx_df_drug <- setDT(data.frame())
	for(t in seq(from=0, to=30*60, by=30)){
		cumrx_df_drug <- rbind(cumrx_df_drug,
			rx_df_drug[AGE_DAYS<=t,.(NOBS=.N, NRX=sum(WEIGHT_INDIV),AGE_DAYS_ROUNDED=t),by=.(BIRTH_YEAR)])
	}
	cumrx_df_drug <- cumrx_df_drug[cohort_sizes,nomatch=0,on=.(BIRTH_YEAR)]
	cumrx_df_drug[BIRTH_YEAR==2013 & AGE_DAYS_ROUNDED>(30*12*4), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2014 & AGE_DAYS_ROUNDED>(30*12*3), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2015 & AGE_DAYS_ROUNDED>(30*12*2), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2016 & AGE_DAYS_ROUNDED>(30*12*1), NRX:=NA_integer_]

	cumrx_df_drug <- as_tibble(cumrx_df_drug)
	fig_cumrx_drug <- cumrx_df_drug %>% 
		mutate(rawmean=NOBS/NMEMB) %>% 
		mutate(rawlwr=qgamma(alphasig/2,NOBS,1)/NMEMB) %>%
		mutate(rawupr=qgamma(1-alphasig/2,NOBS+1,1)/NMEMB) %>%
		mutate(proplwr=(rawmean-rawlwr)/rawmean) %>% 
		mutate(propupr=(rawupr-rawmean)/rawmean) %>% 
		mutate(lwr=NRX-proplwr*NRX) %>% 
		mutate(upr=NRX+propupr*NRX) %>% 
		# select(BIRTH_YEAR, NOBS, NRX, rawmean, lwr, upr) %>% 
		ggplot(aes(x=AGE_DAYS_ROUNDED/30, y=NRX, col=factor(BIRTH_YEAR),fill=factor(BIRTH_YEAR))) + 
			geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.4) + 
			geom_line() + 
			scale_x_continuous(breaks=seq(from=0, to=60, by=12)) + 
			theme_minimal() + 
			labs(x="Months from birth", y=paste0("Cumulative ",drug_class," prescriptions"), col="Birth Year", fill="Birth Year") + 
			scale_color_brewer(palette="Set1") + 
			scale_fill_brewer(palette="Set1") +
			theme(text=element_text(size=16))

	return(fig_cumrx_drug)
}


make_cumrx_comp_fig <- function(drug_class_1, drug_class_2){

	rx_df_drug <- rx_df[
		CLASS%in%c(drug_class_1,drug_class_2)][
		,AGE_DAYS:=as.numeric(difftime(DATE,BIRTH_DATE,units="days"))]

	cumrx_df_drug <- setDT(data.frame())
	for(t in seq(from=0, to=30*60, by=30)){
		cumrx_df_drug <- rbind(cumrx_df_drug,
			rx_df_drug[AGE_DAYS<=t,.(NOBS=.N, NRX=sum(WEIGHT_INDIV),AGE_DAYS_ROUNDED=t),by=.(BIRTH_YEAR,CLASS)])
	}
	cumrx_df_drug <- cumrx_df_drug[cohort_sizes,nomatch=0,on=.(BIRTH_YEAR)]
	cumrx_df_drug[BIRTH_YEAR==2013 & AGE_DAYS_ROUNDED>(30*12*4), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2014 & AGE_DAYS_ROUNDED>(30*12*3), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2015 & AGE_DAYS_ROUNDED>(30*12*2), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2016 & AGE_DAYS_ROUNDED>(30*12*1), NRX:=NA_integer_]

	cumrx_df_drug <- as_tibble(cumrx_df_drug)
	fig_cumrx_drug <- cumrx_df_drug %>% 
		mutate(rawmean=NOBS/NMEMB) %>% 
		mutate(rawlwr=qgamma(alphasig/2,NOBS,1)/NMEMB) %>%
		mutate(rawupr=qgamma(1-alphasig/2,NOBS+1,1)/NMEMB) %>%
		mutate(proplwr=(rawmean-rawlwr)/rawmean) %>% 
		mutate(propupr=(rawupr-rawmean)/rawmean) %>% 
		mutate(lwr=NRX-proplwr*NRX) %>% 
		mutate(upr=NRX+propupr*NRX) %>% 
		mutate(`Drug class`=factor(CLASS,levels=c(drug_class_1,drug_class_2))) %>% 
		ggplot(aes(x=AGE_DAYS_ROUNDED/30, y=NRX, col=factor(BIRTH_YEAR), fill=factor(BIRTH_YEAR), lty=`Drug class`)) + 
			geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.4) + 
			geom_line() + 
			scale_x_continuous(breaks=seq(from=0, to=60, by=12)) + 
			theme_minimal() + 
			labs(x="Months from birth", y=paste0("Cumulative prescriptions"), col="Birth Year", fill="Birth Year") + 
			scale_color_brewer(palette="Set1") + 
			scale_fill_brewer(palette="Set1")  +
			theme(text=element_text(size=16))

	return(fig_cumrx_drug)
}

make_ttfirstrx_comp_fig <- function(drug_class_1, drug_class_2){

	first_rx_df_drug <- rx_df[
		CLASS%in%c(drug_class_1,drug_class_2)][
		,AGE_DAYS:=as.numeric(difftime(DATE,BIRTH_DATE,units="days"))][
		,.(BIRTH_YEAR=first(BIRTH_YEAR), WEIGHT_INDIV=first(WEIGHT_INDIV), AGE_DAYS=min(AGE_DAYS)), by=.(ENROLID,CLASS)]

	first_cumrx_df_drug <- setDT(data.frame())
	for(t in seq(from=0, to=30*60, by=30)){
		first_cumrx_df_drug <- rbind(first_cumrx_df_drug,
			first_rx_df_drug[AGE_DAYS<=t,.(NOBS=.N,NRX=sum(WEIGHT_INDIV),AGE_DAYS_ROUNDED=t),by=.(BIRTH_YEAR,CLASS)])
	}
	first_cumrx_df_drug <- first_cumrx_df_drug[cohort_sizes,nomatch=0,on=.(BIRTH_YEAR)]
	first_cumrx_df_drug[BIRTH_YEAR==2013 & AGE_DAYS_ROUNDED>(30*12*4), NRX:=NA_integer_]
	first_cumrx_df_drug[BIRTH_YEAR==2014 & AGE_DAYS_ROUNDED>(30*12*3), NRX:=NA_integer_]
	first_cumrx_df_drug[BIRTH_YEAR==2015 & AGE_DAYS_ROUNDED>(30*12*2), NRX:=NA_integer_]
	first_cumrx_df_drug[BIRTH_YEAR==2016 & AGE_DAYS_ROUNDED>(30*12*1), NRX:=NA_integer_]

	first_cumrx_df_drug <- as_tibble(first_cumrx_df_drug)
	fig_first_cumrx_drug <- first_cumrx_df_drug %>% 
		mutate(rawmean=NOBS/NMEMB) %>% 
		mutate(rawlwr=qgamma(alphasig/2,NOBS,1)/NMEMB) %>%
		mutate(rawupr=qgamma(1-alphasig/2,NOBS+1,1)/NMEMB) %>%
		mutate(proplwr=(rawmean-rawlwr)/rawmean) %>% 
		mutate(propupr=(rawupr-rawmean)/rawmean) %>% 
		mutate(lwr=NRX-proplwr*NRX) %>% 
		mutate(upr=NRX+propupr*NRX) %>% 
		mutate(`Drug class`=factor(CLASS,levels=c(drug_class_1,drug_class_2))) %>% 
		ggplot(aes(x=AGE_DAYS_ROUNDED/30, y=NRX, col=factor(BIRTH_YEAR),fill=factor(BIRTH_YEAR), lty=`Drug class`)) + 
			geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.4) + 
			geom_line() + 
			scale_x_continuous(breaks=seq(from=0, to=60, by=12)) + 
			theme_minimal() + 
			labs(x="Months from birth", y=paste0("Proportion who have received a prescription"), col="Birth Year", fill="Birth Year") + 
			scale_color_brewer(palette="Set1") + 
			scale_fill_brewer(palette="Set1") + 
			theme(text=element_text(size=16))

	return(fig_first_cumrx_drug)
}



make_cumrx_HHS_fig <- function(drug_class){

	rx_df_drug <- rx_df[
		CLASS==drug_class][
		,AGE_DAYS:=as.numeric(difftime(DATE,BIRTH_DATE,units="days"))] %>% 
		as_tibble() %>% 
		makeHHS() %>% 
		setDT()

	cumrx_df_drug <- setDT(data.frame())
	for(t in seq(from=0, to=30*60, by=30)){
		cumrx_df_drug <- rbind(cumrx_df_drug,
			rx_df_drug[AGE_DAYS<=t,.(NOBS=.N, NRX=sum(WEIGHT_INDIV),AGE_DAYS_ROUNDED=t),by=.(BIRTH_YEAR,HHS)])
	}

	# Account for splitting by region: 
	cumrx_df_drug <- (newweights %>% 
		as_tibble() %>% 
		makeHHS() %>% 
		group_by(BIRTH_YEAR, HHS) %>% 
		summarise(POPFRAC=sum(POPFRAC)) %>%
		setDT())[cumrx_df_drug,on=.(BIRTH_YEAR,HHS)][
		,NRX:=NRX/POPFRAC]

	cumrx_df_drug <- cumrx_df_drug[cohort_sizes_HHS,nomatch=0,on=.(BIRTH_YEAR,HHS)]
	cumrx_df_drug[BIRTH_YEAR==2013 & AGE_DAYS_ROUNDED>(30*12*4), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2014 & AGE_DAYS_ROUNDED>(30*12*3), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2015 & AGE_DAYS_ROUNDED>(30*12*2), NRX:=NA_integer_]
	cumrx_df_drug[BIRTH_YEAR==2016 & AGE_DAYS_ROUNDED>(30*12*1), NRX:=NA_integer_]

	cumrx_df_drug <- as_tibble(cumrx_df_drug)
	fig_cumrx_drug <- cumrx_df_drug %>% 
		mutate(rawmean=NOBS/NMEMB) %>% 
		mutate(rawlwr=qgamma(alphasig/2,NOBS,1)/NMEMB) %>%
		mutate(rawupr=qgamma(1-alphasig/2,NOBS+1,1)/NMEMB) %>%
		mutate(proplwr=(rawmean-rawlwr)/rawmean) %>% 
		mutate(propupr=(rawupr-rawmean)/rawmean) %>% 
		mutate(lwr=NRX-proplwr*NRX) %>% 
		mutate(upr=NRX+propupr*NRX) %>% 
		ggplot(aes(x=AGE_DAYS_ROUNDED/30, y=NRX, col=factor(BIRTH_YEAR),fill=factor(BIRTH_YEAR))) + 
			geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.4) + 
			geom_line() + 
			scale_x_continuous(breaks=seq(from=0, to=60, by=12)) + 
			theme_minimal() + 
			labs(x="Months from birth", y=paste0("Cumulative ",drug_class," prescriptions"), col="Birth Year", fill="Birth Year") + 
			scale_color_brewer(palette="Set1") + 
			scale_fill_brewer(palette="Set1") + 
			facet_wrap(~factor(HHS), nrow=2)  +
			theme(text=element_text(size=16))

	return(fig_cumrx_drug)
}


make_ttfirstrx_fig <- function(drug_class){

	first_rx_df_drug <- rx_df[
		CLASS==drug_class][
		,AGE_DAYS:=as.numeric(difftime(DATE,BIRTH_DATE,units="days"))][
		,.(BIRTH_YEAR=first(BIRTH_YEAR), WEIGHT_INDIV=first(WEIGHT_INDIV), AGE_DAYS=min(AGE_DAYS)), by=.(ENROLID)]

	first_cumrx_df_drug <- setDT(data.frame())
	for(t in seq(from=0, to=30*60, by=30)){
		first_cumrx_df_drug <- rbind(first_cumrx_df_drug,
			first_rx_df_drug[AGE_DAYS<=t,.(NOBS=.N,NRX=sum(WEIGHT_INDIV),AGE_DAYS_ROUNDED=t),by=.(BIRTH_YEAR)])
	}
	first_cumrx_df_drug <- first_cumrx_df_drug[cohort_sizes,nomatch=0,on=.(BIRTH_YEAR)]
	first_cumrx_df_drug[BIRTH_YEAR==2013 & AGE_DAYS_ROUNDED>(30*12*4), NRX:=NA_integer_]
	first_cumrx_df_drug[BIRTH_YEAR==2014 & AGE_DAYS_ROUNDED>(30*12*3), NRX:=NA_integer_]
	first_cumrx_df_drug[BIRTH_YEAR==2015 & AGE_DAYS_ROUNDED>(30*12*2), NRX:=NA_integer_]
	first_cumrx_df_drug[BIRTH_YEAR==2016 & AGE_DAYS_ROUNDED>(30*12*1), NRX:=NA_integer_]

	first_cumrx_df_drug <- as_tibble(first_cumrx_df_drug)
	fig_first_cumrx_drug <- first_cumrx_df_drug %>% 
		mutate(rawmean=NOBS/NMEMB) %>% 
		mutate(rawlwr=qgamma(alphasig/2,NOBS,1)/NMEMB) %>%
		mutate(rawupr=qgamma(1-alphasig/2,NOBS+1,1)/NMEMB) %>%
		mutate(proplwr=(rawmean-rawlwr)/rawmean) %>% 
		mutate(propupr=(rawupr-rawmean)/rawmean) %>% 
		mutate(lwr=NRX-proplwr*NRX) %>% 
		mutate(upr=NRX+propupr*NRX) %>% 
		ggplot(aes(x=AGE_DAYS_ROUNDED/30, y=NRX, col=factor(BIRTH_YEAR),fill=factor(BIRTH_YEAR))) + 
			geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.4) + 
			geom_line() + 
			scale_x_continuous(breaks=seq(from=0, to=60, by=12)) + 
			theme_minimal() + 
			labs(x="Months from birth", y=paste0("Proportion who have received an ",drug_class," prescription"), col="Birth Year", fill="Birth Year") + 
			scale_color_brewer(palette="Set1") + 
			scale_fill_brewer(palette="Set1") + 
			theme(text=element_text(size=16))

	return(fig_first_cumrx_drug)
}

make_cumvisit_fig <- function(cond_name){

if(is.na(cond_name)){

	visit_df_slice <- visit_df[is.na(COND)][
		,AGE_DAYS:=as.numeric(difftime(DATE,BIRTH_DATE,units="days"))]

} else {

	visit_df_slice <- visit_df[COND==cond_name][
		,AGE_DAYS:=as.numeric(difftime(DATE,BIRTH_DATE,units="days"))]	
}


cumvisit_df_slice <- setDT(data.frame())
for(t in seq(from=0, to=30*60, by=30)){
	cumvisit_df_slice <- rbind(cumvisit_df_slice,
		visit_df_slice[AGE_DAYS<=t,.(NOBS=.N, NVISITS=sum(WEIGHT_INDIV),AGE_DAYS_ROUNDED=t),by=.(BIRTH_YEAR)])
}
cumvisit_df_slice <- cumvisit_df_slice[cohort_sizes,nomatch=0,on=.(BIRTH_YEAR)]
cumvisit_df_slice[BIRTH_YEAR==2013 & AGE_DAYS_ROUNDED>(30*12*4), NVISITS:=NA_integer_]
cumvisit_df_slice[BIRTH_YEAR==2014 & AGE_DAYS_ROUNDED>(30*12*3), NVISITS:=NA_integer_]
cumvisit_df_slice[BIRTH_YEAR==2015 & AGE_DAYS_ROUNDED>(30*12*2), NVISITS:=NA_integer_]
cumvisit_df_slice[BIRTH_YEAR==2016 & AGE_DAYS_ROUNDED>(30*12*1), NVISITS:=NA_integer_]

cumvisit_df_slice <- as_tibble(cumvisit_df_slice)
fig_cumvisit_slice <- cumvisit_df_slice %>% 
	mutate(rawmean=NOBS/NMEMB) %>% 
	mutate(rawlwr=qgamma(alphasig/2,NOBS,1)/NMEMB) %>%
	mutate(rawupr=qgamma(1-alphasig/2,NOBS+1,1)/NMEMB) %>%
	mutate(proplwr=(rawmean-rawlwr)/rawmean) %>% 
	mutate(propupr=(rawupr-rawmean)/rawmean) %>% 
	mutate(lwr=NVISITS-proplwr*NVISITS) %>% 
	mutate(upr=NVISITS+propupr*NVISITS) %>% 
	ggplot(aes(x=AGE_DAYS_ROUNDED/30, y=NVISITS, col=factor(BIRTH_YEAR),fill=factor(BIRTH_YEAR))) + 
		geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.4) + 
		geom_line() + 
		scale_x_continuous(breaks=seq(from=0, to=60, by=12)) + 
		theme_minimal() + 
		labs(x="Months from birth", y=paste0("Cumulative ",cond_name," visits"), col="Birth Year", fill="Birth Year") + 
		scale_color_brewer(palette="Set1") + 
		scale_fill_brewer(palette="Set1")  +
		theme(text=element_text(size=16))

return(fig_cumvisit_slice)

}

# make_cumvisit_fig <- function(cond_name){

# if(is.na(cond_name)){

# 	visit_df_slice <- visit_df[is.na(COND)][
# 		# Inner-join to on 'memb_df' to grab demographic information:
# 		memb_df,nomatch=0,on=.(ENROLID)][
# 		# Extract ENROLID, birth year, and age (in days) for each rx:
# 		,.(ENROLID, 
# 			BIRTH_YEAR=year(BIRTH_DATE), 
# 			AGE_DAYS=as.numeric(difftime(DATE,BIRTH_DATE,units="days")))]

# } else {

# 	visit_df_slice <- visit_df[COND==cond_name][
# 		# Inner-join to on 'memb_df' to grab demographic information:
# 		memb_df,nomatch=0,on=.(ENROLID)][
# 		# Extract ENROLID, birth year, and age (in days) for each rx:
# 		,.(ENROLID, 
# 			BIRTH_YEAR=year(BIRTH_DATE), 
# 			AGE_DAYS=as.numeric(difftime(DATE,BIRTH_DATE,units="days")))]	
# }


# cumvisit_df_slice <- setDT(data.frame())
# for(t in seq(from=0, to=30*60, by=30)){
# 	cumvisit_df_slice <- rbind(cumvisit_df_slice,
# 		visit_df_slice[AGE_DAYS<=t,.(NVISITS=.N,AGE_DAYS_ROUNDED=t),by=.(BIRTH_YEAR)])
# }
# cumvisit_df_slice <- cumvisit_df_slice[cohort_sizes,nomatch=0,on=.(BIRTH_YEAR)]
# cumvisit_df_slice[BIRTH_YEAR==2013 & AGE_DAYS_ROUNDED>(30*12*4), NVISITS:=NA_integer_]
# cumvisit_df_slice[BIRTH_YEAR==2014 & AGE_DAYS_ROUNDED>(30*12*3), NVISITS:=NA_integer_]
# cumvisit_df_slice[BIRTH_YEAR==2015 & AGE_DAYS_ROUNDED>(30*12*2), NVISITS:=NA_integer_]
# cumvisit_df_slice[BIRTH_YEAR==2016 & AGE_DAYS_ROUNDED>(30*12*1), NVISITS:=NA_integer_]

# cumvisit_df_slice <- as_tibble(cumvisit_df_slice)
# fig_cumvisit_slice <- cumvisit_df_slice %>% 
# 	mutate(lwr=qgamma(alphasig/2,NVISITS,1)/NMEMB) %>%
# 	mutate(upr=qgamma(1-alphasig/2,NVISITS+1,1)/NMEMB) %>%
# 	ggplot(aes(x=AGE_DAYS_ROUNDED/30, y=NVISITS/NMEMB, col=factor(BIRTH_YEAR),fill=factor(BIRTH_YEAR))) + 
# 		geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.4) + 
# 		geom_line() + 
# 		scale_x_continuous(breaks=seq(from=0, to=60, by=12)) + 
# 		theme_minimal() + 
# 		labs(x="Months from birth", y=paste0("Cumulative ",cond_name," visits"), col="Birth Year", fill="Birth Year") + 
# 		scale_color_brewer(palette="Set1") + 
# 		scale_fill_brewer(palette="Set1") 

# return(fig_cumvisit_slice)

# # cumvisit_df_slice %>%
# # 	filter(!is.na(NVISITS)) %>% 
# # 	group_by(BIRTH_YEAR) %>% 
# # 	arrange(BIRTH_YEAR, desc(AGE_DAYS_ROUNDED)) %>% 
# # 	slice(1) %>% 
# # 	mutate(NVPP=NVISITS/NMEMB)

# }



make_cumvisit_vax_fig <- function(cond_name){
if(is.na(cond_name)){

	visit_df_slice <- visit_df[is.na(COND)][
	# Inner-join to on 'memb_df' to grab demographic information:
	pcvcohorts,nomatch=0,on=.(ENROLID)][
	# Extract ENROLID, birth year, and age (in days) for each rx:
	,.(ENROLID, 
		BIRTH_YEAR=year(BIRTH_DATE), 
		AGE_DAYS=as.numeric(difftime(DATE,BIRTH_DATE,units="days")),
		PCVCOHORT=PCVCOHORT)]	
	} else {

	visit_df_slice <- visit_df[COND==cond_name][
		# Inner-join to on 'memb_df' to grab demographic information:
	pcvcohorts,nomatch=0,on=.(ENROLID)][
	# Extract ENROLID, birth year, and age (in days) for each rx:
	,.(ENROLID, 
		BIRTH_YEAR=year(BIRTH_DATE), 
		AGE_DAYS=as.numeric(difftime(DATE,BIRTH_DATE,units="days")),
		PCVCOHORT=PCVCOHORT)]	
}


cumvisit_df_slice <- setDT(data.frame())
for(t in seq(from=0, to=30*60, by=30)){
	cumvisit_df_slice <- rbind(cumvisit_df_slice,
		visit_df_slice[AGE_DAYS<=t,.(NVISITS=.N,AGE_DAYS_ROUNDED=t),by=.(BIRTH_YEAR, PCVCOHORT)])
}
cumvisit_df_slice <- cumvisit_df_slice[pcv_cohort_sizes,nomatch=0,on=.(PCVCOHORT, BIRTH_YEAR)]
cumvisit_df_slice[BIRTH_YEAR==2013 & AGE_DAYS_ROUNDED>(30*12*4), NVISITS:=NA_integer_]
cumvisit_df_slice[BIRTH_YEAR==2014 & AGE_DAYS_ROUNDED>(30*12*3), NVISITS:=NA_integer_]
cumvisit_df_slice[BIRTH_YEAR==2015 & AGE_DAYS_ROUNDED>(30*12*2), NVISITS:=NA_integer_]
cumvisit_df_slice[BIRTH_YEAR==2016 & AGE_DAYS_ROUNDED>(30*12*1), NVISITS:=NA_integer_]

cumvisit_df_slice <- as_tibble(cumvisit_df_slice)
fig_cumvisit_slice <- cumvisit_df_slice %>% 
	filter(!(BIRTH_YEAR==2008 & PCVCOHORT=="13")) %>% 
	filter(!(BIRTH_YEAR==2009 & PCVCOHORT=="13")) %>% 
	mutate(lwr=qgamma(alphasig/2,NVISITS,1)/NMEMB) %>%
	mutate(upr=qgamma(1-alphasig/2,NVISITS+1,1)/NMEMB) %>%
	ggplot(aes(x=AGE_DAYS_ROUNDED/30, y=NVISITS/NMEMB, col=factor(BIRTH_YEAR),fill=factor(BIRTH_YEAR))) + 
		geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.4) + 
		geom_line() + 
		scale_x_continuous(breaks=seq(from=0, to=60, by=12)) + 
		theme_minimal() + 
		labs(x="Months from birth", y=paste0("Cumulative ",cond_name," visits"), col="Birth Year", fill="Birth Year") + 
		scale_color_brewer(palette="Set1") + 
		scale_fill_brewer(palette="Set1") + 
		facet_wrap(~PCVCOHORT, nrow=2)

return(fig_cumvisit_slice)

}

clean_carolinas <- function(x){
	x %>% 
		mutate(STATE=case_when(
			STATE=="North Carolin"~"North Carolina",
			STATE=="South Carolin"~"South Carolina",
			TRUE~STATE
			))
}