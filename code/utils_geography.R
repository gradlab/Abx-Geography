library(tidyverse) 
library(lubridate) 
library(tidycensus)
library(sf)
library(scales) 
library(rmapshaper)
source('code/utils.R')
source('code/utils_private.R')
census_api_key(censuskey)

state_abbrev_df <- read_csv("data/state_abbrev.csv")

us_states <- get_acs(geography = "state", variables = "B01003_001", geometry = TRUE) %>% 
	filter(!(NAME %in% c("Alaska","Hawaii","Puerto Rico"))) %>% 
	mutate(NAME=case_when(NAME=="District of Columbia"~"Washington DC", TRUE~NAME)) %>% 
	select(NAME,geometry) %>% 
	rename(STATE=NAME)

us_states <- ms_simplify(us_states, keep=0.10) # Default keep = 0.10
# ggplot(us_states, aes(geometry=geometry)) + 
# 	geom_sf() + 
# 	theme_minimal() 

get_rxpkp <- function(geotable, cond="Overall", grouper=c("STATE")){

	# Filter by condition:
	if(cond=="Overall"){
		geotable <- geotable %>% 
			group_by(STATE,MSA,SEX,AGEGRP) %>% 
			summarise(NRX=sum(NRX,na.rm=TRUE), NMEMB=first(NMEMB))
	} else {
		geotable <- geotable %>% 
			filter(PRIMARYCOND==cond)
	}

	# Summarize by grouper:
	geotable <- geotable %>% 
		ungroup() %>% 
		group_by(across(all_of(grouper))) %>% 
		summarise(NRX=sum(NRX,na.rm=TRUE), NMEMB=sum(NMEMB)) %>% 
		mutate(RXPKP=NRX/NMEMB*1000) %>% 
		mutate(PRIMARYCOND=cond) %>% 
		select(PRIMARYCOND, STATE, NRX, NMEMB, RXPKP)

	return(geotable)
}

plot_rxpkp_state <- function(df_rxpkp){
	fig_rxpkp <- df_rxpkp %>% 
		filter(!(STATE %in% c("Alaska","Hawaii"))) %>%
		inner_join(us_states, by="STATE") %>% 
		ggplot(aes(geometry=geometry, fill=RXPKP)) + 
			geom_sf(linetype=1,size=0.1) + 
			coord_sf(crs=sf::st_crs(4269)) + 
			scale_fill_distiller(direction=1, palette="Blues") + 
			theme_void() + 
			labs(fill=paste0(df_rxpkp$PRIMARYCOND[1],"\nantibiotic prescriptions\nper 1000 people"))
	return(fig_rxpkp)
}


get_vpkp <- function(geotable, cond="Overall", grouper=c("STATE")){

	# Filter by condition:
	if(cond=="Overall"){
		geotable <- geotable %>% 
			group_by(STATE,MSA,SEX,AGEGRP) %>% 
			summarise(NVISITS=sum(NVISITS,na.rm=TRUE), NMEMB=first(NMEMB))
	} else {
		geotable <- geotable %>% 
			filter(PRIMARYCOND==cond)
	}

	# Summarize by grouper:
	geotable <- geotable %>% 
		ungroup() %>% 
		group_by(across(all_of(grouper))) %>% 
		summarise(NVISITS=sum(NVISITS,na.rm=TRUE), NMEMB=sum(NMEMB)) %>% 
		mutate(VPKP=NVISITS/NMEMB*1000) %>% 
		mutate(PRIMARYCOND=cond) %>% 
		select(PRIMARYCOND, STATE, NVISITS, NMEMB, VPKP)

	return(geotable)
}

plot_vpkp_state <- function(df_vpkp){
	fig_vpkp <- df_vpkp %>% 
		filter(!(STATE %in% c("Alaska","Hawaii"))) %>%
		inner_join(us_states, by="STATE") %>% 
		ggplot(aes(geometry=geometry, fill=VPKP)) + 
			geom_sf(linetype=1,size=0.1) + 
			coord_sf(crs=sf::st_crs(4269)) + 
			scale_fill_distiller(direction=1, palette="Blues") + 
			theme_void() + 
			labs(fill=paste0(df_vpkp$PRIMARYCOND[1]," visits\nper 1000 people"))
	return(fig_vpkp)
}


get_prx <- function(geotable, cond="Overall", grouper=c("STATE")){

	# Filter by condition:
	if(cond=="Overall"){
		geotable <- geotable %>% 
			group_by(STATE,MSA,SEX,AGEGRP) %>% 
			summarise(NVISITS=sum(NVISITS,na.rm=TRUE), NMEMB=first(NMEMB))
	} else {
		geotable <- geotable %>% 
			filter(PRIMARYCOND==cond)
	}

	# Summarize by grouper:
	geotable <- geotable %>% 
		ungroup() %>% 
		group_by(across(all_of(grouper))) %>% 
		summarise(NRX=sum(NRX,na.rm=TRUE), NVISITS=sum(NVISITS,na.rm=TRUE)) %>% 
		mutate(PRX=NRX/NVISITS) %>% 
		mutate(PRIMARYCOND=cond) %>% 
		select(PRIMARYCOND, STATE, NRX, NVISITS, PRX)

	return(geotable)
}

plot_prx_state <- function(df_prx){
	fig_prx <- df_prx %>% 
		filter(!(STATE %in% c("Alaska","Hawaii"))) %>%
		inner_join(us_states, by="STATE") %>% 
		ggplot(aes(geometry=geometry, fill=PRX)) + 
			geom_sf(linetype=1,size=0.1) + 
			coord_sf(crs=sf::st_crs(4269)) + 
			scale_fill_distiller(direction=1, palette="Blues") + 
			theme_void() + 
			labs(fill=paste0(df_prx$PRIMARYCOND[1],"\nprescriptions per visit"))
	return(fig_prx)
}

make_state_abbrev <- function(df){
	df <- df %>% left_join(state_abbrev_df, by="STATE")
	return(df)
}

filterbottom <- function(df, pct){
	cutval <- quantile(df$NMEMB, pct)
	df <- filter(df, NMEMB>cutval)
	return(df)
}

plot_contours_msa <- function(df, ccscond, plot.xlim=c(0,NA), plot.ylim=c(0,1), contour.limits=c(0,350), contour.breaks=50, contour.breaks.minor=10, vpkp.limits=c(0,500), droppct=0, regline=FALSE, incllegend=FALSE, statelabs=FALSE){

	contourdf <- tibble(RXPKP=seq(from=contour.limits[1],to=contour.limits[2],by=contour.breaks)) %>% 
		split(.$RXPKP) %>% 
		map(~ tibble(
			RXPKP=.[[1,1]], 
			VPKP=seq(from=vpkp.limits[1],to=vpkp.limits[2],by=(vpkp.limits[2]-vpkp.limits[1])/100))) %>% 
		map(~ mutate(., PRX=RXPKP/VPKP)) %>% 
		bind_rows()

	contourdf_minor <- tibble(RXPKP=seq(from=contour.limits[1],to=contour.limits[2],by=contour.breaks.minor)) %>% 
		split(.$RXPKP) %>% 
		map(~ tibble(
			RXPKP=.[[1,1]], 
			VPKP=seq(from=vpkp.limits[1],to=vpkp.limits[2],by=(vpkp.limits[2]-vpkp.limits[1])/100))) %>% 
		map(~ mutate(., PRX=RXPKP/VPKP)) %>% 
		bind_rows()

	out <- df %>% 
		filter(PRIMARYCOND==ccscond) %>% 
		makeHHS() %>% 
		mutate(REGION=case_when(
			HHS%in%c(1,2,3)~"Northeast",
			HHS%in%c(4,6)~"South",
			HHS%in%c(5,7,8)~"Mid/Mountain West",
			HHS%in%c(9,10)~"Pacific West"
			)) %>% 
		left_join(state_abbrev_df,by="STATE") %>% 
		group_by(MSA) %>% 
		summarise(NRX=sum(NRX), NVISITS=sum(NVISITS), NMEMB=sum(NMEMB), REGION=first(REGION), ABBREV=first(ABBREV)) %>% 
		filter(MSA!="00000") %>% 
		filterbottom(droppct) %>% 
		mutate(PRX=NRX/NVISITS, VPKP=NVISITS/NMEMB*1000, RXPKP=NRX/NMEMB*1000) %>% 
		ggplot() + 
			geom_point(aes(x=VPKP,y=PRX,size=NMEMB/1000,col=REGION), alpha=0.8)

		if(statelabs==TRUE){
			out <- out + geom_text(aes(x=VPKP,y=PRX,label=ABBREV), size=2)
		}

		out <- out + 
			geom_line(data=contourdf_minor, aes(x=VPKP, y=PRX, fill=factor(RXPKP)), col="lightgray", size=0.2) + 
			geom_line(data=contourdf, aes(x=VPKP, y=PRX, fill=factor(RXPKP)), col="darkgray") + 
			geom_segment(x=0, xend=plot.xlim[2], y=0, yend=0, col="darkgray") + 
			geom_segment(x=0, xend=0, y=0, yend=Inf, col="darkgray") + 
			theme_classic() + 
			theme(text=element_text(size=10)) + 
			# coord_cartesian(xlim=c(50,400), ylim=c(0.2, 1)) + 
			coord_cartesian(xlim=c(plot.xlim[1],plot.xlim[2]), ylim=c(plot.ylim[1], plot.ylim[2])) + 
			labs(title=ccscond, x="Visits per 1000 people", y="Prescriptions per visit", col="Region", size="Number of\nenrollees (thousands)") + 
			scale_color_brewer(type="qual", palette=3)

	# if(statelabs==TRUE){
	# 	out <- out + text(aes(label=ABBREV))
	# }

	if(regline==TRUE){
		out <- out + stat_smooth(aes(x=VPKP, y=PRX), method="loess", span=1, size=0.5, col="darkgrey")
	}

	if(incllegend==FALSE){
		out <- out + theme(legend.position="none")
	}

	return(out)

}


















