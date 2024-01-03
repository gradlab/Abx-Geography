# =============================================================================
# Import
# =============================================================================

library(tidyverse) 
library(broom)
library(viridis)
library(lubridate) 
library(tidycensus)
library(sf)
library(scales) 
library(rmapshaper)
library(tigris)
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)
source('code/utils.R')
source('code/utils_geography.R')
source('code/utils_private.R')
census_api_key(censuskey)
msadat <- setDT(read_csv("data/msadat.csv"))
msadat$MSA <- as.character(msadat$MSA)

figwidth <- 3.2
figaspect <- 7/10
figres <- 600

geotable17 <- read_csv("output/GeoTable17_2022-10-05.csv")

nmemb17 <- geotable17 %>% 
	group_by(STATE, MSA, SEX, AGEGRP) %>% 
	summarise(NMEMB=first(NMEMB))

nmemb17_state <- nmemb17 %>% 
	group_by(STATE) %>% 
	summarise(NMEMB=sum(NMEMB))

nmemb17_MSA <- nmemb17 %>% 
	group_by(MSA) %>% 
	summarise(NMEMB=sum(NMEMB))

# Maybe restrict to MSAs with at least 5000 people and to MSA/Conditions with at least 100 visits? 

geotable17 %>% 
	group_by(PRIMARYCOND, MSA) %>% 
	summarise(NRX=sum(NRX),NVISITS=sum(NVISITS)) %>% 
	split(.$PRIMARYCOND) %>% 
	map(~ filter(., NVISITS>100))

load("data/state_boundaries.RData") 

HHS_boundaries_lwr48 <- state_boundaries %>% 
	make_lwr48 %>%
	makeHHS %>%
	group_by(HHS) %>%
	summarise(geometry=st_union(geometry))

metros <- core_based_statistical_areas(cb = TRUE) %>%
  select(MSA=CBSAFP, metro_name = NAME)

# Read vaccination data: 
vaxdat <- read_csv("data/Vaccination_Coverage_among_Young_Children__0-35_Months_.csv") %>% 
	filter(Geography %in% c("Alabama",
		"Alaska",
		"Arizona",
		"Arkansas",
		"California",
		"Colorado",
		"Connecticut",
		"Delaware",
		"District of Columbia",
		"Florida",
		"Georgia",
		"Hawaii",
		"Idaho",
		"Illinois",
		"Indiana",
		"Iowa",
		"Kansas",
		"Kentucky",
		"Louisiana",
		"Maine",
		"Maryland",
		"Massachusetts",
		"Michigan",
		"Minnesota",
		"Mississippi",
		"Missouri",
		"Montana",
		"Nebraska",
		"Nevada",
		"New Hampshire",
		"New Jersey",
		"New Mexico",
		"New York",
		"North Carolina",
		"North Dakota",
		"Ohio",
		"Oklahoma",
		"Oregon",
		"Pennsylvania",
		"Rhode Island",
		"South Carolina",
		"South Dakota",
		"Tennessee",
		"Texas",
		"Utah",
		"Vermont",
		"Virginia",
		"Washington",
		"West Virginia",
		"Wisconsin",
		"Wyoming")) %>% 
	filter(Dimension=="Overall") %>%
	select(-`Geography Type`, -`Dimension Type`, -Dimension)

vaxdat_formatted <- vaxdat %>% 
	mutate(VAXTYPE=paste(Vaccine,Dose)) %>% 
	select(VAXTYPE, STATE=Geography, VAXCOV=`Estimate (%)`) %>% 
	filter(VAXTYPE %in% c(
		"MMR ≥1 Dose",
		"DTaP ≥3 Doses",
		"PCV ≥3 Doses",
		"Hib Full Series",
		"Influenza NA")) 



# =============================================================================
# Initial mapping
# =============================================================================

df_rxpkp17 <- get_rxpkp(geotable17)
fig_rxpkp17 <- plot_rxpkp_state(df_rxpkp17)
ggsave(fig_rxpkp17, file="figures/rxpkp17_state.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_rxpkp17, file="figures/rxpkp17_state.png", width=figwidth, height=figwidth*figaspect)

df_vpkp17_sinusitis <- get_vpkp(geotable17,cond="Sinusitis")
fig_vpkp17_sinusitis <- plot_vpkp_state(df_vpkp17_sinusitis)
ggsave(fig_vpkp17_sinusitis, file="figures/vpkp17_sinusitis_state.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_vpkp17_sinusitis, file="figures/vpkp17_sinusitis_state.png", width=figwidth, height=figwidth*figaspect)

df_prx17_sinusitis <- get_prx(geotable17,cond="Sinusitis")
fig_prx17_sinusitis <- plot_prx_state(df_prx17_sinusitis)
ggsave(fig_prx17_sinusitis, file="figures/prx17_sinusitis_state.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_prx17_sinusitis, file="figures/prx17_sinusitis_state.png", width=figwidth, height=figwidth*figaspect)

geotable17 %>% 
	mutate(NRX_UNLINKED=case_when(PRIMARYCOND=="Unlinked"~NRX,TRUE~0)) %>% 
	summarise(
		NRX=sum(NRX,na.rm=TRUE), 
		NRX_UNLINKED=sum(NRX_UNLINKED,na.rm=TRUE)) %>% 
	mutate(PUNLINKED=NRX_UNLINKED/NRX)

#        NRX NRX_UNLINKED PUNLINKED
#      <dbl>        <dbl>     <dbl>
# 1 11252508      3000627     0.267

# unlinked by age group: 
geotable17 %>% 
	mutate(NRX_UNLINKED=case_when(PRIMARYCOND=="Unlinked"~NRX,TRUE~0)) %>% 
	group_by(AGEGRP) %>% 
	summarise(
		NRX=sum(NRX,na.rm=TRUE), 
		NRX_UNLINKED=sum(NRX_UNLINKED,na.rm=TRUE)) %>% 
	mutate(PUNLINKED=NRX_UNLINKED/NRX)

#    AGEGRP     NRX NRX_UNLINKED PUNLINKED
#    <chr>    <dbl>        <dbl>     <dbl>
#  1 00_04   933364        58700    0.0629
#  2 05_09   667494        67771    0.102 
#  3 10_14   564956        76961    0.136 
#  4 15_19   731294       177277    0.242 
#  5 20_24   702749       180362    0.257 
#  6 25_29   498244       133632    0.268 
#  7 30_34   677840       179673    0.265 
#  8 35_39   796647       209001    0.262 
#  9 40_44   807620       211882    0.262 
# 10 45_49   899166       240599    0.268 
# 11 50_54   960695       262205    0.273 
# 12 55_59  1049166       297559    0.284 
# 13 60_64   886505       255733    0.288

# =============================================================================
# Population summary
# =============================================================================

nmembdf17 <- geotable17 %>% 
	group_by(STATE,MSA,SEX,AGEGRP) %>% 
	slice(1) %>% 
	select(STATE,MSA,SEX,AGEGRP,NMEMB)

nmembdf17_TOTAL <- nmembdf17 %>% 
	ungroup() %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	mutate(CAT="TOTAL") %>% 
	mutate(VAR="TOTAL") %>% 
	select(VAR, CAT, NMEMB)

nmembdf17_HHS <- nmembdf17 %>% 
	makeHHS %>% 
	group_by(HHS) %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	left_join(HHSnames, by="HHS") %>% 
	mutate(CAT=HHSname) %>% 
	mutate(VAR="HHS") %>% 
	select(VAR, CAT, NMEMB)

nmembdf17_SEX <- nmembdf17 %>% 
	group_by(SEX) %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	mutate(SEXNM=case_when(SEX==1~"Male", TRUE~"Female")) %>% 
	mutate(CAT=SEXNM) %>% 
	mutate(VAR="SEX") %>% 
	select(VAR, CAT, NMEMB)

nmembdf17_AGEGRP <- nmembdf17 %>% 
	group_by(AGEGRP) %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	mutate(CAT=AGEGRP) %>% 
	mutate(VAR="AGEGRP") %>% 
	select(VAR, CAT, NMEMB)

membsummary_df <- bind_rows(nmembdf17_TOTAL, nmembdf17_SEX, nmembdf17_HHS, nmembdf17_AGEGRP) %>% 
	group_by(VAR) %>% 
	mutate(PCT=round(100*NMEMB/sum(NMEMB),1)) %>% 
	mutate(NMEMB=as.character(NMEMB), PCT=as.character(PCT)) %>% 
	mutate(NMEMBLAB=paste0(NMEMB," (",PCT,"%)"))

popsize_total <- nmembdf17_TOTAL$NMEMB

# =============================================================================
# Basic counts 
# =============================================================================

abxtable <- geotable17 %>% 
	group_by(PRIMARYCOND) %>% 
	summarise(NRX=sum(NRX), NVISITS=sum(NVISITS)) %>% 
	mutate(NMEMB=sum(nmembdf17$NMEMB)) %>% 
	mutate(PRX=NRX/NVISITS, RXPKP=NRX/NMEMB*1000, VPKP=NVISITS/NMEMB*1000) %>% 
	arrange(desc(RXPKP)) %>% 
	ungroup() %>% 
	mutate(PCTRX=round(RXPKP/sum(RXPKP)*100,2)) %>% 
	select(PRIMARYCOND, RXPKP, PCTRX, VPKP, PRX)

abxtable_total <- abxtable %>% 
	ungroup() %>% 
	summarise(RXPKP=sum(RXPKP), PCTRX=sum(PCTRX), VPKP=sum(VPKP), PRX=RXPKP/VPKP)

# =============================================================================
# Automate contours
# =============================================================================

sort(unique(geotable17$PRIMARYCOND)) 

# plot_contours_msa(geotable17, ccscond="Viral infection")

fig_contours_asthma <- plot_contours_msa(geotable17, ccscond="Asthma", plot.xlim=c(0,150), plot.ylim=c(0,0.25), contour.limits=c(0,50), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,500))
fig_contours_bactinf <- plot_contours_msa(geotable17, ccscond="Bacterial infection (other)", plot.xlim=c(0,150), plot.ylim=c(0,0.75), contour.limits=c(0,100), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,500))
fig_contours_bronchitis <- plot_contours_msa(geotable17, ccscond="Bronchitis (acute)", plot.xlim=c(0,150), plot.ylim=c(0,1), contour.limits=c(0,100), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,500))
fig_contours_hypertension <- plot_contours_msa(geotable17, ccscond="Essential hypertension", plot.xlim=c(0,1250), plot.ylim=c(0,0.1), contour.limits=c(0,200), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,1250))
fig_contours_influenza <- plot_contours_msa(geotable17, ccscond="Influenza", plot.xlim=c(0,100), plot.ylim=c(0,0.5), contour.limits=c(0,100), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,200))
fig_contours_intestinalinfection <- plot_contours_msa(geotable17, ccscond="Intestinal infection", plot.xlim=c(0,250), plot.ylim=c(0,0.5), contour.limits=c(0,100), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,300))
fig_contours_other <- plot_contours_msa(geotable17, ccscond="Other", plot.xlim=c(0,15000), plot.ylim=c(0,0.05), contour.limits=c(0,10000), contour.breaks=100, contour.breaks.minor=20, vpkp.limits=c(0,15000))
fig_contours_otitismedia <- plot_contours_msa(geotable17, ccscond="Otitis media", plot.xlim=c(0,200), plot.ylim=c(0,1), contour.limits=c(0,200), contour.breaks=50, contour.breaks.minor=10, vpkp.limits=c(0,200))
fig_contours_pneumonia <- plot_contours_msa(geotable17, ccscond="Pneumonia", plot.xlim=c(0,100), plot.ylim=c(0,1), contour.limits=c(0,100), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,100))
fig_contours_sinusitis <- plot_contours_msa(geotable17, ccscond="Sinusitis", plot.xlim=c(0,400), plot.ylim=c(0,1), contour.limits=c(0,350), contour.breaks=50, contour.breaks.minor=10, vpkp.limits=c(0,500))
fig_contours_ssti <- plot_contours_msa(geotable17, ccscond="SSTI", plot.xlim=c(0,150), plot.ylim=c(0,0.75), contour.limits=c(0,200), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,200))
fig_contours_strepphar <- plot_contours_msa(geotable17, ccscond="Strep pharyngitis", plot.xlim=c(0,150), plot.ylim=c(0,1), contour.limits=c(0,200), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,200))
fig_contours_tonsillitis <- plot_contours_msa(geotable17, ccscond="Tonsillitis", plot.xlim=c(0,100), plot.ylim=c(0,0.75), contour.limits=c(0,200), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,200))
fig_contours_uri <- plot_contours_msa(geotable17, ccscond="URI (other)", plot.xlim=c(0,500), plot.ylim=c(0,0.75), contour.limits=c(0,500), contour.breaks=50, contour.breaks.minor=10, vpkp.limits=c(0,500))
fig_contours_uti <- plot_contours_msa(geotable17, ccscond="UTI", plot.xlim=c(0,200), plot.ylim=c(0,0.75), contour.limits=c(0,500), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,500))
fig_contours_viralinfection <- plot_contours_msa(geotable17, ccscond="Viral infection", plot.xlim=c(0,200), plot.ylim=c(0,0.15), contour.limits=c(0,250), contour.breaks=5, contour.breaks.minor=1, vpkp.limits=c(0,200))


ggsave(fig_contours_asthma, file="figures/fig_contours_asthma.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_bactinf, file="figures/fig_contours_bactinf.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_bronchitis, file="figures/fig_contours_bronchitis.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_hypertension, file="figures/fig_contours_hypertension.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_influenza, file="figures/fig_contours_influenza.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_intestinalinfection, file="figures/fig_contours_intestinalinfection.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_other, file="figures/fig_contours_other.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_otitismedia, file="figures/fig_contours_otitismedia.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_pneumonia, file="figures/fig_contours_pneumonia.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_sinusitis, file="figures/fig_contours_sinusitis.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_ssti, file="figures/fig_contours_ssti.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_strepphar, file="figures/fig_contours_strepphar.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_tonsillitis, file="figures/fig_contours_tonsillitis.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_uri, file="figures/fig_contours_uri.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_uti, file="figures/fig_contours_uti.pdf", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_viralinfection, file="figures/fig_contours_viralinfection.pdf", width=figwidth, height=figwidth*figaspect)

ggsave(fig_contours_asthma, file="figures/fig_contours_asthma.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_bactinf, file="figures/fig_contours_bactinf.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_bronchitis, file="figures/fig_contours_bronchitis.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_hypertension, file="figures/fig_contours_hypertension.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_influenza, file="figures/fig_contours_influenza.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_intestinalinfection, file="figures/fig_contours_intestinalinfection.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_other, file="figures/fig_contours_other.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_otitismedia, file="figures/fig_contours_otitismedia.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_pneumonia, file="figures/fig_contours_pneumonia.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_sinusitis, file="figures/fig_contours_sinusitis.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_ssti, file="figures/fig_contours_ssti.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_strepphar, file="figures/fig_contours_strepphar.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_tonsillitis, file="figures/fig_contours_tonsillitis.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_uri, file="figures/fig_contours_uri.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_uti, file="figures/fig_contours_uti.png", width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_viralinfection, file="figures/fig_contours_viralinfection.png", width=figwidth, height=figwidth*figaspect)


# Formatted figures for the paper: 

(fig_contours_sinusitis + labs(title=element_blank(), tag="A)")) %>% 
	ggsave(file="figures/fig_contours_sinusitis_formatted.png", width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_sinusitis + labs(title=element_blank(), tag="A)")) %>% 
	ggsave(file="figures/fig_contours_sinusitis_formatted.pdf", width=figwidth, height=figwidth*8/10, dpi=figres)

(fig_contours_uri + labs(title=element_blank(), tag="B)")) %>% 
	ggsave(file="figures/fig_contours_uri_formatted.png", width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_uri + labs(title=element_blank(), tag="B)")) %>% 
	ggsave(file="figures/fig_contours_uri_formatted.pdf", width=figwidth, height=figwidth*8/10, dpi=figres)

(fig_contours_otitismedia + labs(title=element_blank(), tag="C)")) %>% 
	ggsave(file="figures/fig_contours_otitismedia_formatted.png", width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_otitismedia + labs(title=element_blank(), tag="C)")) %>% 
	ggsave(file="figures/fig_contours_otitismedia_formatted.pdf", width=figwidth, height=figwidth*8/10, dpi=figres)

(fig_contours_strepphar + labs(title=element_blank(), tag="D)")) %>% 
	ggsave(file="figures/fig_contours_strepphar_formatted.png", width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_strepphar + labs(title=element_blank(), tag="D)")) %>% 
	ggsave(file="figures/fig_contours_strepphar_formatted.pdf", width=figwidth, height=figwidth*8/10, dpi=figres)

(fig_contours_bronchitis + labs(title=element_blank(), tag="E)")) %>% 
	ggsave(file="figures/fig_contours_bronchitis_formatted.png", width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_bronchitis + labs(title=element_blank(), tag="E)")) %>% 
	ggsave(file="figures/fig_contours_bronchitis_formatted.pdf", width=figwidth, height=figwidth*8/10, dpi=figres)

# =============================================================================
# A retreat to R2
# =============================================================================

regdf <- geotable17 %>% 
	filter(MSA!="00000") %>% 
	group_by(PRIMARYCOND, STATE, MSA) %>% 
	summarise(NRX=sum(NRX), NVISITS=sum(NVISITS)) %>% 
	ungroup() %>% 
	left_join(
		(geotable17 %>% 
			filter(MSA!="00000") %>% 
			group_by(STATE,MSA,SEX,AGEGRP) %>% 
			slice(1) %>% 
			select(STATE,MSA,SEX,AGEGRP,NMEMB) %>% 
			group_by(STATE,MSA) %>% 
			summarise(NMEMB=sum(NMEMB))),
		by=c("STATE","MSA")) %>% 
	mutate(p=NRX/NVISITS, v=NVISITS/NMEMB*1000, n=NRX/NMEMB*1000) %>% 
	select(PRIMARYCOND, STATE, MSA, p,v,n, NMEMB) %>% 
	split(.$PRIMARYCOND) %>% 
	map(~ list(
			tibble(var="p",R2=(summary(lm(n~p, data=.))$r.squared)),
			tibble(var="v",R2=(summary(lm(n~v, data=.))$r.squared))
			)) %>% 
	map(~ bind_rows(.)) %>% 
	bind_rows(.id="PRIMARYCOND")


fig_regdf <- regdf %>% 
	filter(PRIMARYCOND!="Motor vehicle traffic") %>% 
	(function(x){
		orderednames <- x %>% filter(var=="v") %>% arrange(R2) %>% pull(PRIMARYCOND)
		x <- x %>% mutate(PRIMARYCOND=factor(PRIMARYCOND, levels=orderednames))	
		return(x)
		}) %>% 
	ggplot(aes(x=R2, y=PRIMARYCOND)) +
		geom_line(aes(group=PRIMARYCOND), col="black") + 
		geom_point(aes(col=var)) + 
		scale_color_manual(values=c("blue","red"), labels=c("Prescriptions","Visits")) + 
		theme_classic() + 
		scale_x_continuous(limits=c(0,1)) + 
		labs(x="R-squared", y=element_blank()) + 
		theme(legend.position="bottom") 
		

fig_regdf_neworder <- regdf %>% 
	filter(PRIMARYCOND!="Motor vehicle traffic") %>% 
	mutate(PRIMARYCOND=factor(PRIMARYCOND, levels=rev(c("Sinusitis",
		"URI (other)",
		"Otitis media",
		"Strep pharyngitis",
		"Bronchitis (acute)",
		"SSTI",
		"UTI",
		"Pneumonia",
		"Essential hypertension",
		"Bacterial infection (other)",
		"Asthma",
		"Tonsillitis",
		"Viral infection",
		"Influenza",
		"Intestinal infection",
		"Motor vehicle traffic",
		"Other")))) %>% 
	ggplot(aes(x=R2, y=PRIMARYCOND)) +
		geom_line(aes(group=PRIMARYCOND), col="black") + 
		geom_point(aes(col=var)) + 
		scale_color_manual(values=c("blue","red"), labels=c("Prescriptions","Visits")) + 
		theme_classic() + 
		scale_x_continuous(limits=c(0,1)) + 
		labs(x="R-squared", y=element_blank()) + 
		theme(legend.position="bottom") 


fig_regdf_flipped <- regdf %>% 
	filter(PRIMARYCOND!="Motor vehicle traffic") %>% 
	left_join(select(abxtable, PRIMARYCOND, RXPKP), by="PRIMARYCOND") %>% 
	ggplot(aes(x=log(RXPKP), y=R2)) +
		geom_line(aes(group=PRIMARYCOND), col="lightgray") + 
		geom_point(aes(col=var)) + 
		geom_line(aes(col=var), stat="smooth", method="lm") + 
		scale_color_manual(values=c("blue","red"), labels=c("Prescriptions","Visits")) + 
		theme_classic() + 
		scale_y_continuous(limits=c(0,1)) + 
		labs(x="log(Prescrptions per 1000 people)", y="R2") + 
		theme(legend.position="bottom", text=element_text(size=10))

fig_regdf_flipped %>% 
	ggsave(file="figures/fig_regdf_flipped.pdf", width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_regdf_flipped + theme(legend.position="none") + labs(tag="F)")) %>% 
	ggsave(file="figures/fig_regdf_flipped_nolegend.pdf", width=figwidth, height=figwidth*8/10, dpi=figres)

# Try a sigmoidal regression? 
sigmoiddf <- regdf %>% 
	filter(PRIMARYCOND!="Motor vehicle traffic") %>% 
	left_join(select(abxtable, PRIMARYCOND, RXPKP), by="PRIMARYCOND") %>% 
	mutate(logRXPKP=log(RXPKP)) %>%
	mutate(logitR2=qlogis(R2)) %>% 
	split(.$var) %>% 
	map(~ coefficients(lm(logitR2~logRXPKP, data=.))) %>% 
	map(~ tibble(
		logRXPKP=seq(
			from=min(log(abxtable$RXPKP)), 
			to = max(log(abxtable$RXPKP)),
			by = 0.01), 
		intercept=.["(Intercept)"],
		slope=.["logRXPKP"]
		)) %>% 
	map(~ mutate(., logitR2pred=slope*logRXPKP+intercept)) %>% 
	map(~ select(., logRXPKP, logitR2pred)) %>% 
	bind_rows(.id="var") %>% 
	mutate(R2pred=plogis(logitR2pred)) %>% 
	select(logRXPKP, R2pred, var) 

figsigmoid <- regdf %>% 
	filter(PRIMARYCOND!="Motor vehicle traffic") %>% 
	left_join(select(abxtable, PRIMARYCOND, RXPKP), by="PRIMARYCOND") %>% 
	ggplot(aes(x=log(RXPKP), y=R2)) +
		geom_line(aes(group=PRIMARYCOND), col="lightgray") + 
		geom_point(aes(col=var)) + 
		geom_line(data=sigmoiddf, aes(x=logRXPKP,y=R2pred, col=var)) + 
		scale_color_manual(values=c("blue","red"), labels=c("Prescriptions","Visits")) + 
		theme_classic() + 
		scale_y_continuous(limits=c(0,1)) + 
		scale_x_continuous(limits=c(0,5)) + 
		labs(x="log(Prescrptions per 1000 people)", y="R2") + 
		theme(legend.position="bottom")


# =============================================================================
# Regressions 
# =============================================================================

temp <- geotable17 %>% 
	mutate(RXPKP=NRX/NMEMB*1000) %>% 
	group_by(PRIMARYCOND, STATE, MSA) %>% 
	summarise(NRX=sum(NRX), NVISITS=sum(NVISITS), NMEMB=sum(NMEMB)) %>% 
	select(PRIMARYCOND, STATE, MSA, NRX, NVISITS, NMEMB) %>% 
	mutate(RXPKP=NRX/NMEMB*1000, VPKP=NVISITS/NMEMB*1000) %>% 
	left_join(vaxdat_formatted, by="STATE") %>% 
	split(.$VAXTYPE) %>% 
	map(~ split(., .$PRIMARYCOND)) %>% 
	map(~ imap(., ~ (ggplot(data=.x, aes(x=VAXCOV, y=VPKP)) + 
				geom_point() + 
				geom_line(stat="smooth", method="lm") + 
				labs(title=.y) + 
				theme_classic())))

tempfits <- geotable17 %>% 
	filter(PRIMARYCOND!="Other") %>% 
	filter(PRIMARYCOND!="Motor vehicle traffic") %>% 
	mutate(RXPKP=NRX/NMEMB*1000) %>% 
	group_by(PRIMARYCOND, STATE, MSA) %>% 
	summarise(NRX=sum(NRX), NVISITS=sum(NVISITS), NMEMB=sum(NMEMB)) %>% 
	select(PRIMARYCOND, STATE, MSA, NRX, NVISITS, NMEMB) %>% 
	mutate(RXPKP=NRX/NMEMB*1000, VPKP=NVISITS/NMEMB*1000) %>% 
	left_join(vaxdat_formatted, by="STATE") %>% 
	split(.$VAXTYPE) %>% 
	map(~ split(., .$PRIMARYCOND)) %>% 
	map(~ map(., ~ tidy(lm(VPKP~VAXCOV, data=.)))) %>% 
	map(~ bind_rows(., .id="PRIMARYCOND")) %>% 
	bind_rows(.id="VAXTYPE") %>% 
	filter(term!="(Intercept)") 

fig_tempfits <- tempfits %>% 
	ggplot(aes(x=-log(p.value), y=estimate)) + 
			geom_point() + 
			geom_vline(xintercept=-log(.05), col="grey", lty="dashed") + 
			geom_hline(yintercept=1, col="grey", lty="dashed") + 
			geom_hline(yintercept=-1, col="grey", lty="dashed") + 
			theme_classic() 
























