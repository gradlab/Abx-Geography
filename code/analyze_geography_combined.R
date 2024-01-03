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

geotable15 <- mutate(read_csv("output/GeoTable15_2023-03-14.csv"),YEAR=2015)
geotable16 <- mutate(read_csv("output/GeoTable16_2023-03-14.csv"),YEAR=2016)
geotable17 <- mutate(read_csv("output/GeoTable17_2023-03-14.csv"),YEAR=2017)
geotable18 <- mutate(read_csv("output/GeoTable18_2023-03-14.csv"),YEAR=2018)

# unlinked:
geotable <- bind_rows(
	filter(geotable15,PRIMARYCOND!="Unlinked"),
	filter(geotable16,PRIMARYCOND!="Unlinked"),
	filter(geotable17,PRIMARYCOND!="Unlinked"),
	filter(geotable18,PRIMARYCOND!="Unlinked"))

# with unlinked
geotable_wunlinked <- bind_rows(geotable15, geotable16, geotable17, geotable18)



metros <- core_based_statistical_areas(cb = TRUE) %>%
  select(MSA=CBSAFP, metro_name = NAME)





# =============================================================================
# Population summary
# =============================================================================

nmembdf <- geotable %>% 
	group_by(STATE,MSA,SEX,AGEGRP,YEAR) %>% 
	slice(1) %>% 
	select(STATE,MSA,SEX,AGEGRP,YEAR,NMEMB)

nmembdf_TOTAL <- nmembdf %>% 
	ungroup() %>% 
	group_by(YEAR) %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	mutate(CAT="TOTAL") %>% 
	mutate(VAR="TOTAL") %>% 
	select(YEAR, VAR, CAT, NMEMB)

nmembdf_MSA <- nmembdf %>% 
	ungroup() %>% 
	group_by(YEAR,MSA) %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	mutate(CAT=MSA) %>% 
	mutate(VAR="MSA") %>% 
	select(YEAR, VAR, CAT, NMEMB)

nmembdf_HHS <- nmembdf %>% 
	makeHHS %>% 
	group_by(YEAR, HHS) %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	left_join(HHSnames, by="HHS") %>% 
	mutate(CAT=HHSname) %>% 
	mutate(VAR="HHS") %>% 
	select(YEAR, VAR, CAT, NMEMB)

nmembdf_SEX <- nmembdf %>% 
	group_by(YEAR, SEX) %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	mutate(SEXNM=case_when(SEX==1~"Male", TRUE~"Female")) %>% 
	mutate(CAT=SEXNM) %>% 
	mutate(VAR="SEX") %>% 
	select(YEAR, VAR, CAT, NMEMB)

nmembdf_AGEGRP <- nmembdf %>% 
	group_by(YEAR, AGEGRP) %>% 
	summarise(NMEMB=sum(NMEMB)) %>% 
	mutate(CAT=AGEGRP) %>% 
	mutate(VAR="AGEGRP") %>% 
	select(YEAR, VAR, CAT, NMEMB)

membsummary_df <- bind_rows(nmembdf_TOTAL, nmembdf_SEX, nmembdf_HHS, nmembdf_AGEGRP) %>% 
	group_by(YEAR, VAR) %>% 
	mutate(PCT=round(100*NMEMB/sum(NMEMB),1)) %>% 
	mutate(NMEMB=as.character(NMEMB), PCT=as.character(PCT)) %>% 
	mutate(NMEMBLAB=paste0(NMEMB," (",PCT,"%)")) %>% 
	select(YEAR,VAR,CAT,NMEMBLAB) %>% 
	pivot_wider(names_from=YEAR, values_from=NMEMBLAB)

write_csv(membsummary_df, file="figures/combined/membsummary.csv")

# popsize_total <- nmembdf_TOTAL$NMEMB

# =============================================================================
# Basic counts 
# =============================================================================

abxtable <- geotable %>% 
	group_by(PRIMARYCOND) %>% 
	summarise(NRX=sum(NRX), NVISITS=sum(NVISITS)) %>% 
	mutate(NMEMB=sum(nmembdf$NMEMB)) %>% 
	mutate(PRX=NRX/NVISITS, RXPKP=NRX/NMEMB*1000, VPKP=NVISITS/NMEMB*1000) %>% 
	arrange(desc(RXPKP)) %>% 
	ungroup() %>% 
	mutate(PCTRX=round(RXPKP/sum(RXPKP)*100,2)) %>% 
	select(PRIMARYCOND, RXPKP, PCTRX, VPKP, PRX)

condlevels <- c((abxtable$PRIMARYCOND)[2:nrow(abxtable)],"Other")
abxtable_formatted <- abxtable %>% 
	mutate(PRIMARYCOND=factor(PRIMARYCOND,levels=condlevels)) %>% 
	arrange(PRIMARYCOND) %>% 
	mutate(RXPKP=round(RXPKP,1)) %>% 
	mutate(PCTRX=round(PCTRX,1)) %>% 
	mutate(VPKP=round(VPKP,1)) %>% 
	mutate(PRX=round(PRX*100,1))

abxtable_total <- abxtable %>% 
	ungroup() %>% 
	summarise(RXPKP=round(sum(RXPKP),1), PCTRX=round(sum(PCTRX),1), VPKP=round(sum(VPKP),1), PRX=round(100*RXPKP/VPKP,1)) %>% 
	mutate(PRIMARYCOND="Total") %>% 
	select(PRIMARYCOND,RXPKP,PCTRX,VPKP,PRX) 

write_csv(bind_rows(abxtable_formatted,abxtable_total), file=paste0("figures/combined/abxtable.csv"))


abxtable_regional <- geotable %>% 
	makeHHS() %>% 
	mutate(REGION=case_when(
		HHS%in%c(1,2,3)~"Northeast",
		HHS%in%c(4,6)~"South",
		HHS%in%c(5,7,8)~"Mid/Mountain West",
		HHS%in%c(9,10)~"Pacific West"
		)) %>% 
	group_by(PRIMARYCOND,REGION) %>% 
	summarise(NRX=sum(NRX), NVISITS=sum(NVISITS)) %>% 
	left_join((nmembdf %>% 
			makeHHS() %>% 
			mutate(REGION=case_when(
				HHS%in%c(1,2,3)~"Northeast",
				HHS%in%c(4,6)~"South",
				HHS%in%c(5,7,8)~"Mid/Mountain West",
				HHS%in%c(9,10)~"Pacific West"
				)) %>% 
			group_by(REGION) %>% 
			summarise(NMEMB=sum(NMEMB))), by="REGION") %>% 
	mutate(NMEMB_NAT=sum(nmembdf$NMEMB)) %>% 
	mutate(PRX=NRX/NVISITS, RXPKP=NRX/NMEMB*1000, VPKP=NVISITS/NMEMB*1000) %>% 
	mutate(RXPKP_NAT=NRX/NMEMB_NAT*1000) %>% 
	ungroup() %>% 
	mutate(PCTRX=round(RXPKP/sum(RXPKP)*100,2)) %>% 
	mutate(PCTRX_NAT=round(RXPKP_NAT/sum(RXPKP_NAT)*100,2)) %>% 
	filter(PRIMARYCOND %in% c("Sinusitis",
		"URI (other)",
		"Otitis media",
		"Strep pharyngitis",
		"Bronchitis (acute)")) %>% 
	mutate(PRIMARYCOND=factor(PRIMARYCOND, levels=c("Sinusitis",
		"URI (other)",
		"Otitis media",
		"Strep pharyngitis",
		"Bronchitis (acute)"))) %>% 
	mutate(REGION=factor(REGION,levels=c("Northeast",
		"South",
		"Mid/Mountain West",
		"Pacific West"
		))) %>% 
	arrange(PRIMARYCOND, REGION) %>% 
	select(PRIMARYCOND, REGION, RXPKP, PCTRX, PCTRX_NAT, VPKP, PRX)

abxtable_regional_formatted <- abxtable_regional %>% 
	mutate(RXPKP=round(RXPKP,1), PCTRX_NAT=round(PCTRX_NAT,1), VPKP=round(VPKP,1), PRX=round(100*PRX,1)) %>% 
	select(PRIMARYCOND, REGION, RXPKP, PCTRX_NAT, VPKP, PRX)

write_csv(abxtable_regional_formatted, file=paste0("figures/combined/abxtable_regional.csv"))



# regional totals
abxtable_regional_total <- geotable %>% 
  makeHHS() %>% 
  mutate(REGION=case_when(
    HHS%in%c(1,2,3)~"Northeast",
    HHS%in%c(4,6)~"South",
    HHS%in%c(5,7,8)~"Mid/Mountain West",
    HHS%in%c(9,10)~"Pacific West"
  )) %>% 
  group_by(REGION) %>% 
  summarise(NRX=sum(NRX), NVISITS=sum(NVISITS)) %>% 
  left_join((nmembdf %>% 
               makeHHS() %>% 
               mutate(REGION=case_when(
                 HHS%in%c(1,2,3)~"Northeast",
                 HHS%in%c(4,6)~"South",
                 HHS%in%c(5,7,8)~"Mid/Mountain West",
                 HHS%in%c(9,10)~"Pacific West"
               )) %>% 
               group_by(REGION) %>% 
               summarise(NMEMB=sum(NMEMB))), by="REGION") %>% 
  mutate(NMEMB_NAT=sum(nmembdf$NMEMB)) %>% 
  mutate(PRX=NRX/NVISITS, RXPKP=NRX/NMEMB*1000, VPKP=NVISITS/NMEMB*1000) %>% 
  mutate(RXPKP_NAT=NRX/NMEMB_NAT*1000) %>% 
  ungroup() %>% 
  mutate(PCTRX=round(RXPKP/sum(RXPKP)*100,2)) %>% 
  mutate(PCTRX_NAT=round(RXPKP_NAT/sum(RXPKP_NAT)*100,2))

# =============================================================================
# Automate contours
# =============================================================================

sort(unique(geotable$PRIMARYCOND)) 

# plot_contours_msa(geotable, ccscond="Viral infection")

fig_contours_asthma_legend <- plot_contours_msa(geotable, ccscond="Asthma", plot.xlim=c(0,150), plot.ylim=c(0,0.25), contour.limits=c(0,50), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,500),incllegend=TRUE)
fig_contours_asthma_statelabs <- plot_contours_msa(geotable, ccscond="Asthma", plot.xlim=c(0,150), plot.ylim=c(0,0.25), contour.limits=c(0,50), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,500),statelabs=TRUE)
fig_contours_asthma <- plot_contours_msa(geotable, ccscond="Asthma", plot.xlim=c(0,150), plot.ylim=c(0,0.25), contour.limits=c(0,50), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,500))
fig_contours_bactinf <- plot_contours_msa(geotable, ccscond="Bacterial infection (other)", plot.xlim=c(0,150), plot.ylim=c(0,0.75), contour.limits=c(0,100), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,500))
fig_contours_bronchitis <- plot_contours_msa(geotable, ccscond="Bronchitis (acute)", plot.xlim=c(0,150), plot.ylim=c(0,1), contour.limits=c(0,100), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,500))
fig_contours_hypertension <- plot_contours_msa(geotable, ccscond="Essential hypertension", plot.xlim=c(0,1250), plot.ylim=c(0,0.1), contour.limits=c(0,200), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,1250))
fig_contours_influenza <- plot_contours_msa(geotable, ccscond="Influenza", plot.xlim=c(0,100), plot.ylim=c(0,0.5), contour.limits=c(0,100), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,200))
fig_contours_intestinalinfection <- plot_contours_msa(geotable, ccscond="Intestinal infection", plot.xlim=c(0,250), plot.ylim=c(0,0.5), contour.limits=c(0,100), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,300))
fig_contours_other <- plot_contours_msa(geotable, ccscond="Other", plot.xlim=c(0,15000), plot.ylim=c(0,0.05), contour.limits=c(0,10000), contour.breaks=100, contour.breaks.minor=20, vpkp.limits=c(0,15000))
fig_contours_otitismedia <- plot_contours_msa(geotable, ccscond="Otitis media", plot.xlim=c(0,200), plot.ylim=c(0,1), contour.limits=c(0,200), contour.breaks=50, contour.breaks.minor=10, vpkp.limits=c(0,200))
fig_contours_pneumonia <- plot_contours_msa(geotable, ccscond="Pneumonia", plot.xlim=c(0,100), plot.ylim=c(0,1), contour.limits=c(0,100), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,100))
fig_contours_sinusitis <- plot_contours_msa(geotable, ccscond="Sinusitis", plot.xlim=c(0,400), plot.ylim=c(0,1), contour.limits=c(0,350), contour.breaks=50, contour.breaks.minor=10, vpkp.limits=c(0,500))
fig_contours_ssti <- plot_contours_msa(geotable, ccscond="SSTI", plot.xlim=c(0,150), plot.ylim=c(0,0.75), contour.limits=c(0,200), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,200))
fig_contours_strepphar <- plot_contours_msa(geotable, ccscond="Strep pharyngitis", plot.xlim=c(0,150), plot.ylim=c(0,1), contour.limits=c(0,200), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,200))
fig_contours_tonsillitis <- plot_contours_msa(geotable, ccscond="Tonsillitis", plot.xlim=c(0,100), plot.ylim=c(0,0.75), contour.limits=c(0,200), contour.breaks=10, contour.breaks.minor=2, vpkp.limits=c(0,200))
fig_contours_uri <- plot_contours_msa(geotable, ccscond="URI (other)", plot.xlim=c(0,500), plot.ylim=c(0,0.75), contour.limits=c(0,500), contour.breaks=50, contour.breaks.minor=10, vpkp.limits=c(0,500))
fig_contours_uti <- plot_contours_msa(geotable, ccscond="UTI", plot.xlim=c(0,200), plot.ylim=c(0,0.75), contour.limits=c(0,500), contour.breaks=20, contour.breaks.minor=5, vpkp.limits=c(0,500))
fig_contours_viralinfection <- plot_contours_msa(geotable, ccscond="Viral infection", plot.xlim=c(0,200), plot.ylim=c(0,0.15), contour.limits=c(0,250), contour.breaks=5, contour.breaks.minor=1, vpkp.limits=c(0,200))


ggsave(fig_contours_asthma, file=paste0("figures/combined/fig_contours_asthma.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_bactinf, file=paste0("figures/combined/fig_contours_bactinf.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_bronchitis, file=paste0("figures/combined/fig_contours_bronchitis.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_hypertension, file=paste0("figures/combined/fig_contours_hypertension.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_influenza, file=paste0("figures/combined/fig_contours_influenza.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_intestinalinfection, file=paste0("figures/combined/fig_contours_intestinalinfection.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_other, file=paste0("figures/combined/fig_contours_other.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_otitismedia, file=paste0("figures/combined/fig_contours_otitismedia.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_pneumonia, file=paste0("figures/combined/fig_contours_pneumonia.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_sinusitis, file=paste0("figures/combined/fig_contours_sinusitis.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_ssti, file=paste0("figures/combined/fig_contours_ssti.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_strepphar, file=paste0("figures/combined/fig_contours_strepphar.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_tonsillitis, file=paste0("figures/combined/fig_contours_tonsillitis.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_uri, file=paste0("figures/combined/fig_contours_uri.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_uti, file=paste0("figures/combined/fig_contours_uti.pdf"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_viralinfection, file=paste0("figures/combined/fig_contours_viralinfection.pdf"), width=figwidth, height=figwidth*figaspect)

ggsave(fig_contours_asthma, file=paste0("figures/combined/fig_contours_asthma.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_bactinf, file=paste0("figures/combined/fig_contours_bactinf.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_bronchitis, file=paste0("figures/combined/fig_contours_bronchitis.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_hypertension, file=paste0("figures/combined/fig_contours_hypertension.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_influenza, file=paste0("figures/combined/fig_contours_influenza.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_intestinalinfection, file=paste0("figures/combined/fig_contours_intestinalinfection.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_other, file=paste0("figures/combined/fig_contours_other.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_otitismedia, file=paste0("figures/combined/fig_contours_otitismedia.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_pneumonia, file=paste0("figures/combined/fig_contours_pneumonia.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_sinusitis, file=paste0("figures/combined/fig_contours_sinusitis.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_ssti, file=paste0("figures/combined/fig_contours_ssti.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_strepphar, file=paste0("figures/combined/fig_contours_strepphar.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_tonsillitis, file=paste0("figures/combined/fig_contours_tonsillitis.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_uri, file=paste0("figures/combined/fig_contours_uri.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_uti, file=paste0("figures/combined/fig_contours_uti.png"), width=figwidth, height=figwidth*figaspect)
ggsave(fig_contours_viralinfection, file=paste0("figures/combined/fig_contours_viralinfection.png"), width=figwidth, height=figwidth*figaspect)


# Formatted figures for the paper: 

(fig_contours_sinusitis + labs(title=element_blank(), tag="A)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_sinusitis_formatted.png"), width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_sinusitis + labs(title=element_blank(), tag="A)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_sinusitis_formatted.pdf"), width=figwidth, height=figwidth*8/10, dpi=figres)

(fig_contours_uri + labs(title=element_blank(), tag="B)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_uri_formatted.png"), width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_uri + labs(title=element_blank(), tag="B)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_uri_formatted.pdf"), width=figwidth, height=figwidth*8/10, dpi=figres)

(fig_contours_otitismedia + labs(title=element_blank(), tag="C)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_otitismedia_formatted.png"), width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_otitismedia + labs(title=element_blank(), tag="C)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_otitismedia_formatted.pdf"), width=figwidth, height=figwidth*8/10, dpi=figres)

(fig_contours_strepphar + labs(title=element_blank(), tag="D)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_strepphar_formatted.png"), width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_strepphar + labs(title=element_blank(), tag="D)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_strepphar_formatted.pdf"), width=figwidth, height=figwidth*8/10, dpi=figres)

(fig_contours_bronchitis + labs(title=element_blank(), tag="E)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_bronchitis_formatted.png"), width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_contours_bronchitis + labs(title=element_blank(), tag="E)")) %>% 
	ggsave(file=paste0("figures/combined/fig_contours_bronchitis_formatted.pdf"), width=figwidth, height=figwidth*8/10, dpi=figres)




# =============================================================================
# R2
# =============================================================================

regdf <- geotable %>% 
	filter(MSA!="00000") %>% 
	group_by(PRIMARYCOND, STATE, MSA) %>% 
	summarise(NRX=sum(NRX), NVISITS=sum(NVISITS)) %>% 
	ungroup() %>% 
	left_join(
		(geotable %>% 
			filter(MSA!="00000") %>% 
			group_by(STATE,MSA,SEX,YEAR,AGEGRP) %>% 
			slice(1) %>% 
			select(STATE,MSA,SEX,YEAR,AGEGRP,NMEMB) %>% 
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


# overall (summing over all conditions)
regdf_all <- geotable %>% 
  filter(MSA!="00000") %>% 
  group_by(STATE, MSA) %>% 
  summarise(NRX=sum(NRX), NVISITS=sum(NVISITS)) %>% 
  ungroup() %>% 
  left_join(
    (geotable %>% 
       filter(MSA!="00000") %>% 
       group_by(STATE,MSA,SEX,YEAR,AGEGRP) %>% 
       slice(1) %>% 
       select(STATE,MSA,SEX,YEAR,AGEGRP,NMEMB) %>% 
       group_by(STATE,MSA) %>% 
       summarise(NMEMB=sum(NMEMB))),
    by=c("STATE","MSA")) %>% 
  mutate(p=NRX/NVISITS, v=NVISITS/NMEMB*1000, n=NRX/NMEMB*1000) %>% 
  select(STATE, MSA, p,v,n, NMEMB)


# with unlinked claims
regdf_all_unlinked <- geotable_wunlinked %>% 
  filter(MSA!="00000") %>% 
  group_by(STATE, MSA) %>% 
  summarise(NRX=sum(NRX), NVISITS=sum(NVISITS)) %>% 
  ungroup() %>% 
  left_join(
    (geotable %>% 
       filter(MSA!="00000") %>% 
       group_by(STATE,MSA,SEX,YEAR,AGEGRP) %>% 
       slice(1) %>% 
       select(STATE,MSA,SEX,YEAR,AGEGRP,NMEMB) %>% 
       group_by(STATE,MSA) %>% 
       summarise(NMEMB=sum(NMEMB))),
    by=c("STATE","MSA")) %>% 
  mutate(p=NRX/NVISITS, v=NVISITS/NMEMB*1000, n=NRX/NMEMB*1000) %>% 
  select(STATE, MSA, p,v,n, NMEMB)




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
	filter(PRIMARYCOND!="Other") %>% 
	left_join(select(abxtable, PRIMARYCOND, RXPKP), by="PRIMARYCOND") %>% 
	ggplot(aes(x=RXPKP, y=R2)) +
		geom_line(aes(group=PRIMARYCOND), col="lightgray") + 
		geom_point(aes(col=var)) + 
		geom_line(aes(col=var), stat="smooth", method="lm") + 
		scale_color_manual(values=c("blue","red"), labels=c("Prescriptions per visit","Visits per capita")) + 
		theme_classic() + 
		scale_y_continuous(limits=c(0,1)) + 
		scale_x_log10(breaks=c(2,5,10,25,50,100)) + 
		labs(x="Prescriptions per 1000 people", y=bquote(R^2)) + 
		theme(legend.position="bottom", legend.title=element_blank(), text=element_text(size=10))

fig_regdf_flipped %>% 
	ggsave(file=paste0("figures/combined/fig_regdf_flipped.pdf"), width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_regdf_flipped + theme(legend.position="none") + labs(tag="F)")) %>% 
	ggsave(file=paste0("figures/combined/fig_regdf_flipped_nolegend.pdf"), width=figwidth, height=figwidth*8/10, dpi=figres)

fig_regdf_flipped %>% 
	ggsave(file=paste0("figures/combined/fig_regdf_flipped.png"), width=figwidth, height=figwidth*8/10, dpi=figres)
(fig_regdf_flipped + theme(legend.position="none") + labs(tag="F)")) %>% 
	ggsave(file=paste0("figures/combined/fig_regdf_flipped_nolegend.png"), width=figwidth, height=figwidth*8/10, dpi=figres)



fool <- regdf %>% 
  filter(PRIMARYCOND!="Motor vehicle traffic") %>% 
  filter(PRIMARYCOND!="Other") %>% 
  left_join(select(abxtable, PRIMARYCOND, RXPKP), by="PRIMARYCOND")

write.csv(fool, 'regdf_with_rxpkp.csv')
