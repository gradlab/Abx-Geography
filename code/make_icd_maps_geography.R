# =============================================================================
# Import
# =============================================================================

library(tidyverse) 
library(data.table) 
library(icd)
source('code/utils.R')
source('code/utils_private.R')

# =============================================================================
# Define AHRQ ICD mappings: 
# =============================================================================

cond_priority <- tibble(
	COND=c("Pneumonia",
		"UTI",
		"Intestinal infection",
		"Strep pharyngitis",
		"Sinusitis",
		"URI (other)",
		"SSTI",
		"Otitis media",
		"Asthma",
		"Bacterial infection (other)",
		"Tonsillitis",
		"Bronchitis (acute)",
		"Influenza",
		"Viral infection",
		"Essential hypertension",
		"Motor vehicle traffic")) %>% 
	mutate(PRIORITY=1:n())

# Pneumonia
# Urinary tract infections
# Skin and subcutaneous tissue infections
# Sinusitis *****
# Other upper respiratory infections
# Asthma
# Acute bronchitis
# Influenza
# Intestinal infection
# Otitis media
# Bacterial infections
# Viral infection *****
# Essential hypertension
# Motor vehicle traffic

# ICD9:
ccs9 <- read_csv("data/temp_ccs_dxref2015.csv",
	col_types=list(
		col_character(), 
		col_double(), 
		col_character(), 
		col_skip(), 
		col_skip(), 
		col_skip()))
names(ccs9) <- c("ICD9","CCS","CCS_DESC")
ccs9 <- ccs9 %>% 
	filter(CCS %in% c(
		122,	# Pneumonia
		159,	# Urinary tract infections
		197,	# Skin and subcutaneous tissue infections
		126,	# Other upper respiratory infections
		128,	# Asthma
		125,	# Acute bronchitis
		123,	# Influenza
		124,	# Tonsillitis
		135,	# Intestinal infection
		92,		# Otitis media
		3,		# Bacterial infection, unspecified site
		7,		# Viral infection
		98,		# Essential hypertension
		2607	# Motor vehicle traffic
		)) %>% 
	mutate(COND=case_when(
		CCS==122~"Pneumonia",
		CCS==159~"UTI",
		CCS==197~"SSTI",
		CCS==126~"URI (other)",
		CCS==128~"Asthma",
		CCS==125~"Bronchitis (acute)",
		CCS==123~"Influenza",
		CCS==124~"Tonsillitis",
		CCS==135~"Intestinal infection",
		CCS==92~"Otitis media",
		CCS==3~"Bacterial infection (other)",
		CCS==7~"Viral infection",
		CCS==98~"Essential hypertension",
		CCS==2607~"Motor vehicle traffic")) %>%
	mutate(COND=case_when(
		ICD9=="0340"~"Strep pharyngitis",
		ICD9%in%c("4610",
			"4611",
			"4612",
			"4613",
			"4618",
			"4619",
			"4730",
			"4731",
			"4732",
			"4733",
			"4738",
			"4739")~"Sinusitis",
		TRUE~COND
		)) %>% 
	mutate(DX=ICD9, ICD="9") %>%
	select(COND,DX,ICD)

# ICD10: 
ccs0 <- read_csv("data/ccsr_dxref2021.csv",
	col_types=list(
		col_character(),
		col_skip(),
		col_character(),
		col_character(),
		col_skip(),
		col_skip(),
		col_skip(),
		col_skip(),
		col_skip()))
names(ccs0) <- c("ICD10","CCS","CCS_DESC") 
ccs0 <- ccs0 %>% 
	filter(CCS %in% c(
		"RSP002",	# Pneumonia
		"GEN004",	# Urinary tract infections
		"SKN001",	# Skin and subcutaneous tissue infections
		"RSP006",	# Other upper respiratory infections
		"RSP009",	# Asthma
		"RSP005",	# Acute bronchitis
		"RSP003",	# Influenza
		"RSP004",	# Tonsillitis
		"DIG001",	# Intestinal infection
		"EAR001",	# Otitis media
		"INF003",	# Bacterial infections
		"INF008",	# Viral infection
		"RSP001",	# Sinusitis
		"CIR007",	# Essential hypertension
		"EXT007"	# Motor vehicle traffic
		)) %>% 
	mutate(COND=case_when(
		CCS=="RSP002"~"Pneumonia",
		CCS=="GEN004"~"UTI",
		CCS=="SKN001"~"SSTI",
		CCS=="RSP006"~"URI (other)",
		CCS=="RSP009"~"Asthma",
		CCS=="RSP005"~"Bronchitis (acute)",
		CCS=="RSP003"~"Influenza",
		CCS=="RSP004"~"Tonsillitis",
		CCS=="DIG001"~"Intestinal infection",
		CCS=="EAR001"~"Otitis media",
		CCS=="INF003"~"Bacterial infection (other)",
		CCS=="INF008"~"Viral infection",
		CCS=="RSP001"~"Sinusitis",
		CCS=="CIR007"~"Essential hypertension",
		CCS=="EXT007"~"Motor vehicle traffic")) %>%
	mutate(COND=case_when(
		ICD10=="J020"~"Strep pharyngitis",
		TRUE~COND)) %>% 
	mutate(DX=ICD10, ICD="0") %>%
	select(COND,DX,ICD) %>% 
	distinct() 

ccs_map <- bind_rows(ccs9,ccs0) %>% 
	left_join(cond_priority, by="COND")

# Ensure ccs_map has only 1:1 diagnosis to condition, based on above priority:
ccs_map <- ccs_map %>% 
	group_by(DX,ICD) %>% 
	mutate(MAX_PRIORITY=min(PRIORITY)) %>% 
	filter(PRIORITY==MAX_PRIORITY) %>% 
	select(COND,DX,ICD,PRIORITY)

write.csv(ccs_map, file="data/geography/ccs_map_geography.csv", row.names=FALSE)
