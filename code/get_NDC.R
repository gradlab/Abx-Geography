library(tidyverse) 

redbook <- read_csv("data/redbook18.csv", col_types=list(col_character(),col_character(),col_character(),col_character()))

toextract <- redbook %>% 
	filter(
		(THRCLDS %in% c(
		"Antibiot, Aminoglycosides", 
		"Antibiot, B-Lactam Antibiotics", 
		"Antibiot, Cephalosporin & Rel.", 
		"Antibiot, Chloramphenicol&Comb", 
		"Antibiot, Erythromycn&Macrolid", 
		"Antibiot, Penicillins", 
		"Antibiot, Tetracyclines", 
		"Antibiotics, Misc")) | 
		(THRDTDS=="Metronidazole & Comb." & ROADS=="Oral") | 
		(THRDTDS=="Metronidazole & Comb." & is.na(ROADS)) | 
		(THRDTDS=="Fosfomycin" & ROADS=="Oral") | 
		(THRDTDS=="Nalidixic Acid" & ROADS=="Oral") | 
		(THRDTDS=="Nalidixic Acid" & is.na(ROADS)) | 
		(THRDTDS=="Nitrofurantoin" & ROADS=="Oral") | 
		(THRDTDS=="Nitrofurantoin" & is.na(ROADS)) 
		)

write_csv(toextract, file="data/ndc_to_extract_geography.csv", quote="all")

