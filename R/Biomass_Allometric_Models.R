########################################################################################################################################
# Title: Above-ground biomass calculations using Updated generalized biomass equations for North American tree species 
#        using species specific generalized geometrical equations   
# Writer: Soren Donisvitch 
# Date 4/09/22
# Cited material: Chojnacky, D.C., Heath, L.S., Jenkins, J.C., 2014. Updated generalized biomass equations for North American tree species. 
#                 Forestry 87 (1), 129,151. https://doi.org/10.1093/forestry/cpt053.USDA Forest Service, Environmental Monitoring and Assessment 
#                 Program Forest Health Monitoring Manual (Eastern and Western), Appendix F
#
# Contact: soren.donisvitch@gmail.com
# Requirements: function requires tree species (FIA codes) and diameter measurements (DBH in cm) 
# Disclaimer: The author is not liable for any use and or application of these scripts and does not claim accuracy or concurrency of its models.
########################################################################################################################################

# Above ground biomass function
# parameters and formula from Table 5 from Chojnacky,2014.
# Where (with biomass in kg and diameter in cm) the biomass equation is: ln(biomass) = b0 + b1 ln(diameter) after mathematical transformation: biomass = e^b0 * diameter^b1
# 
# FIA_species_code = FIA species code 
# dbh = diameter at breast height in cm
# will return above ground live tree biomass in kg

AG_Biomass = function(FIA_species_code,dbh){
  dbh = as.numeric(dbh) # force numeric 
  dbh = ifelse(dbh > 0,dbh,NA) # fills missing and some non-real instances 
  
  if(FIA_species_code %in% c(12,16,19)){ #Abies , 0.35 spg*
    b0 = -2.3123 
    b1 = 2.3482
  } 
  else if(FIA_species_code %in% c(10,11,15,17,20,22)){ #Abies >= 0.35 spg
    b0 = -3.1774	
    b1 = 2.6426
  }
  else if(FIA_species_code %in% c(241)){ #Cupressaceae , 0.30 spg
    b0 = -1.9615	
    b1 = 2.1063
  }
  else if(FIA_species_code %in% c(81,212,242)){ #Cupressaceae 0.30 - 0.39 spg
    b0 = -2.7765	
    b1 = 2.4195
  }
  else if(FIA_species_code %in% c(42,68)){ #Cupressaceae  0.40 spg
    b0 = -2.6327	
    b1 = 2.4757
  }
  else if(FIA_species_code %in% c(71,73,70)){ #Larix
    b0 = -2.3012	
    b1 = 2.3853
  }
  else if(FIA_species_code %in% c(93,98)){ #Picea , 0.35 spg
    b0 = -3.03	
    b1 = 2.5567
  }
  else if(FIA_species_code %in% c(91,94,95,97)){ #Picea >= 0.35 spg
    b0 = -2.1364	
    b1 = 2.3233
  }
  else if(FIA_species_code %in% c(100,101,135,105,108,116,117,118,119,122,125,129)){ #Pinus , 0.45 spg
    b0 = -2.6177	
    b1 = 2.4638
  }
  else if(FIA_species_code %in% c(110,111,121,126,131)){ #Pinus >= 0.45 spg
    b0 = -3.0506	
    b1 = 2.6465
  }
  else if(FIA_species_code %in% c(202)){ #Pseudotsuga
    b0 = -2.4623	
    b1 = 2.4852
  }
  else if(FIA_species_code %in% c(261)){ #Tsuga , 0.40 spg
    b0 = -2.348	
    b1 = 2.3876
  }
  else if(FIA_species_code %in% c(263,264)){ #Tsuga >= 0.40 spg
    b0 = -2.9208	
    b1 = 2.5697
  }
  else if(FIA_species_code %in% c(312,315,316,317,319)){ #Aceraceae , 0.50 spg*
    b0 = -2.047	
    b1 = 2.3852
  }
  else if(FIA_species_code %in% c(318)){ #Aceraceae >= 0.50 spg
    b0 = -1.8011	
    b1 = 2.3852
  }
  else if(FIA_species_code %in% c(351,350,376)){ #Betulaceae , 0.40 spg
    b0 = -2.5932	
    b1 = 2.5349
  }
  else if(FIA_species_code %in% c(375,379,378)){ #Betulaceae 0.40 - 0.49 spg
    b0 = -2.2271	
    b1 = 2.4513
  }
  else if(FIA_species_code %in% c(371,370)){ #Betulaceae 0.50 - 0.59 spg
    b0 = -1.8096	
    b1 = 2.348
  }
  else if(FIA_species_code %in% c(372,701,703)){ #Betulaceae >= 0.60 spg
    b0 = -2.2652	
    b1 = 2.5349
  }
  else if(FIA_species_code %in% c(491)){ #Cornaceae/Ericaceae/
    b0 = -2.2118	
    b1 = 2.4133
  }
  else if(FIA_species_code %in% c(691)){ #Lauraceae/Platanaceae/
    b0 = -2.2118	
    b1 = 2.4133
  }
  else if(FIA_species_code %in% c(693,361,711,981,931,731,356,761,762,763,935,972,970)){ #Rosaceae/Ulmaceae1
    b0 = -2.2118	
    b1 = 2.4133
  }
  else if(FIA_species_code %in% c(404)){ #Fabaceae/Juglandaceae
    b0 = -2.5095	
    b1 = 2.6175
  }
  else if(FIA_species_code %in% c(407,403,402,409,601)){ #Carya
    b0 = -2.5095	
    b1 = 2.6175
  }
  else if(FIA_species_code %in% c(901)){ #Fabaceae/Juglandaceae
    b0 = -2.5095	
    b1 = 2.5437
  }
  else if(FIA_species_code %in% c(421,531,802,806,809,812,823,827,832,833,800,835,837)){ #Fagaceae, deciduous
    b0 = -2.0705	
    b1 = 2.441
  }
  else if(FIA_species_code %in% c(431,631,807,820,841)){ #Fagaceae, evergreen
    b0 = -2.2198	
    b1 = 2.441
  }
  else if(FIA_species_code %in% c(611)){ #Hamamelidaceae
    b0 = -2.639	
    b1 = 2.5466
  }
  else if(FIA_species_code %in% c(332)){ #Hippocastanaceae
    b0 = -2.4108
    b1 = 2.4177
  }
  else if(FIA_species_code %in% c(951,952)){ #Tiliaceae
    b0 = -2.4108	
    b1 = 2.4177
  }
  else if(FIA_species_code %in% c(621,655,653)){ #Magnoliaceae
    b0 = -2.5497
    b1 = 2.5011
  }
  else if(FIA_species_code %in% c(543,544,540)){ #Oleaceae , 0.55 spg
    b0 = -2.0314	
    b1 = 2.3524
  }
  else if(FIA_species_code %in% c(541)){ #Oleaceae >= 0.55 spg
    b0 = -1.8384	
    b1 = 2.3524
  }
  else if(FIA_species_code %in% c(741,747,470)){ #Salicaceae , 0.35 spg
    b0 = -2.6863	
    b1 = 2.4561
  }
  else if(FIA_species_code %in% c(742,743,740,746,927,920)){ #Salicaceae >= 0.35 spg
    b0 = -2.4441	
    b1 = 2.4561
  }
  else if(FIA_species_code %in% c(50,69,64,65)){ #Cupressaceae
    b0 = -2.7096
    b1 = 2.1942
  }
  else if(FIA_species_code %in% c(755,475,477)){ #Fabaceae/Rosaceae
    b0 = -2.9255	
    b1 = 2.4109
  }
  else if(FIA_species_code %in% c(807,814,843,850,804)){ #Fagaceae
    b0 = -3.0304	
    b1 = 2.4982
  }
  else if(FIA_species_code %in% c(140,106,133)){ #Pinaceae
    b0 = -3.2007	
    b1 = 2.5339
  }
  else if(FIA_species_code %in% c(001,998)){ #composite of all hardwoods
    b0 = -22.2702
    b1 = 2.45436
  }
  else if(FIA_species_code %in% c(002)){ #composite of all conifers/softwoods
    b0 = -22.5944
    b1 = 2.4470
  }
  else if(FIA_species_code %in% c(003,999,"NA")){ #composite of all species
    b0 = -2.23982
    b1 = 2.4529
  } else {
    b0 = -2.23982
    b1 = 2.4529
  }
  biomass = ifelse(!is.na(dbh),(exp(1)^b0)*(dbh^b1),NA) #ln(biomass) = b0 + b1 ln(diameter) or biomass = e^b0 * dbh^b1
  return(biomass)
}
