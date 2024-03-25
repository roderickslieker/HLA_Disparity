# Utilts MHC-I

getFrequency <- function(allele)
{
  cat(allele, "\n")
  test <- XML::readHTMLTable(sprintf("http://www.allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=%s&hla_allele2=%s&hla_selection=&hla_pop_selection=&hla_population=&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=equal&hla_sample_year=&hla_sample_year_pattern=equal&hla_level=&hla_level_pattern=equal&hla_show=&hla_order=order_1&standard=a",allele,allele))
  
  
  testx <- test[[3]]
  if(is.null(testx)) return()
  
  testx <- testx[testx$Allele == allele,]
  

  conv <- rio::import("./Data/Conversion.xlsx")
  
  testx$Name <- reshape2::colsplit(testx$Population," ", LETTERS[1:2])[,1]
  conv$Name <- reshape2::colsplit(conv$Country," ", LETTERS[1:2])[,1]
  
  testx$Region <- conv[match(testx$Name, conv$Name),1]
  
  testx$AlleleFrequency <- as.numeric(testx$`Allele Frequency`)
 
  testx
  
}



getFrequency <- function(allele)
{
  cat(allele, "\n")
  test <- XML::readHTMLTable(sprintf("http://www.allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=%s&hla_allele2=%s&hla_selection=&hla_pop_selection=&hla_population=&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=equal&hla_sample_year=&hla_sample_year_pattern=equal&hla_level=&hla_level_pattern=equal&hla_show=&hla_order=order_1&standard=a",allele,allele))
  
  
  testx <- test[[3]]
  if(is.null(testx)) return()
  
  testx <- testx[testx$Allele == allele,]
  

  conv <- rio::import("./Data/Conversion.xlsx")
  
  testx$Name <- reshape2::colsplit(testx$Population," ", LETTERS[1:2])[,1]
  conv$Name <- reshape2::colsplit(conv$Country," ", LETTERS[1:2])[,1]
  
  testx$Region <- conv[match(testx$Name, conv$Name),1]
  
  testx$AlleleFrequency <- as.numeric(testx$`Allele Frequency`)
 
  testx
  
}

getFrequencySpecific <- function(row, data.in)
{
  #cat(row,"\n")
  dx <- data.in[row,]
  
  test <- XML::readHTMLTable(sprintf("http://www.allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=%s&hla_allele2=%s&hla_selection=&hla_pop_selection=&hla_population=&hla_country=&hla_dataset=&hla_region=%s&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=equal&hla_sample_year=&hla_sample_year_pattern=equal&hla_level=&hla_level_pattern=equal&hla_show=&hla_order=order_1&standard=a", dx$Allele, dx$Allele, dx$Continent))
  
  
  
  testx <- test[[3]]
  if(is.null(testx)) return()
  
  #testx <- testx[testx$`Sample Size` >= 1000,]
  
  testx$Region <- dx$Continent
  
  testx$AlleleFrequency <- as.numeric(testx$`Allele Frequency`)
  #nms <- by(testx$AlleleFrequency, testx$Region, mean) %>% sort() %>% as.list() %>% names()
  #testx$Region <- factor(testx$Region, levels=nms)
  testx
  
}








getPopAlleles <- function(population)
{
  cat(population, "\n")
  test <- XML::readHTMLTable(sprintf("http://www.allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=A&hla_allele2=C&hla_selection=&hla_pop_selection=&hla_population=%s&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=equal&hla_sample_year=&hla_sample_year_pattern=equal&hla_level=2&hla_level_pattern=equal&hla_show=&hla_order=order_1&standard=a",population))
  
  
  testx <- test[[3]]
  if(is.null(testx)) return()
  
  
  testx
}


  
  

# A allele
getCoverage <- function(locus, locusName,d.sub, data.in)
{
  ha <- d.sub[grep(locus,d.sub$HLA),1]
  ha <- paste0("HLA-",ha)

  
  fq <- table(rowSums(Database[,ha]) >0) / nrow(Database)
  fq <- round(fq["TRUE"]*100,1)
  
  ha <- gsub("HLA-","",ha)
  fa.a <- data.complete[grep(locus, data.complete$Allele),]
  
  pop.med <- by(fa.a$AlleleFrequency, fa.a$Population, sum, na.rm=T) %>% as.list() %>% do.call(what=c)
  
  pop.scale <- ifelse(pop.med >= 1, pop.med, 1)
  
  fa.a.sub <- fa.a[fa.a$Allele %in% ha,]
  pop.med.sub <- by(fa.a.sub$AlleleFrequency, fa.a.sub$Population, sum, na.rm=T) %>% as.list() %>% do.call(what=c)
  
  is <- intersect(names(pop.scale), names(pop.med.sub))
  
  pop.scale <- pop.scale[match(is, names(pop.scale))]
  pop.med.sub <- pop.med.sub[match(is, names(pop.med.sub))]
  
  fa.a <- fa.a[match(is, fa.a$Population),]
  table(names(pop.med.sub) == fa.a$Population)
  
  
  fa.a$pop.med.sub <- pop.med.sub
  fa.a$pop.scale <- pop.scale
  
  fa.a$Coverage <- pop.med.sub/pop.scale
  
  
  dma.a <- fa.a %>%
    group_by(Region) %>%
    summarise(Median = median(Coverage, na.rm=T))
  dma.a$Locus <- locusName
  
  
  data.in$Median <- dma.a[match(data.in$continent, dma.a$Region),"Median",drop=T]*100
  data.in <- data.in[!is.na(data.in$Median),]
  
  data.txt <- data.in %>%
    group_by(continent) %>%
    summarise(longm = median(long), latm = median(lat), medm = mean(Median, na.rm = T))
  
  px <- ggplot(data.in) + 
    geom_map(data = data.in, map = data.in, aes(x=long, y=lat, map_id=region, fill = Median))+
    scale_fill_gradientn(colours = c("white","darkblue"), limits=c(0,100))+
    ggtitle(sprintf("Coverage of %s-allele frequencies\nIncluded in %s%% of the studies", locusName, fq))+
    theme_void()+
    theme(legend.position="bottom")+
    geom_text(data = data.txt, aes(x=longm, y=latm, label = paste0(format(round(medm, 1), nsmall = 1),"%")))
  
  list(dma.a, px)
}

getCoveragePopBased <- function(locus, locusName,d.sub, data.in)
{
  ha <- d.sub[grep(locus,d.sub$HLA),1]
  ha <- paste0("HLA-",ha)
  
  fq <- table(rowSums(Database[,ha]) >0) / nrow(Database)
  fq <- round(fq["TRUE"]*100,1)
  
  ha <- gsub("HLA-","",ha)
  fa.a <- data.complete[grep(locus, data.complete$Allele),]
  
  pop.med <- by(fa.a$AlleleFrequency, fa.a$Population, sum, na.rm=T) %>% as.list() %>% do.call(what=c)
  
  pop.scale <- ifelse(pop.med >= 1, pop.med, 1)
  
  fa.a.sub <- fa.a[fa.a$Allele %in% ha,]
  pop.med.sub <- by(fa.a.sub$AlleleFrequency, fa.a.sub$Population, sum, na.rm=T) %>% as.list() %>% do.call(what=c)
  
  is <- intersect(names(pop.scale), names(pop.med.sub))
  
  pop.scale <- pop.scale[match(is, names(pop.scale))]
  pop.med.sub <- pop.med.sub[match(is, names(pop.med.sub))]
  
  fa.a <- fa.a[match(is, fa.a$Population),]
  table(names(pop.med.sub) == fa.a$Population)
  
  fa.a$pop.med.sub <- pop.med.sub
  fa.a$pop.scale <- pop.scale
  
  fa.a$Coverage <- pop.med.sub/pop.scale
  
  
  connvx <- rio::import("./Data/Conversion_studies.xlsx")
  
  fa.a$Country <- connvx[match(fa.a$Population, connvx$Population),"Country"]
  
  dma.a <- fa.a %>%
    group_by(Country) %>%
    summarise(Median = median(Coverage, na.rm=T))
  
  fa.a$Locus <- locusName
  
  
  data.in$Median <- dma.a[match(data.in$region, dma.a$Country),"Median",drop=T]*100
  
  table(data.in[is.na(data.in$Median),"region"])
  
  
  data.txt <- data.in %>%
    group_by(region) %>%
    summarise(longm = median(long), latm = median(lat), medm = mean(Median, na.rm = T))
  data.txt <- data.txt[-grep("NaN",data.txt$medm),]
  
  px <- ggplot(data.in) + 
    geom_map(data = data.in, map = data.in, aes(x=long, y=lat, map_id=region, fill = Median))+
    scale_fill_gradientn(colours = c("white","darkblue"), limits=c(0,100))+
    ggtitle(sprintf("Coverage of %s-allele frequencies\nIncluded in %s%% of the studies", locusName, fq))+
    theme_void()+
    theme(legend.position="bottom")
  #geom_text(data = data.txt, aes(x=longm, y=latm, label = paste0(format(round(medm, 1), nsmall = 1),"%")))
  
  list(dma.a, px)
}





# A allele
getCoverageSR <- function(locus, locusName,d.sub, data.in)
{
  ha <- d.sub[grep(locus,d.sub$HLA),1] %>% unique() %>% as.character()
  
  fq <- table(rowSums(DB[,ha]) >0) / nrow(DB)
  fq <- round(fq["TRUE"]*100,1)
  
  fa.a <- data.complete[grep(locus, data.complete$Allele),]
  
  pop.med <- by(fa.a$AlleleFrequency, fa.a$Population, sum, na.rm=T) %>% as.list() %>% do.call(what=c)
  
  pop.scale <- ifelse(pop.med >= 1, pop.med, 1)
  
  fa.a.sub <- fa.a[fa.a$Allele %in% ha,]
  
  
  pop.med.sub <- by(fa.a.sub$AlleleFrequency, fa.a.sub$Population, sum, na.rm=T) %>% as.list() %>% do.call(what=c)
  
  is <- intersect(names(pop.scale), names(pop.med.sub))
  
  pop.scale <- pop.scale[match(is, names(pop.scale))]
  pop.med.sub <- pop.med.sub[match(is, names(pop.med.sub))]
  
  fa.a <- fa.a[match(is, fa.a$Population),]
  table(names(pop.med.sub) == fa.a$Population)
  
  
  fa.a$pop.med.sub <- pop.med.sub
  fa.a$pop.scale <- pop.scale
  
  fa.a$Coverage <- pop.med.sub/pop.scale
  
  
  dma.a <- fa.a %>%
    group_by(Region) %>%
    summarise(Median = median(Coverage, na.rm=T))
  dma.a$Locus <- locusName
  
  
  data.in$Median <- dma.a[match(data.in$continent, dma.a$Region),"Median",drop=T]*100
  data.in <- data.in[!is.na(data.in$Median),]
  
  data.txt <- data.in %>%
    group_by(continent) %>%
    summarise(longm = median(long), latm = median(lat), medm = mean(Median, na.rm = T))
  
  px <- ggplot(data.in) + 
    geom_map(data = data.in, map = data.in, aes(x=long, y=lat, map_id=region, fill = Median))+
    scale_fill_gradientn(colours = c("white","darkblue"), limits=c(0,100))+
    ggtitle(sprintf("Coverage of %s-allele frequencies\nIncluded in %s%% of the studies", locusName, fq))+
    theme_void()+
    theme(legend.position="bottom")+
    geom_text(data = data.txt, aes(x=longm, y=latm, label = paste0(format(round(medm, 1), nsmall = 1),"%")))
  
  list(dma.a, px)
}

getCoveragePopBasedSR <- function(locus, locusName,d.sub, data.in)
{
  ha <- d.sub[grep(locus,d.sub$HLA),1]
  ha <- paste0("HLA-",ha)
  
  fq <- table(rowSums(DB[,ha]) >0) / nrow(DB)
  fq <- round(fq["TRUE"]*100,1)
  
  fa.a <- data.complete[grep(locus, data.complete$Allele),]
  
  ha <- gsub("HLA-","",ha)
  fa.a <- data.complete[grep(locus, data.complete$Allele),]
  
  pop.med <- by(fa.a$AlleleFrequency, fa.a$Population, sum, na.rm=T) %>% as.list() %>% do.call(what=c)
  
  pop.scale <- ifelse(pop.med >= 1, pop.med, 1)
  
  fa.a.sub <- fa.a[fa.a$Allele %in% ha,]
  pop.med.sub <- by(fa.a.sub$AlleleFrequency, fa.a.sub$Population, sum, na.rm=T) %>% as.list() %>% do.call(what=c)
  
  is <- intersect(names(pop.scale), names(pop.med.sub))
  
  pop.scale <- pop.scale[match(is, names(pop.scale))]
  pop.med.sub <- pop.med.sub[match(is, names(pop.med.sub))]
  
  fa.a <- fa.a[match(is, fa.a$Population),]
  table(names(pop.med.sub) == fa.a$Population)
  
  fa.a$pop.med.sub <- pop.med.sub
  fa.a$pop.scale <- pop.scale
  
  fa.a$Coverage <- pop.med.sub/pop.scale
  
  
  connvx <- rio::import("./Conversion_studies.xlsx")
  
  fa.a$Country <- connvx[match(fa.a$Population, connvx$Population),"Country"]
  
  dma.a <- fa.a %>%
    group_by(Country) %>%
    summarise(Median = median(Coverage, na.rm=T))
  
  fa.a$Locus <- locusName
  
  
  data.in$Median <- dma.a[match(data.in$region, dma.a$Country),"Median",drop=T]*100
  
  table(data.in[is.na(data.in$Median),"region"])
  
  
  data.txt <- data.in %>%
    group_by(region) %>%
    summarise(longm = median(long), latm = median(lat), medm = mean(Median, na.rm = T))
  data.txt <- data.txt[-grep("NaN",data.txt$medm),]
  
  px <- ggplot(data.in) + 
    geom_map(data = data.in, map = data.in, aes(x=long, y=lat, map_id=region, fill = Median))+
    scale_fill_gradientn(colours = c("white","darkblue"), limits=c(0,100))+
    ggtitle(sprintf("Coverage of %s-allele frequencies\nIncluded in %s%% of the studies", locusName, fq))+
    theme_void()+
    theme(legend.position="bottom")
  #geom_text(data = data.txt, aes(x=longm, y=latm, label = paste0(format(round(medm, 1), nsmall = 1),"%")))
  
  list(dma.a, px)
}