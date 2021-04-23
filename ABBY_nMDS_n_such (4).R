################################################################################
#   Code written by Joe Nolan, Dept. Biol. Sciences, West Liberty University   #
#   Tested using R ver. 4.0.4 on platform x86_64-w64-mingw32/x64 (64-bit)      #
################################################################################
####   Package installation checks 

if(!require(vegan)) install.packages("vegan", dependencies = TRUE)
if(!require(stringr)) install.packages("stringr", dependencies = TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
library(vegan); library(stringr); library(ggplot2);  library(ggrepel); 
library(plyr); library(rlang); library(dplyr); library(kableExtra);

library(dplyr)

####   Loading, cleaning, and data assembly
df1 <- read.csv("Jan_April.csv")
df2 <- read.csv("May_Aug.csv")
df3 <- read.csv("Sept_Oct.csv")
df4 <- read.csv("Nov_Dec.csv")

df1 <- df1[ ,c(1,4,8)]; df2 <- df2[ ,c(1,4,8)]
df3 <- df3[ ,c(1,4,8)]; df4 <- df4[ ,c(1,4,8)]

df5 <- rbind(df1, df2, df3, df4)
junk <- c("df1", "df2", "df3", "df4")
rm(list = junk); rm(junk)

chem <- read.csv("ABBY_chem.csv", stringsAsFactors = TRUE)

####   nMDS of Months
#   Calculation of cross-tabulation Month × Species
df6 <- subset(df5, barcode == "barcode01" | barcode == "barcode02" |barcode == "barcode03" |
               barcode == "barcode04" |barcode == "barcode05" |barcode == "barcode06" |
               barcode == "barcode07" |barcode == "barcode08" |barcode == "barcode09" |
               barcode == "barcode10" |barcode == "barcode11" |barcode == "barcode12")


df9 <- df6[complete.cases(df6), ]                       # ANYTHING WITH BARCODE "NA" IS REMOVED HERE

df9$moBar <- paste(df9$month, df9$barcode)


#mobarcounts <- as.data.frame(table(df9$species, df9$moBar))
#names(mobarcounts) <- c("species","moBar","count")
#head(mobarcounts)
#moCo10 <- subset(mobarcounts, count >10 & species != "")

#monthcounts <- as.data.frame(table(df9$species, df9$month))
#names(monthcounts) <- c("species","month","count")
#head(monthcounts)
#month10 <- subset(monthcounts, count >10 & species != "")



for(i in 1:length(moCo10$moBar)){
  moCo10$Month[i] <- strsplit(as.character(moCo10$moBar[i]), "\\s+")[[1]][1]
  moCo10$Barcode[i] <- strsplit(as.character(moCo10$moBar[i]), "\\s+")[[1]][2]
}
df <- moCo10



df9 <- ddply(df9,.(month, barcode, moBar, species),nrow)
View(df)

df <- subset(df9, V1 >= 10)

dfsortmonth <- read.csv("dfsortmonth.csv")
dfsortb <- read.csv("dfsortb.csv")
dfsorta <- read.csv("dfsorta.csv")

kable(dfsorta) %>% 
  kable_styling(bootstrap_options = c("striped", "hover"))
  
#col.names = c("Site","Coordinates", "Mean", "Median", "Min", "Max", "Km from mouth"), digits = 1)



#data.count <- ddply(df,.(moBar),nrow)
#head(data.count)
#write.csv(data.count, 'datacount.csv')
#datacount <- read.csv("datacount.csv")
#head(datacount)



spByMonth <- as.matrix(table(df$moBar, df$species))
#months <- row.names(spByMonth)
#   Generation of Bray-Curtis (default) distance matrix
dist_matrix <- vegdist(x = spByMonth); rm(spByMonth)
nMDS_scores <- monoMDS(dist_matrix, k = 2, maxit = 5000)
data.scores <- as.data.frame(scores(nMDS_scores))
data.scores$MoBar <- rownames(data.scores)
head(data.scores)



for(i in 1:length(data.scores$MoBar)){
  #data.scores$Months[i] <- str_to_title(data.scores$Months[i])
  data.scores$Month[i] <- strsplit(data.scores$MoBar[i], " ")[[1]][1]
  data.scores$Barcode[i] <- strsplit(data.scores$MoBar[i], " ")[[1]][2]
}
tail(data.scores)


#   Examination of generic stress plot
stressplot(nMDS_scores)


# adding vectors
anti_join(data.scores, chem, by = "MoBar")
anti_join(chem, data.scores, by = "MoBar")
chem <- chem[chem$MoBar != "april barcode 06", ]
chem <- chem[chem$MoBar != "july barcode 03", ]
chem <- chem[chem$MoBar != "march barcode 03", ]
chem <- chem[chem$MoBar != "march barcode 09", ]
chem <- chem[chem$MoBar != "november barcode 03", ]
chem <- chem[chem$MoBar != "november barcode 08", ]
chem <- chem[chem$MoBar != "october barcode 03", ]
chem <- chem[chem$MoBar != "october barcode 08", ]
chem <- chem[chem$MoBar != "january barcode 03", ]
chem <- chem[chem$MoBar != "february barcode 02", ]

head(chem)

fit <- envfit(nMDS_scores, chem, permutations = 9999)
arrow <- data.frame(scores(fit, "vectors"),
                    R = fit$vectors$r,
                    P = fit$vectors$pvals)
arrow$Measurement <- rownames(arrow)
rownames(arrow) <- NULL
arrow.p <- filter(arrow, P < 0.05)
arrow.r <- filter(arrow, R > 0.1)

View(arrow)

####   Plotting nMDS results
data.scores$Month <- factor(data.scores$Month, 
                            levels = c("january", "february", "march", "april",
                                       "may","june","july","august", 
                                       "september", "october", "november", 
                                       "december"))
data.scores$Barcode <- factor(data.scores$Barcode, 
                               levels = c("barcode01", "barcode02", "barcode03", 
                                          "barcode04", "barcode05", "barcode06", 
                                          "barcode07", "barcode08", "barcode09",
                                          "barcode10", "barcode11", "barcode12"),
                              labels = c("Ranch Bar", "Cricket Hollow", "Grandstaff Run", "Britt's Run", "County Line Bridge", 
                              "Elm Grove Park", "Junior Ave", "Washington Ave", "Kroger", "Fulton", "Tunnel Green", "Mouth"))


ggplot() +
  #geom_segment(data = arrow.p,
        #       aes(x=0, y=0, xend = 4*MDS1, yend = 4*MDS2),
         #      arrow = arrow(length = unit(0.3, "cm")*arrow.p$R,
          #                   type = "open"), size = 0.2) +
  geom_point(data = data.scores, aes(x = MDS1, y = MDS2, col = Barcode), 
             alpha = 0.8, size = 3) +
  #scale_shape_manual(values=c(17, 17, 17, 17, 17, 16, 16, 16, 16, 16, 16, 16)) +
 #geom_point(data = data.scores, aes(x = MDS1, y = MDS2, pch = Barcode) +
    #         alpha = 0.8, pch = 21, 
     #        colour = "black", size = 3)) +
 #geom_text_repel(data = data.scores, aes(x = MDS1, y = MDS2, label = Barcode),
   #       size = 3, hjust = -0.4, vjust = -0.6, alpha = 0.2) +
  #scale_colour_grey()+
coord_equal() +
facet_wrap(~Month) +
  theme_bw() +
  labs(x = "MDS1", y = "MDS2")+
  theme(text = element_text(family = "serif"),
        axis.text.x = element_text(size=14),  #x-axis marker text
        axis.text.y = element_text(size=14), # y-axis marker text
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10))

ggplot() +
  geom_segment(data = arrow.r,
              aes(x=0, y=0, xend = 4*MDS1, yend = 4*MDS2),
               arrow = arrow(length = unit(0.3, "cm")*arrow.r$R,
                          type = "open"), size = 1) +
  geom_point(data = data.scores, aes(x = MDS1, y = MDS2, col = Month, shape = Barcode), 
             alpha = 0.8, size = 3) +
scale_shape_manual(values=c(17, 17, 17, 17, 17, 16, 16, 16, 16, 16, 16, 16)) +
  #geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=Month), size=1, linetype=1)
 geom_text(data = arrow.r, aes(x = 4*MDS1, y = 4*MDS2, label = Measurement),
           size = 6, hjust = 1, vjust = -0, alpha = 0.5) +
  #geom_text_repel(data = data.scores, aes(x = MDS1, y = MDS2, label = Barcode),
        #size = 3, hjust = -0.4, vjust = -0.6, alpha = 0.2) +
  #scale_colour_grey()+
  #coord_equal() +
  #facet_wrap(~Month)+
  theme_bw() +
  labs(x = "MDS1", y = "MDS2")+
  theme(text = element_text(family = "serif"),
        axis.text.x = element_text(size=14),  #x-axis marker text
        axis.text.y = element_text(size=14), # y-axis marker text
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10))

#rm(data.scores); rm(nMDS_scores)


###ellipses

ord <- ordiellipse(nMDS_scores, data.scores$Month, display = "sites", 
                 kind = "se", conf = 0.95, label = T)

df_ell <- data.frame()
for(g in levels(data.scores$Month)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(data.scores[data.scores$Month==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}

#  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=Month), size=1, linetype=1) 

ord <- ordiellipse(nMDS_scores, data.scores$Barcode, display = "sites", 
                  kind = "se", conf = 0.95, label = T)

df_ell <- data.frame()
for(g in levels(data.scores$Barcode)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(data.scores[data.scores$Barcode==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}


#  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=Barcode), size=1, linetype=1) 

####   Generating relative proportions, by barcode
#pathogens <- subset(df, grepl('Acinetobacter baumanii|Aeromonas hydrophyla|Arcobacter butzleri|
                  #  Bacillus cereus|Campylobacter jejuni|Citrobacter freundii|Clostridium perfringens|
                  #  Desulphovibrio|Fusobacterium necrophorum biovar A|Helicobacter pylori|
                  #  Klebsiella pneumoniae|Legionella dumoffii|Legionella pneumophila|
                  #  Mycobacterium avium|Mycobacterium intracellulare|Plesiomonas shigelloides|
                  #  Pseudomonas aeruginosa|Salmonella enterica serovar Enteritidis|
                  #  Salmonella enterica serovar Paratyphi|Salmonella enterica serovar Typhi|
                  #  Salmonella enterica serovar Typhimurium|Salmonella newport|Serratia marcescens|
                  #  Shigella boydii|Shigella dysenteriae|Shigella flexneri|Shigella sonnei|
                  #  Staphylococcus aureus|Streptococcus pneumoniae|Vibrio cholerae|Vibrio parahaemolyticus|
                  #  Vibrio vulnificus|Yersinia enterocolitica$', species))

pathogens <- subset(df, grepl('Bacillus cereus|Campylobacter jejuni|Campylobacter sputorum|
                 |Legionella pneumophila|Mycobacterium avium|Mycobacterium intracellulare|
                  Pseudomonas aeruginosa|Salmonella enterica||Salmonella newport|
                  Shigella boydii|Shigella dysenteriae|Shigella flexneri|Shigella sonnei|
                  Staphylococcus aureus$', species))
                    

path <- df
path$species <- word(path$species, 1, 2)
View(path)

pathogens <- subset(path, species == "Bacillus cereus"|species == "Campylobacter jejuni"| species == "Camplobacter sputorum"|
                      species == "Legionella pneumophila"| 
                      species == "Mycobacterium avium"|species == "Mycobacterium intracellulare"|
                      species == "Pseudomonas aeruginosa"|
                      species == "Salmonella enterica"|
                      species == "Shigella boydii"|species == "Shigella dysenteriae"|species == "Shigella flexneri"|
                      species == "Shigella sonnei"|
                      species == "Staphylococcus aureus"| species == "Escherichia coli")

View(pathogens)
pathogens$moBar <- paste(pathogens$month, pathogens$barcode)

#barpath <- prop.table(table(pathogens$species, pathogens$barcode), 2)

                  
pathogens$month <- factor(pathogens$month, levels = c("january", "february", 
                            "march","april", "may", "june", "july","august", 
                            "september", "october", "november", "december"))

pathogens$barcode <- factor(pathogens$barcode, levels = c("barcode01", "barcode02", "barcode03", 
                                                          "barcode04", "barcode05", "barcode06", 
                                                          "barcode07", "barcode08", "barcode09",
                                                          "barcode10", "barcode11", "barcode12"),
                            labels = c("Ranch Bar", "Cricket Hollow", "Grandstaff Run", "Britt's Run", 
                                       "County Line Bridge", "Elm Grove Park", "Junior Ave", "Washington Ave", 
                                       "Kroger", "Fulton", "Tunnel Green", "Mouth"))



pathogens$moBar <- factor(pathogens$moBar, levels = c("january barcode01", "january barcode02", "january barcode03", "january barcode04", "january barcode05", "january barcode06", "january barcode07", "january barcode08", "january barcode09", "january barcode10", "january barcode11", "january barcode12", 
                                                      "february barcode01", "february barcode02", "february barcode03", "february barcode04", "february barcode05", "february barcode06",  "february barcode07", "february barcode08", "february barcode09", "february barcode10", "february barcode11", "february barcode12",
                                                      "march barcode01", "march barcode02", "march barcode03", "march barcode04", "march barcode05", "march barcode06", "march barcode07", "march barcode08", "march barcode09", "march barcode10", "march barcode11", "march barcode12", 
                                                      "april barcode01", "april barcode02", "april barcode03", "april barcode04", "april barcode05", "april barcode06", "april barcode07", "april barcode08", "april barcode09", "april barcode10", "april barcode11", "april barcode12", 
                                                      "may barcode01", "may barcode02", "may barcode03", "may barcode04", "may barcode05", "may barcode06", "may barcode07", "may barcode08", "may barcode09", "may barcode10", "may barcode11", "may barcode12", 
                                                      "june barcode01", "june barcode02", "june barcode03", "june barcode04", "june barcode05", "june barcode06", "june barcode07", "june barcode08", "june barcode09", "june barcode10", "june barcode11", "june barcode12", 
                                                      "july barcode01", "july barcode02", "july barcode03", "july barcode04", "july barcode05", "july barcode06", "july barcode07", "july barcode08", "july barcode09", "july barcode10", "july barcode11", "july barcode12",   
                                                      "august barcode01", "august barcode02", "august barcode03", "august barcode04", "august barcode05", "august barcode06", "august barcode07", "august barcode08", "august barcode09", "august barcode10", "august barcode11", "august barcode12",   
                                                      "september barcode01", "september barcode02", "september barcode03", "september barcode04", "september barcode05", "september barcode06", "september barcode07", "september barcode08", "september barcode09", "september barcode10", "september barcode11", "september barcode12", 
                                                      "october barcode01", "october barcode02", "october barcode03", "october barcode04", "october barcode05", "october barcode06", "october barcode07", "october barcode08", "october barcode09", "october barcode10", "october barcode11", "october barcode12", 
                                                      "november barcode01", "november barcode02", "november barcode03", "november barcode04", "november barcode05", "november barcode06", "november barcode07", "november barcode08", "november barcode09", "november barcode10", "november barcode11", "november barcode12", 
                                                      "december barcode01", "december barcode02", "december barcode03", "december barcode04", "december barcode05", "december barcode06", "december barcode07", "december barcode08", "december barcode09", "december barcode10", "december barcode11", "december barcode12"),
                          
                          labels = c ("January Ranch Bar", "January Cricket Hollow", "January Grandstaff Run", "January Britt's Run", "January County Line Bridge", "January Elm Grove Park", "January Junior Ave", "January Washington Ave", "January Kroger", "January Fulton", "January Tunnel Green", "January Mouth", 
                                      "February Ranch Bar", "February Cricket Hollow", "February Grandstaff Run", "February Britt's Run", "February County Line Bridge", "February Elm Grove Park",  "February Junior Ave", "February Washington Ave", "February Kroger", "February Fulton", "February Tunnel Green", "February Mouth",
                                      "March Ranch Bar", "March Cricket Hollow", "March Grandstaff Run", "March Britt's Run", "March County Line Bridge", "March Elm Grove Park", "March Junior Ave", "March Washington Ave", "March Kroger", "March Fulton", "March Tunnel Green", "March Mouth", 
                                      "April Ranch Bar", "April Cricket Hollow", "April Grandstaff Run", "April Britt's Run", "April County Line Bridge", "April Elm Grove Park", "April Junior Ave", "April Washington Ave", "April Kroger", "April Fulton", "April Tunnel Green", "April Mouth", 
                                      "May Ranch Bar", "May Cricket Hollow", "May Grandstaff Run", "May Britt's Run", "May County Line Bridge", "May Elm Grove Park", "May Junior Ave", "May Washington Ave", "May Kroger", "May Fulton", "May Tunnel Green", "May Mouth", 
                                      "June Ranch Bar", "June Cricket Hollow", "June Grandstaff Run", "June Britt's Run", "June County Line Bridge", "June Elm Grove Park", "June Junior Ave", "June Washington Ave", "June Kroger", "June Fulton", "June Tunnel Green", "June Mouth", 
                                      "July Ranch Bar", "July Cricket Hollow", "July Grandstaff Run", "July Britt's Run", "July County Line Bridge", "July Elm Grove Park", "July Junior Ave", "July Washington Ave", "July Kroger", "July Fulton", "July Tunnel Green", "July Mouth",   
                                      "August Ranch Bar", "August Cricket Hollow", "August Grandstaff Run", "August Britt's Run", "August County Line Bridge", "August Elm Grove Park", "August Junior Ave", "August Washington Ave", "August Kroger", "August Fulton", "August Tunnel Green", "August Mouth",   
                                      "September Ranch Bar", "September Cricket Hollow", "September Grandstaff Run", "September Britt's Run", "September County Line Bridge", "September Elm Grove Park", "September Junior Ave", "September Washington Ave", "September Kroger", "September Fulton", "September Tunnel Green", "September Mouth", 
                                      "October Ranch Bar", "October Cricket Hollow", "October Grandstaff Run", "October Britt's Run", "October County Line Bridge", "October Elm Grove Park", "October Junior Ave", "October Washington Ave", "October Kroger", "October Fulton", "October Tunnel Green", "October Mouth", 
                                      "november Ranch Bar", "november Cricket Hollow", "november Grandstaff Run", "november Britt's Run", "november County Line Bridge", "november Elm Grove Park", "november Junior Ave", "november Washington Ave", "november Kroger", "november Fulton", "november Tunnel Green", "november Mouth", 
                                      "december Ranch Bar", "december Cricket Hollow", "december Grandstaff Run", "december Britt's Run", "december County Line Bridge", "december Elm Grove Park", "december Junior Ave", "december Washington Ave", "december Kroger", "december Fulton", "december Tunnel Green", "december Mouth"))
                                                      
pathogens$species <- word(pathogens$species, 1)
#View(pathogens)



p <- ggplot(pathogens, aes(x= moBar, y= V1, group= species)) +
  scale_y_continuous(trans='log10') +
  geom_line(aes(color=species), size = 1)+
  geom_point(aes(color=species), size = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
p


ggplot(pathogens, aes(x = barcode, fill = factor(species)))+
  geom_bar(position = "fill") +
  facet_wrap( ~ month) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8))

ggplot(pathogens, aes(x = month, fill = factor(species)))+
  geom_bar(position = "fill") +
  facet_wrap( ~ barcode) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8))

#ggplot(df, aes(x = month, fill = factor(species)))+
 # geom_bar(position = "fill") +
  #facet_wrap( ~ barcode) +
  #theme_bw() +
  #theme(axis.text.x = element_text(angle = 90),
   #     panel.grid = element_blank(),
    #    legend.title = element_text(size = 8), 
     #   legend.text = element_text(size = 8))


write.csv(pathogens, 'pathogens1.csv')
pathogens <- read.csv("pathogens1.csv")
head(pathogens)

path <- ddply(pathogens,.(moBar,species),nrow)

path <- subset(path, V1 >= 10)

head(path)
                    
pathcount <- ddply(path,.(moBar),nrow)

View(pathcount)
write.csv(pathcount, 'pathcounter1.csv')


monthave <- read.csv("month.csv")

pc <- read.csv("pathcounter1.csv")
head(pc)

ggplot(data = pc, aes(y = count, x = E_coli_MPN), stat = "identity", fill = "grey") + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("MPN") +
  ylab("Pathogen Count") + 
  scale_x_continuous(trans='log10') +
  theme_bw() +
  theme_classic()

ggplot(data = pc, aes(y = count, x = Total_Coliforms), stat = "identity", fill = "grey") + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("MPN") +
  ylab("Pathogen Count") + 
  scale_x_continuous(trans='log10') +
  theme_bw() +
  theme_classic()

plot(pc$count ~ pc$E_coli_MPN, xlab = "E. coli MPN (CFU/100mL)", ylab = "Count of Pathogenic Species")
reg_model <- lm(count ~ E_coli_MPN, data = pc)
abline(lm(pc$count ~ pc$E_coli_MPN))


plot(pc$count ~ pc$Total_Coliforms, xlab = "Total Coliform MPN (CFU/100mL)", ylab = "Count of Pathogenic Species")
reg_model <- lm(count ~ Total_Coliforms, data = pc)
abline(lm(pc$count ~ pc$Total_Coliforms))

plot(pc$count ~ pc$Development, xlab = "Percent Development", ylab = "Count of Pathogenic Species")
reg_model <- lm(count ~ Total_Coliforms, data = pc)
abline(lm(pc$count ~ pc$Total_Coliforms))


RsquareAdj(lm(count ~ E_coli_MPN, data = pc))
RsquareAdj(lm(count ~ Total_Coliforms, data = pc))



coliforms <- subset(df, grepl('Escherichia|Klebsiella|Enterobacter|Serratia|Citrobacter|Raoultella|Shigella|Hafnia$', species))
barTab <- prop.table(table(coliforms$species, coliforms$barcode), 2)
barTab

coliforms$month <- factor(coliforms$month, levels = c("january", "february", 
                          "march","april", "may", "june", "july","august", 
                          "september", "october", "november", "december"))

coliforms$moBar <- paste(coliforms$month, coliforms$barcode)

coliforms$moBar <- factor(coliforms$moBar, levels = c("january barcode01", "january barcode02", "january barcode03", "january barcode04", "january barcode05", "january barcode06", "january barcode07", "january barcode08", "january barcode09", "january barcode10", "january barcode11", "january barcode12", 
                                                      "february barcode01", "february barcode02", "february barcode03", "february barcode04", "february barcode05", "february barcode06",  "february barcode07", "february barcode08", "february barcode09", "february barcode10", "february barcode11", "february barcode12",
                                                      "march barcode01", "march barcode02", "march barcode03", "march barcode04", "march barcode05", "march barcode06", "march barcode07", "march barcode08", "march barcode09", "march barcode10", "march barcode11", "march barcode12", 
                                                      "april barcode01", "april barcode02", "april barcode03", "april barcode04", "april barcode05", "april barcode06", "april barcode07", "april barcode08", "april barcode09", "april barcode10", "april barcode11", "april barcode12", 
                                                      "may barcode01", "may barcode02", "may barcode03", "may barcode04", "may barcode05", "may barcode06", "may barcode07", "may barcode08", "may barcode09", "may barcode10", "may barcode11", "may barcode12", 
                                                      "june barcode01", "june barcode02", "june barcode03", "june barcode04", "june barcode05", "june barcode06", "june barcode07", "june barcode08", "june barcode09", "june barcode10", "june barcode11", "june barcode12", 
                                                      "july barcode01", "july barcode02", "july barcode03", "july barcode04", "july barcode05", "july barcode06", "july barcode07", "july barcode08", "july barcode09", "july barcode10", "july barcode11", "july barcode12",   
                                                      "august barcode01", "august barcode02", "august barcode03", "august barcode04", "august barcode05", "august barcode06", "august barcode07", "august barcode08", "august barcode09", "august barcode10", "august barcode11", "august barcode12",   
                                                      "september barcode01", "september barcode02", "september barcode03", "september barcode04", "september barcode05", "september barcode06", "september barcode07", "september barcode08", "september barcode09", "september barcode10", "september barcode11", "september barcode12", 
                                                      "october barcode01", "october barcode02", "october barcode03", "october barcode04", "october barcode05", "october barcode06", "october barcode07", "october barcode08", "october barcode09", "october barcode10", "october barcode11", "october barcode12", 
                                                      "november barcode01", "november barcode02", "november barcode03", "november barcode04", "november barcode05", "november barcode06", "november barcode07", "november barcode08", "november barcode09", "november barcode10", "november barcode11", "november barcode12", 
                                                      "december barcode01", "december barcode02", "december barcode03", "december barcode04", "december barcode05", "december barcode06", "december barcode07", "december barcode08", "december barcode09", "december barcode10", "december barcode11", "december barcode12"),
                          
                          labels = c ("January Ranch Bar", "January Cricket Hollow", "January Grandstaff Run", "January Britt's Run", "January County Line Bridge", "January Elm Grove Park", "January Junior Ave", "January Washington Ave", "January Kroger", "January Fulton", "January Tunnel Green", "January Mouth", 
                                      "February Ranch Bar", "February Cricket Hollow", "February Grandstaff Run", "February Britt's Run", "February County Line Bridge", "February Elm Grove Park",  "February Junior Ave", "February Washington Ave", "February Kroger", "February Fulton", "February Tunnel Green", "February Mouth",
                                      "March Ranch Bar", "March Cricket Hollow", "March Grandstaff Run", "March Britt's Run", "March County Line Bridge", "March Elm Grove Park", "March Junior Ave", "March Washington Ave", "March Kroger", "March Fulton", "March Tunnel Green", "March Mouth", 
                                      "April Ranch Bar", "April Cricket Hollow", "April Grandstaff Run", "April Britt's Run", "April County Line Bridge", "April Elm Grove Park", "April Junior Ave", "April Washington Ave", "April Kroger", "April Fulton", "April Tunnel Green", "April Mouth", 
                                      "May Ranch Bar", "May Cricket Hollow", "May Grandstaff Run", "May Britt's Run", "May County Line Bridge", "May Elm Grove Park", "May Junior Ave", "May Washington Ave", "May Kroger", "May Fulton", "May Tunnel Green", "May Mouth", 
                                      "June Ranch Bar", "June Cricket Hollow", "June Grandstaff Run", "June Britt's Run", "June County Line Bridge", "June Elm Grove Park", "June Junior Ave", "June Washington Ave", "June Kroger", "June Fulton", "June Tunnel Green", "June Mouth", 
                                      "July Ranch Bar", "July Cricket Hollow", "July Grandstaff Run", "July Britt's Run", "July County Line Bridge", "July Elm Grove Park", "July Junior Ave", "July Washington Ave", "July Kroger", "July Fulton", "July Tunnel Green", "July Mouth",   
                                      "August Ranch Bar", "August Cricket Hollow", "August Grandstaff Run", "August Britt's Run", "August County Line Bridge", "August Elm Grove Park", "August Junior Ave", "August Washington Ave", "August Kroger", "August Fulton", "August Tunnel Green", "August Mouth",   
                                      "September Ranch Bar", "September Cricket Hollow", "September Grandstaff Run", "September Britt's Run", "September County Line Bridge", "September Elm Grove Park", "September Junior Ave", "September Washington Ave", "September Kroger", "September Fulton", "September Tunnel Green", "September Mouth", 
                                      "October Ranch Bar", "October Cricket Hollow", "October Grandstaff Run", "October Britt's Run", "October County Line Bridge", "October Elm Grove Park", "October Junior Ave", "October Washington Ave", "October Kroger", "October Fulton", "October Tunnel Green", "October Mouth", 
                                      "november Ranch Bar", "november Cricket Hollow", "november Grandstaff Run", "november Britt's Run", "november County Line Bridge", "november Elm Grove Park", "november Junior Ave", "november Washington Ave", "november Kroger", "november Fulton", "november Tunnel Green", "november Mouth", 
                                      "december Ranch Bar", "december Cricket Hollow", "december Grandstaff Run", "december Britt's Run", "december County Line Bridge", "december Elm Grove Park", "december Junior Ave", "december Washington Ave", "december Kroger", "december Fulton", "december Tunnel Green", "december Mouth"))

coliforms$barcode <- factor(coliforms$barcode, levels = c("barcode01", "barcode02", "barcode03", 
                                                         "barcode04", "barcode05", "barcode06", 
                                                         "barcode07", "barcode08", "barcode09",
                                                         "barcode10", "barcode11", "barcode12"),
                            labels = c("Ranch Bar", "Cricket Hollow", "Grandstaff Run", "Britt's Run", 
           "County Line Bridge", "Elm Grove Park", "Junior Ave", "Washington Ave", 
           "Kroger", "Fulton", "Tunnel Green", "Mouth"))


coliforms$species <- word(coliforms$species, 1)

p <- ggplot(coliforms, aes(x= moBar, y= V1, group= species)) +
  scale_y_continuous(trans='log10') +
  geom_line(aes(color=species), size = 1)+
  geom_point(aes(color=species), size = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
p


ggplot(coliforms, aes(x = barcode, fill = factor(species)))+
  geom_bar(position = "fill") +
  facet_wrap( ~ month) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8))


# Diversity Indices by Month
monthlist <- levels(as.factor(df$month))
monthlist <- factor(monthlist, levels = c("january", "february", "march",
                                          "april", "may", "june", "july", 
                                          "august", "september", "october",
                                          "november", "december"))
for(i in 1:length(monthlist)){
  print(
    paste(
  levels(monthlist)[i], "Shannon-Weiner:", 
  round(diversity(table(#df$barcode[df$month == monthlist[i]], 
                  df$species[df$month == levels(monthlist)[i]])), 2)
    )
  )
}       
# Diversity Indices by Barcode
barcodelist <- levels(as.factor(df$barcode))
barcodelist <- factor(barcodelist, levels = c("barcode01", "barcode02", "barcode03", 
                                              "barcode04", "barcode05", "barcode06", 
                                              "barcode07", "barcode08", "barcode09",
                                              "barcode10", "barcode11", "barcode12"))
for(i in 1:length(barcodelist)){
  print(
    paste(
      levels(barcodelist)[i], "Shannon-Weiner:", 
      round(diversity(table(#df$month[df$barcode == barcodelist[i]], 
        df$species[df$barcode == levels(barcodelist)[i]])), 2)
    )
  )
}       


divdev <- read.csv("diversity_development.csv")
head(divdev)

divdev$Site <- factor(divdev$Site, levels = c("Ranch Bar","Cricket Hollow", "Grandstaff Run", "Britt's Run", "County Line Bridge", "Elm Grove Park", "Junior Ave", "Washington Ave", "Kroger", "Fulton", "Tunnel Green", "Mouth"))


ggplot(divdev, aes(x = Site, y = Shannon)) +
  geom_bar(stat = "identity") +
  labs(title = "Diversity") +
  xlab("Site") +
  ylab("Shannon-Weiner") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45))



# Plotting diversity vs development
divdev <- read.csv("diversity_development.csv")

plot(divdev$Shannon.Weiner ~ divdev$Developed_LMH)
plot(divdev$Shannon.Weiner ~ divdev$CSO)
plot(divdev$Shannon.Weiner ~ divdev$CSO_total)
plot(divdev$Shannon.Weiner ~ divdev$Mean)

plot(divdev$Shannon ~ divdev$Developed_total, xlab = "Total developed land use (%)", ylab = "Shannon-Weiner diversity score")
reg_model <- lm(Shannon ~ Developed_total, data = divdev)
abline(lm(divdev$Shannon ~ divdev$Developed_total))
summary(reg_model)

RsquareAdj(lm(Shannon.Weiner ~ Developed_LMH, data = divdev))
RsquareAdj(lm(Shannon ~ Developed_total, data = divdev))
RsquareAdj(lm(Shannon ~ Developed_LMH, data = divdev))

lm()
