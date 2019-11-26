library(ggpubr)
#the res is the LM output for Indian dataset
#the bb is the LM output for UK dataset 

bb_res_match <- rownames(res)[match(rownames(bb), rownames(res))]
bb_res_match <- as.data.frame(bb_res_match)
bb_res_match <-bb_res_match[which(bb_res_match$bb_res_match != "NA"),]
bb_res_match <- as.data.frame(bb_res_match)

#extrapolate the beta, sd and p values for these probes(from 
#bb_res_match) from both datasets
res_m <- res[which(bb_res_match$bb_res_match %in% rownames(res)),c(1,2,3,4,5,6)]
bb_m <- bb[which(bb_res_match$bb_res_match %in% rownames(bb)),
           c('Age.Coeff', 'Age.SE', 'Age.P', 'Sex.Coeff.M', 'Sex.SE', 'Sex.P')]

#merge the two lm data as they should have the same probes
bb_res_comb <- merge(res_m, bb_m, by  = "row.names")
rownames(bb_res_comb) <- bb_res_comb[,1]
bb_res_comb <- bb_res_comb[,-1]

bb_res_age <- bb_res_comb[which(bb_res_comb$Age_P <= 0.05),] #filter based on significant Indian age probes
bb_res_sex <- bb_res_comb[which(bb_res_comb$Sex_P <= 0.05),] #filter based on significant Indian sex probes

#The UK dataset had a tenth decrease in their coeffs
bb_res_age2 <- bb_res_age
bb_res_age2$Age.Coeff <- bb_res_age$Age.Coeff*100
bb_res_sex2 <- bb_res_sex
bb_res_sex2$Sex.Coeff.M <- bb_res_sex2$Sex.Coeff.M*100

pdf('effect magnitute.pdf', onefile = TRUE, w = 10, h= 10) #change the plot title 
ggscatter(bb_res_age2, x = 'Age_Beta', y = 'Age.Coeff',
          add = 'reg.line',
          add.params = list(color = 'blue', fill = 'lightgray'),
          conf.int = T) +
  stat_cor(method = 'pearson') +
  labs(title = "Effect size of significant Indian sites against UK sites on age (frontal)",
       x = 'India', y = 'UK') +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(bb_res_sex2, x = 'Sex_Beta', y = 'Sex.Coeff.M',
          add = 'reg.line',
          add.params = list(color = 'blue', fill = 'lightgray'),
          conf.int = T) +
  stat_cor(method = 'pearson') +
  labs(title = "Effect size of significant Indian sites against UK sites on sex (frontal)",
       x = 'India', y = 'UK') +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
