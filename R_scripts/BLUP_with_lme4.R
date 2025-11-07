# Load packages
library(lme4)

#Load data
data = read.table("prt_training_lines.txt",header = T, sep = ",")

#Remove missing values
data = data %>% filter(!is.na(PRT))

#make a new variable location_year
data$ENV = paste0(data$Location,"_",data$Year)

plot(
  density(data$PRT),
  main = "Phenotype value distribution",
  xlab = "Phenotype values",
  ylab = "Density",
  col = "red",
  lwd = 2
)

#Run model

model = lmer(PRT ~ Loc_Year + (1 | Line), data = data)

fixed_effects <- fixef(model)

print(fixed_effects)

line_effects = ranef(model)$Line # Here line effect is the BLUPs

print(line_effects)

length(line_effects[[1]])


plot(density(line_effects[[1]]),
     main="BLUP value distribution",
     col = "red",
     lwd= 2
)



#save file 
line_effects$lines <- rownames(line_effects)

line_effects <- line_effects[,c(2,1)]

write.table(line_effects, "blup_vals_PRT.txt", quote = FALSE, row.names = FALSE)
