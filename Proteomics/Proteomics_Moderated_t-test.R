library(MKmisc)

setwd('dir')
df <- read.csv('Processing_Foldchange.csv')
df = subset(df, select = -c(X) )


# M1 vs M2a
dfM1vsM2a = subset(df, select = c('M1_R1', 'M1_R2', 'M1_R3', 'M1_R4', 'M2a_R1', 'M2a_R2', 'M2a_R3', 'M2a_R4'))
X = data.matrix(dfM1vsM2a)
g <- factor(c(rep("M1", 4), rep("M2a", 4)))
testM1vsM2a = mod.t.test(X, group = g, adjust.method = "BH")
pM1vsM2a = testM1vsM2a$p.value
FDRM1vsM2a = testM1vsM2a$adj.p.value


# M1 vs M2c
dfM1vsM2c = subset(df, select = c('M1_R1', 'M1_R2', 'M1_R3', 'M1_R4', 'M2c_R1', 'M2c_R2', 'M2c_R3', 'M2c_R4'))
X = data.matrix(dfM1vsM2c)
g <- factor(c(rep("M1", 4), rep("M2c", 4)))
testM1vsM2c = mod.t.test(X, group = g, adjust.method = "BH")
pM1vsM2c = testM1vsM2c$p.value
FDRM1vsM2c = testM1vsM2c$adj.p.value

# M2a vs M2c
dfM2avsM2c = subset(df, select = c('M2a_R1', 'M2a_R2', 'M2a_R3', 'M2a_R4', 'M2c_R1', 'M2c_R2', 'M2c_R3', 'M2c_R4'))
X = data.matrix(dfM2avsM2c)
g <- factor(c(rep("M2a", 4), rep("M2c", 4)))
testM2avsM2c = mod.t.test(X, group = g, adjust.method = "BH")
pM2avsM2c = testM2avsM2c$p.value
FDRM2avsM2c = testM2avsM2c$adj.p.value

# Create df with values
values <- data.frame(pM1vsM2a,FDRM1vsM2a,pM1vsM2c,FDRM1vsM2c,pM2avsM2c, FDRM2avsM2c)
# Save df
write.csv(values, "Processing_Moderated_t-test.csv", row.names=TRUE)
