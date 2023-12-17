# Question 3.1.----------------------------------------------------------------------------------------

#Download the sequences of Chlamydia trachomatis and Chlamydia psittaci, and read them into
#R. You can use the library seqinr (recommended commands: read.fasta and unlist; or for
#python users, SeqIO in the Biopython library). Treat each sequence as one sample and make
#some plots to explore the distributions of the four nucleotide contents for the two strains, as well
#s the GC content. (10%)
# -------------------------------------------------------------------------------------------------------



library(phylotools)
library(stringr)

# read the fasta file as a df
PSI <- read.fasta("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/FOUNDATIONS_IN_MATH_AND_STATISTICS/LAB/chlamydia_psittaci.ffn")
TRA <- read.fasta("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/FOUNDATIONS_IN_MATH_AND_STATISTICS/LAB/chlamydia_trachomatis.ffn")

#write.csv(PSI, "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/FOUNDATIONS_IN_MATH_AND_STATISTICS/LAB/PSI.csv", row.names=FALSE)
#write.csv(TRA, "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/FOUNDATIONS_IN_MATH_AND_STATISTICS/LAB/TRA.csv", row.names=FALSE)

# number of sequences
NumSeq_PSI <- nrow(PSI)
NumSeq_TRA <- nrow(TRA)

# vecotr of sequences in PSI
PSI_bases_freq <- data.frame(matrix(nrow = NumSeq_PSI, ncol = 4))
colnames(PSI_bases_freq) <- c('A_freq', 'C_freq', 'G_freq', 'T_freq')

for (i in 1:NumSeq_PSI){
  sequence_current <- PSI[i, 2]
  sequence_current <- strsplit(sequence_current, "")[[1]]
  sequence_current <- sort(sequence_current)                 # this ensures the order A --> C --> G --> T
  bases_count <- unique(as.vector(rep(table(sequence_current), table(sequence_current))))
  bases_freq <- c()
  
  for(j in bases_count){
    bases_freq <- c(bases_freq, j/sum(bases_count))
  }
  
  PSI_bases_freq$A_freq[i] <- bases_freq[1]
  PSI_bases_freq$C_freq[i] <- bases_freq[2]
  PSI_bases_freq$G_freq[i] <- bases_freq[3]
  PSI_bases_freq$T_freq[i] <- bases_freq[4]
}

# vecotr of sequences in TRA
TRA_bases_freq <- data.frame(matrix(nrow = NumSeq_TRA, ncol = 4))
colnames(TRA_bases_freq) <- c('A_freq', 'C_freq', 'G_freq', 'T_freq')

for (i in 1:NumSeq_TRA){
  sequence_current <- TRA[i, 2]
  sequence_current <- strsplit(sequence_current, "")[[1]]
  sequence_current <- sort(sequence_current)                 # this ensures the order A --> C --> G --> T
  bases_count <- unique(as.vector(rep(table(sequence_current), table(sequence_current))))
  bases_freq <- c()
  
  for(j in bases_count){
    bases_freq <- c(bases_freq, j/sum(bases_count))
  }
  
  TRA_bases_freq$A_freq[i] <- bases_freq[1]
  TRA_bases_freq$C_freq[i] <- bases_freq[2]
  TRA_bases_freq$G_freq[i] <- bases_freq[3]
  TRA_bases_freq$T_freq[i] <- bases_freq[4]
}

# histograms of base feequencies in psitacci
par(mfrow = c(2,2))
ylimit <- c(0, 0.6)
boxplot(PSI_bases_freq$A_freq, main = 'A frequency psittaci', ylim = ylimit, outline = F)
boxplot(PSI_bases_freq$C_freq, main = 'C frequency psittaci', ylim = ylimit, outline = F)
boxplot(PSI_bases_freq$G_freq, main = 'G frequency psittaci', ylim = ylimit, outline = F)
boxplot(PSI_bases_freq$T_freq, main = 'T frequency psittaci', ylim = ylimit, outline = F)

# histograms of base feequencies in trachomatis
par(mfrow = c(2,2))
ylimit <- c(0, 0.6)
boxplot(TRA_bases_freq$A_freq, main = 'A frequency trachomatis', ylim = ylimit, outline = F)
boxplot(TRA_bases_freq$C_freq, main = 'C frequency trachomatis', ylim = ylimit, outline = F)
boxplot(TRA_bases_freq$G_freq, main = 'G frequency trachomatis', ylim = ylimit, outline = F)
boxplot(TRA_bases_freq$T_freq, main = 'T frequency trachomatis', ylim = ylimit, outline = F)



PSI_sequence_concat <- c()
for (i in 1:NumSeq_PSI){
  PSI_sequence_concat <- paste(PSI_sequence_concat, PSI[i,2], sep = "")
}
PSI_sequence_concat <- strsplit(PSI_sequence_concat, "")[[1]]

TRA_sequence_concat <- c()
for (i in 1:NumSeq_TRA){
  TRA_sequence_concat <- paste(TRA_sequence_concat, TRA[i,2], sep = "")
}
TRA_sequence_concat <- strsplit(TRA_sequence_concat, "")[[1]]



GC_content <- function(sequence_in){
  GC_content <- 0
  
  for(bp in sequence_in){
    if (bp == "A" || bp == "T" || bp == "C" || bp == "G"){
      if(bp == "G"){GC_content <- GC_content + 1}
      else if (bp == "C"){GC_content <- GC_content + 1}
      else{GC_content <- GC_content}
    }
  }
  seq_in_len <- length(sequence_in)
  return(GC_content/seq_in_len)
}

PSI_GC_content <- GC_content(PSI_sequence_concat)
TRA_GC_content <- GC_content(TRA_sequence_concat)

PSI_GC_content
TRA_GC_content



# -----------------------------------------------------------------------------------------------------
#Compute a 99% empirical bootstrap interval for the GC content of Chlamydia trachomatis based
#on resampling (recommended command: sample) 200 times from the original data. Repeat the
#same for Chlamydia psittaci. (10%)
# ----------------------------------------------------------------------------------------------------


PSI_sequences <- data.frame(matrix(nrow = NumSeq_PSI, ncol = 1))
  
for(i in 1:NumSeq_PSI){
  sequence_current <- PSI[i, 2]
  PSI_sequences[i,1] <- sequence_current
}


TRA_sequences <- data.frame(matrix(nrow = NumSeq_TRA, ncol = 1))

for(i in 1:NumSeq_TRA){
  sequence_current <- TRA[i, 2]
  TRA_sequences[i,1] <- sequence_current
}



# we want to have 200 boot samples each of 1119 sequences
#PSI_boot <- replicate(nrow(PSI_sequences),sample(PSI_sequences,200,replace = TRUE))

# first compute the GC content for all sequences and then make bootstrap samples from these GC contents
# for every bootstrap sample calculate the mean GC content
# then sort the means and get the bootstrap estimator

PSI_GC_contents <- PSI_sequences

for (i in 1:nrow(PSI_GC_contents)){
  x <- strsplit(PSI_GC_contents[i,1], "")[[1]]
  PSI_GC_contents[i,1] <- GC_content(x)
}

colnames(PSI_GC_contents) <- c("GC_content")
head(PSI_GC_contents)


TRA_GC_contents <- TRA_sequences

for (i in 1:nrow(TRA_GC_contents)){
  x <- strsplit(TRA_GC_contents[i,1], "")[[1]]
  TRA_GC_contents[i,1] <- GC_content(x)
}

colnames(TRA_GC_contents) <- c("GC_content")
head(TRA_GC_contents)


PSI_GC_contents <- as.numeric(as.vector(PSI_GC_contents$GC_content))
TRA_GC_contents <- as.numeric(as.vector(TRA_GC_contents$GC_content))

# now we need to make 200 bootstrap samples of these GC contents
# then for every bootstrap sample we need to compute the mean GC content which will be the GC content for the bacterium



PSI_GC_boot <- replicate(length(PSI_GC_contents), sample(PSI_GC_contents, 200, replace = TRUE))
# each row of PSI_GC_boot is a bootstrap sample (so we have 200 rows)
# for every row compute the mean GC content
PSI_GC_boot_estimators <- 1:200

for(i in 1:200){
  PSI_GC_boot_estimators[i] <- mean(PSI_GC_boot[i,])
}

# sort the bootstrap estimators
PSI_GC_boot_estimators <- sort(PSI_GC_boot_estimators)

# to get a 99% bootstrap confidence interval of the GC content get the 2nd and 199th element of the PSI_GC_boot_estimatros sorted array
PSI_boot_lower <- PSI_GC_boot_estimators[2]
PSI_boot_upper <- PSI_GC_boot_estimators[199]
print(c("We are 99% sure that the true value of the GC content of psitacchi is in the range", PSI_boot_lower ,"-", PSI_boot_upper))




# now to the same but for the other bacterium
TRA_GC_boot <- replicate(length(TRA_GC_contents), sample(TRA_GC_contents, 200, replace = TRUE))
# each row of PSI_GC_boot is a bootstrap sample (so we have 200 rows)
# for every row compute the mean GC content
TRA_GC_boot_estimators <- 1:200

for(i in 1:200){
  TRA_GC_boot_estimators[i] <- mean(TRA_GC_boot[i,])
}

# sort the bootstrap estimators
TRA_GC_boot_estimators <- sort(TRA_GC_boot_estimators)

# to get a 99% bootstrap confidence interval of the GC content get the 2nd and 199th element of the TRA_GC_boot_estimatros sorted array
TRA_boot_lower <- TRA_GC_boot_estimators[2]
TRA_boot_upper <- TRA_GC_boot_estimators[199]
print(c("We are 99% sure that the true value of the GC content of trachomatis is in the range", TRA_boot_lower ,"-", TRA_boot_upper))





# --------------------------------------------------------------------------------------------------------------------
#Based on the bootstrapping results, can GC content be used to distinguish the two strains? (5%)
# ------------------------------------------------------------------------------------------------------------------

# Yes, because the ranges do not overlap and the confidence is high.




# -------------------------------------------------------------------------------------------------------------------
#Next, we want to implement a simple desicion tree to classify a sequence into one of the two
#strains. Implement the gini criterion with weighted average and compute the average gini impurity at a GC content of 40%. (10%)
# --------------------------------------------------------------------------------------------------------------------

# 1 --> Psitacchi
# 0 --> Trachomatis
tree_data <- data.frame(matrix(nrow = nrow(PSI_sequences) + nrow(TRA_sequences), ncol = 2))
colnames(tree_data) <- c("GC.content", "species")

for (i in 1:nrow(PSI_sequences)){
  tree_data$GC.content[i] <- GC_content(strsplit(PSI_sequences[i,1], "")[[1]])
  tree_data$species[i] <- 1
}


j <- 1
for (i in nrow(PSI_sequences):nrow(tree_data)){
  if (j <= nrow(TRA_sequences)){
    tree_data$GC.content[i] <- GC_content(strsplit(TRA_sequences[j,1], "")[[1]])
    tree_data$species[i] <- 0
  }
  j <- j + 1
}

tree_data <- na.omit(tree_data)

Gini <- function(dataset, split_variable_idx, outcome_idx, split_point){
  
    L <- dataset[dataset[, split_variable_idx] < split_point, ]
    R <- dataset[dataset[, split_variable_idx] > split_point, ]
    
    # number of 0s in L divided by the length of L and squared
    L_gini_1 <- (nrow(L[L[, outcome_idx] == 0, ])/nrow(L))^2
    # number of 1s in L divided by the length of L and squared
    L_gini_2 <- (nrow(L[L[, outcome_idx] == 1, ])/nrow(L))^2
    # sum of the above
    L_gini <- 1 - (L_gini_1 + L_gini_2)
    
    # number of os in R divided by the length of R and squared
    R_gini_1 <- (nrow(R[R[, outcome_idx] == 0, ])/nrow(R))^2
    # number of 1s in R divided by the length of R and sqaured
    R_gini_2 <- (nrow(R[R[, outcome_idx] == 1, ])/nrow(R))^2
    # sum of the above
    R_gini <- 1 - (R_gini_1 + R_gini_2)
    
    # add the weights
    L_gini <- L_gini*(nrow(L)/nrow(dataset))
    R_gini <- R_gini*(nrow(R)/nrow(dataset))
    # sum of the above
    gini_at_split_point <- L_gini + R_gini
    
    return(gini_at_split_point)
  
}

avg_Gini_at_40 <- Gini(dataset = tree_data, split_variable_idx = 1, outcome_idx = 2, split_point = 0.4)

print(c("The average Gini impurity at a GC content of 40% is", avg_Gini_at_40))




# ------------------------------------------------------------------------------------------------------------------
#Plot the the average gini impurity as a function of GC content. What is the optimal GC content
#for the first split of a decision tree?
# ------------------------------------------------------------------------------------------------------------------

GCs <- seq(1:100)
GCs2 <- 1:length(GCs)

for (i in GCs){
  GCs2[i] <- i/100.00
}

ginis <- c()
for(i in GCs2){
  res <- Gini(dataset = tree_data, split_variable_idx = 1, outcome_idx = 2, split_point = i)
  ginis <- c(ginis, res)
}

df <- data.frame(GCs2, ginis)
df <- na.omit(df)
par(mfrow = c(1,1))
plot(df$GCs2, df$ginis, main = "Impurity by GC content split point", col = "blue", xlab = "GC content split point", ylab = "Gini impurity", type = 'p')
