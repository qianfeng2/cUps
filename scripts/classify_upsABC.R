#######################################################################
# Copyright (C) 2023  Qian Feng

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Usuage instruction
# This is used in R.
# This code is to calculate the posterior probs of the query sequence.
# The only input argument is your preferred output directory 
# Usage example: cd your_project_dir \\
#                Rscript scripts/classify_upsABC.R results
#######################################################################

args = commandArgs(trailingOnly=TRUE)
output_dir = args[1]
library(data.table)


#setwd("/Users/fengqian/Downloads/upsABC_global/algorithm/algorithm2/algorithm2_finalmodel/Evaluation_MunHua/MH_update_without_mix")#set working dir

Mun_Hua_data <- fread("reference_data/MH_846_reference.csv", header=TRUE, data.table = FALSE)
visua_data=cbind(as.matrix(table(Mun_Hua_data$`DBLa domain`)),as.matrix(table(Mun_Hua_data$`DBLa domain`,Mun_Hua_data$UpsType_MCL)))
colnames(visua_data)=c("count","upsA","upsB","upsC")

prob_subclass <- visua_data[,1]/(colSums(visua_data)[1])
upsA_prob <- visua_data[,2]/visua_data[,1]
upsB_prob <- visua_data[,3]/visua_data[,1]
upsC_prob <- visua_data[,4]/visua_data[,1]
MH_data_analyze <- as.data.frame(visua_data)
MH_data_analyze$prob_subclass=prob_subclass
MH_data_analyze$upsA_prob=upsA_prob
MH_data_analyze$upsB_prob=upsB_prob
MH_data_analyze$upsC_prob=upsC_prob
#MH_data_analyze


P_xAd=fread(paste0(output_dir,"/PAD.csv"), skip = 1,data.table = TRUE)
P_xAd_matrix <- as.matrix(P_xAd[,2:ncol(P_xAd)]);rownames(P_xAd_matrix) <- P_xAd$V1

P_xBd=fread(paste0(output_dir,"/PBD.csv"), skip = 1, data.table = TRUE)
P_xBd_matrix <- as.matrix(P_xBd[,2:ncol(P_xBd)]);rownames(P_xBd_matrix) <- P_xBd$V1

P_xCd=fread(paste0(output_dir,"/PCD.csv"), skip = 1, data.table = TRUE)
P_xCd_matrix <- as.matrix(P_xCd[,2:ncol(P_xCd)]);rownames(P_xCd_matrix) <- P_xCd$V1

seqnum=dim(P_xAd_matrix)[1]# seqnum is total number of your query sequences
soft_matrix <- matrix(0,seqnum,3);
rownames(soft_matrix) <- rownames(P_xAd_matrix);colnames(soft_matrix) <- c("A","B","C")
hard_matrix <- matrix(0,2,3);
rownames(hard_matrix) <- c("Count","Normalized Prob");colnames(hard_matrix) <- c("A","B","C")
sum_soft_matrix <- matrix(0,2,3)
rownames(sum_soft_matrix) <- c("Sum Prob","Normalized Prob");colnames(sum_soft_matrix) <- c("A","B","C")
temp <- rep(0,seqnum)


for (i in 1:seqnum){
  P_xAd=rep(0,33);P_xBd=rep(0,33);P_xCd=rep(0,33);
  P_d=MH_data_analyze$prob_subclass;
  P_Ad=MH_data_analyze$upsA_prob;P_Bd=MH_data_analyze$upsB_prob;P_Cd=MH_data_analyze$upsC_prob
  
  
  P_xAd=exp(as.numeric(P_xAd_matrix[i,]));P_xBd=exp(as.numeric(P_xBd_matrix[i,]));P_xCd=exp(as.numeric(P_xCd_matrix[i,]));
  prob_A <- sum(P_xAd*P_d*P_Ad)
  prob_B <- sum(P_xBd*P_d*P_Bd)
  prob_C <- sum(P_xCd*P_d*P_Cd)
  soft_matrix[i,]=c(prob_A/(prob_A+prob_B+prob_C),prob_B/(prob_A+prob_B+prob_C),prob_C/(prob_A+prob_B+prob_C))
  max_index=which(soft_matrix[i,]==max(soft_matrix[i,]))
  if (length(max_index)>1){max_index <- sample(max_index,1)}
  temp[i]=max_index
}

hard_matrix[1,]=c(length(which(temp==1)),length(which(temp==2)),length(which(temp==3)))
hard_matrix[2,]=hard_matrix[1,]/sum(hard_matrix[1,])
sum_soft_matrix[1,]=colSums(soft_matrix)
sum_soft_matrix[2,]=sum_soft_matrix[1,]/sum(sum_soft_matrix[1,])
write.csv(soft_matrix,paste0(output_dir,"/classification_result.csv"))
