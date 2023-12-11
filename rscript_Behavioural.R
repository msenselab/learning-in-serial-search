##Learning regular cross-trial shifts of the target location in serial search involves awareness – an eye-tracking study #####
# by Hao Yu, Ph.D student &  Fredrik Allenmark, Ph.D
# Ludwig Maximilian University of Munich
# This code is the intellectual property of Hao Yu and Fredrik Allenmark
# please do not redistribute or reuse this code for commercial purposes.  
# If you use this code for scientific purposes, please be sure to cite our original paper.

library(R.matlab)
library(tidyverse)
library(ez)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)

datafiles <- c('01_JW_dynamic target_2022-11-22_17h07.03.843.csv',
               '02_MI_dynamic target_2022-11-24_11h11.49.780.csv',
               '03_SY_dynamic target_2022-11-24_14h20.32.429.csv',
               '04_UA_dynamic target_2022-11-24_17h28.31.222.csv',
               '05_YY_dynamic target_2022-11-24_18h41.38.992.csv',
               '06_JM_dynamic target_2022-11-25_13h46.28.630.csv',
               '07_DS_dynamic target_2022-11-25_17h19.18.478.csv',
               '08_YY_dynamic target_2022-11-29_09h22.23.107.csv',
               '09_JXL_dynamic target_2022-11-29_14h15.30.527.csv',
               '10_FF_dynamic target_2022-11-29_17h44.26.418.csv',
               '11_JH_dynamic target_2022-11-30_11h37.31.576.csv',
               '12_ZXF_dynamic target_2022-12-01_14h50.17.264.csv',
               '13_AT_dynamic target_2022-12-01_16h55.14.827.csv',
               '14_CS_dynamic target_2022-12-06_13h43.56.132.csv',
               '15_ZTY_dynamic target_2022-12-07_12h16.05.019.csv',
               '16_Lizz_dynamic target_2022-12-08_13h10.24.388.csv',
               '17_YYF_dynamic target_2022-12-09_09h48.52.719.csv',
               '18_GV_dynamic target_2022-12-12_09h43.59.205.csv',
               '19_DA_dynamic target_2022-12-12_12h18.49.973.csv',
               '20_NK_dynamic target_2022-12-12_15h17.40.434.csv',
               '21_HXY_dynamic target_2022-12-12_16h43.20.623.csv',
               '22_SU_dynamic target_2022-12-13_16h27.04.500.csv',
               '23_SN_dynamic target_2022-12-14_15h04.43.222.csv',
               '24_BJR_dynamic target_2022-12-17_10h10.39.568.csv')
 




data <- do.call('rbind',lapply(datafiles, readData))

qdata <- filter(data, is.na(trials.thisTrialN))

seqno <- seq(1,24)
order <- rep(1:4,24/4)
dyndirections <- as.numeric((seqno+1)%%8 < 4)*2 - 1
slider_responses <- seq(0,100,10)
qdata <- mutate(qdata, q1 = ifelse(qdata$key_resp1.keys > 3, 
                                   "Pattern", "Random"), 
                correct_response = ifelse(dyndirections[group] > 0, 1, 2),    
                correct_q2 = ifelse(key_resp_2.keys == correct_response, 1, 0),
                correct_q4 = ifelse(key_resp1_3.keys== 1,1,0)) %>% group_by(participant) %>% summarize(q1_resp = ifelse("Pattern" %in% q1, "Pattern", "Random"), 
                                                                                                       q2_corr = ifelse(1 %in% correct_q2, 1, 0), q4_corr = ifelse(1 %in% correct_q4, 1, 0), 
                                                                                                       awareness=ifelse("Pattern" %in% q1 & 1 %in% correct_q2 & 1 , "Aware", "Unaware"))

slider_data <- filter(data, !is.na(slider.response)) %>%  mutate(q3_resp = slider_responses[slider.response]) %>% select(participant, q3_resp)
rating_data <- filter(data, !is.na(key_resp1.keys)) %>% select(participant, q1_rating = key_resp1.keys)
qdata <- full_join(qdata, slider_data)
qdata <- full_join(qdata, rating_data)

ndata <- filter(data_new, !is.na(trials.thisTrialN)) %>% 
    mutate (block_order = ifelse(order[group] > 2, "mixlast","misxfirst"))

data %>% group_by(participant) %>% 
  mutate(search_resp.rt = ifelse(search_resp.rt > mean(search_resp.rt) + 3*sd(search_resp.rt) |  search_resp.rt< mean(search_resp.rt) - 3*sd(search_resp.rt), NA, search_resp.rt)) %>%
  ungroup() %>% filter(!is.na( search_resp.rt))-> data_new 


sdata <- group_by( data_new,block_type,tar_pos_type,awareness,participant) %>% filter(correct==1) %>%    
                                                              summarize(rt=mean(search_resp.rt)*1000)


ssdata <- sdata %>% summarize(mrt=mean(rt), sert=sd(rt)/sqrt(n()))
ssdata$block_type<-factor(ssdata$block_type,levels = c("fix","mix"),labels = c("fixed","mixed"))

sdata_awareness <- group_by(ndata, tar_pos_type, awareness, participant) %>% 
  filter(search_resp.rt<5.5, search_resp.rt>0.2, tar_pos_type!="Random") %>%
  summarize(rt=mean(search_resp.rt)*1000)

sdata_awareness %>% ezANOVA(dv=rt, wid=participant, within=tar_pos_type, between=awareness)

### bar plot###
Rtime<-ggplot(ssdata, aes(x=tar_pos_type, y=mrt,fill= tar_pos_type)) + 
  geom_bar(stat="identity") + 
  scale_fill_grey(start = 0.4,end = 0.8)+
  geom_errorbar(aes(ymin=mrt-sert, ymax=mrt+sert), width=0.2) + theme_classic() + 
  labs(x="Target position", y="Mean response time (ms)") + 
  coord_cartesian(ylim=c(1500,3000))+ 
  theme(legend.position = "none")+
  facet_grid(block_type~awareness)


#####line plots###

Rtime<-ggplot(ssdata, aes(x=tar_pos_type,y=mrt,fill= tar_pos_type,shape=block_type,color=block_type,group=block_type)) + 
  geom_line()+
  geom_point(size=3)+
  scale_fill_grey(start = 0.4,end = 0.8)+
  geom_errorbar(aes(ymin=mrt-sert, ymax=mrt+sert), width=0.2) + 
  theme_light() + 
  labs(x="Target position", y="Mean response time (ms)") + 
  coord_cartesian(ylim=c(1700,3100))+
  theme(legend.position = "none")+
  scale_color_manual(name="",values = c('black', 'dimgrey'))+
  facet_grid(~awareness)


###### violin plot with dot plot####
sidata <- group_by(ndata, tar_pos_type, block_type,awareness,participant) %>% filter(correct==1) %>% summarize(rt=mean(search_resp.rt)*1000) %>%
  spread(tar_pos_type,rt) %>% mutate( pc = `Infreq. neighbour` - `Freq. neighbour`)  #

sidata1 <- sidata %>% group_by(awareness,block_type) %>% #,
  ggplot(aes(x=awareness,y=pc))+
  geom_violin()+
  geom_point(data= sidata,aes(y=pc))+
  theme(legend.position = "none")+
  theme_classic()+
  facet_grid(~block_type)+
  labs(x="Awareness", y="Mean probability cueing effect (ms)")

#########  box plot#######

sidata <- group_by(ndata, tar_pos_type, block_type,awareness,participant) %>% filter(correct==1) %>% summarize(rt=mean(search_resp.rt)*1000) %>%
  spread(tar_pos_type,rt) %>% mutate( pc = `Infreq. neighbour` - `Freq. neighbour`) 

sidata$block_type<-factor(sidata$block_type,levels = c("fix","mix"),labels = c("Fixed","Mixed"))

sidata1_box <- sidata %>% group_by(awareness,block_type) %>% #,
  ggplot(aes(x=awareness,y=pc))+
  geom_boxplot()+
  geom_point(data= sidata,aes(y=pc))+
  theme(legend.position = "none")+
  theme_bw()+ #theme_classic()+
  facet_grid(~block_type)+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) +
 labs(x="Awareness group", y="Mean probability-cueing effect (ms)")


######## correlation plot#########

pc_Q1<-awareness_data %>% ggplot(aes(x=q1_rating, y=pc, color= awareness)) + 
  geom_point(size=2)+ 
  stat_poly_line(aes(fill= awareness))+
  stat_poly_eq() +
  labs(x="Q1 confidence rating", y="Probability-cueing effect (ms)")+
  theme_bw() +
  theme(legend.position = "top")+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) +
  labs(x="Awareness group", y="Mean probability-cueing effect (ms)")

pc_Q3<-awareness_data %>% ggplot(aes(x=q3_resp, y=pc,color= awareness)) + 
  geom_point(size=2)+
  stat_poly_line(aes(fill= awareness))+
  stat_poly_eq()+              
  labs(x="Q3 frequency rating", y="Probability-cueing effect (ms)")+ 
  theme_bw() + 
  theme(legend.position = "top")+ 
  #theme(legend.position = "none")+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) +
  theme(legend.title = element_blank())


# Mean error rates 
edata <- group_by(ndata,participant) %>%
  summarize(err=1-mean(correct))

sedata <- edata  %>% summarize(merr=mean(err), se_err=sd(err)/sqrt(n()))

####Mean outlier#####
data %>% group_by(participant)%>%summarize(mean(outlier))
mean(ndata$outlier)
######Mean incorrect####
ndata %>% group_by(participant)%>%summarize(mean(correct))

###########  inter trial priming1 ###########

## remember to change SSdata####
tar_pos_inttrial_plot3_mean <- tar_pos_inttrial3 %>%summarize(mRT=mean(rt), se_rt = sd(rt)/sqrt(n()-1)) %>%           ####此处用mean不用median
  ggplot(aes(t_t_dist, y=mRT, group=block_type)) + 
  geom_point(position=pj, size=2.8) + 
  geom_line(position=pj) + 
  theme_bw() + #calssic,bw,mimimal_grid 
  geom_errorbar(aes(ymin=mRT-se_rt, ymax=mRT+se_rt), width=0.4, position=pj) +
  theme(legend.position = "top")+
  geom_point(data=filter(ssdata, tar_pos_type!="Random"),aes(x=1,  shape=tar_pos_type,color=tar_pos_type, y=mrt), size=3.2)+
  geom_errorbar(data=filter(ssdata, tar_pos_type!="Random"), aes(x=1, y=mrt,color=tar_pos_type,ymin=mrt-sert, ymax=mrt+sert),width=0.4)+
  labs(x="Target position inter-trial change", y="Mean response time (ms)")+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) 

### order####

tar_pos_inttrial_plot4<- tar_pos_inttrial3 %>%  summarize(mRT=mean(rt), se_rt = sd(rt)/sqrt(n()-1)) %>% 
  ggplot(aes(t_t_dist, y=mRT, group=block_type, color = block_type)) + 
  geom_point(position=pj, size=4) + geom_line(position=pj) + 
  theme_classic() + 
  geom_errorbar(aes(ymin=mRT-se_rt, ymax=mRT+se_rt), width=0.4, position=pj) +
  labs(x="Target position inter-trial change", y="Mean response time (ms)")+
  theme(legend.position = "top")+
  theme(legend.title=element_blank())+
  facet_grid(~block_order)


# ---- target learning effect across blocks ----

TLE_block_data <- rt_block_data %>% spread(tar_pos_type,RT) %>%
  mutate(TLE = `Infreq. neighbour` - `Freq. neighbour`) 

# Compute average distractor interference 
sint_block_data <- summarize(TLE_block_data, mint=mean(TLE), 
                             se_int=sd(TLE)/sqrt(n()-1))

# Plot figure of average target interference across blocks
block_learning<-sint_block_data %>% 
  ggplot(aes(x=as.integer(block), y=mint,color = block_order)) +        #scale_color_grey(start = 0,end = 0.4,labels = c("Random", "Regular"))+
  geom_point(size=3) +  
  geom_line() + 
  geom_errorbar(aes(ymin=mint-se_int, ymax=mint+se_int), width=0.2) + 
  theme_classic() + 
  theme(legend.position = "top")+
  theme(legend.title = element_blank())+
  labs(x = 'Block', y = "Target learning effect (ms)") + 
  scale_x_continuous("Block", 1:10)+
  facet_wrap(awareness~.)
block_learning


########################intertrial priming2 ##############

sidata <- group_by(data_new, tar_pos_type,awareness,participant) %>% filter(correct==1) %>% summarize(rt=mean(search_resp.rt)*1000) %>% #
  spread(tar_pos_type,rt) %>% mutate( pc = `Infreq. neighbour` - `Freq. neighbour`)  

sdata3 <- group_by(data_new,participant) %>% mutate(previous_tar_pos_type = lag(tar_pos_type)) %>% group_by(block_type,tar_pos_type,previous_tar_pos_type, awareness, participant) %>% 
  filter(trials.thisTrialN %% 60 != 0,correct==1,previous_tar_pos_type!= 'Random') %>%  
  summarise(rt= mean(search_resp.rt)*1000) %>%  spread(tar_pos_type,rt) %>% mutate( pc = `Infreq. neighbour` - `Freq. neighbour`)  

Rtime3_mean<-ggplot(ssdata3, aes(x= previous_tar_pos_type, y=mpc,fill= previous_tar_pos_type)) + 
  geom_bar( stat="identity") + 
  scale_fill_grey(start = 0.4,end = 0.8)+ 
  geom_errorbar(aes(ymin=mpc-sert, ymax=mpc+sert), width=0.2)+
  theme_bw() + 
  scale_y_continuous(breaks = seq(-400, 700, 50), limits = c(-400,700))+
  theme(legend.position = "none")+
  facet_nested(~awareness+block_type)+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) +
  scale_y_continuous()+
  labs(x="Previous target position", y="Mean probability-cueing effect (ms)") 

###learning effect across blocks #####

sdata_block <- mutate(data_new, block=floor(trials.thisTrialN/120)+1,block_order = ifelse(order[group] > 2, "mixlast","misxfirst")) %>% 
  mutate(session_nr = ifelse(block<5, 1,2)) %>% 
group_by(session_nr,block_order, tar_pos_type,awareness,participant) %>% filter(correct==1,awareness =="Aware") %>% 
  summarize(rt=mean(search_resp.rt)*1000)

sdata_block  <- sdata_block  %>% spread(tar_pos_type,rt) %>%
  mutate(TLE = `Infreq. neighbour` - `Freq. neighbour`) 

###Compute average distractor interference ###
sint_block_data <- summarize(sdata_block, mint=mean(TLE), 
                             se_int=sd(TLE)/sqrt(n()-1))
sint_block_data$session_nr<- factor(sint_block_data$session_nr)
# Plot figure of average target interference across blocks
block_learning<-sint_block_data %>% 
  ggplot(aes(x=session_nr, y= mint, group= block_order,color = block_order)) +   
  geom_point(size=3) +  
  geom_line() + 
  geom_errorbar(aes(ymin=mint-se_int, ymax=mint+se_int), width=0.2) + 
  theme_classic() + 
  theme(legend.position = "top")+
  theme(legend.title = element_blank())+
  labs(x = 'session_nr', y = "target learning effect (ms)") 

####### 8 blocks
sdata_rbo <- mutate(data_new,block_order = ifelse(order[group] > 2, "mixlast","misxfirst")) %>% group_by(block_type,block_order, tar_pos_type,awareness,participant) %>% filter(correct==1,awareness =="Aware") %>%      
  summarize(rt=mean(search_resp.rt)*1000)
sdata_rbo  <- sdata_rbo  %>% spread(tar_pos_type,rt) %>%
  mutate(TLE = `Infreq. neighbour` - `Freq. neighbour`) 

########4 blocks
sdata_block4 <- mutate(data_new, block=floor(trials.thisTrialN/120)+1,block_order = ifelse(order[group] > 2, "mixlast","misxfirst")) %>% group_by(block_type,block_order, tar_pos_type,awareness,participant) %>% filter(correct==1,awareness =="Aware", block %in% c(1,2,3,4)) %>%      
  summarize(rt=mean(search_resp.rt)*1000)
sdata_block4  <- sdata_block4 %>% spread(tar_pos_type,rt) %>% mutate(TLE = `Infreq. neighbour` - `Freq. neighbour`) 
