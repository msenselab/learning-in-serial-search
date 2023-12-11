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
library(hdf5r)
library(RColorBrewer)
library(viridis)
library(lattice)
library(cowplot)
library(rstatix)
library(ggpmisc)
library(BayesFactor)


# ---- Convert data ----

files <- c('01_JW_dynamic target_2022-11-22_17h07.03.843',
           '02_MI_dynamic target_2022-11-24_11h11.49.780',
           '03_SY_dynamic target_2022-11-24_14h20.32.429',
           '04_UA_dynamic target_2022-11-24_17h28.31.222',
           '05_YY_dynamic target_2022-11-24_18h41.38.992',
           '06_JM_dynamic target_2022-11-25_13h46.28.630',
           '07_DS_dynamic target_2022-11-25_17h19.18.478',
           '08_YY_dynamic target_2022-11-29_09h22.23.107',
           '09_JXL_dynamic target_2022-11-29_14h15.30.527',
           '10_FF_dynamic target_2022-11-29_17h44.26.418',
           '11_JH_dynamic target_2022-11-30_11h37.31.576',
           '12_ZXF_dynamic target_2022-12-01_14h50.17.264',
           '13_AT_dynamic target_2022-12-01_16h55.14.827',
           '14_CS_dynamic target_2022-12-06_13h43.56.132',
           '15_ZTY_dynamic target_2022-12-07_12h16.05.019',
           '16_Lizz_dynamic target_2022-12-08_13h10.24.388',
           '17_YYF_dynamic target_2022-12-09_09h48.52.719',
           '18_GV_dynamic target_2022-12-12_09h43.59.205',
           '19_DA_dynamic target_2022-12-12_12h18.49.973',
           '20_NK_dynamic target_2022-12-12_15h17.40.434',
           '21_HXY_dynamic target_2022-12-12_16h43.20.623',
           '22_SU_dynamic target_2022-12-13_16h27.04.500',
           '23_SN_dynamic target_2022-12-14_15h04.43.222',
           '24_BJR_dynamic target_2022-12-17_10h10.39.568')

# ---- Classify fixations ----

all_fixations <- readRDS('data/fixation_summary.rds')
all_saccades <- readRDS('data/saccade_summary.rds')

loc_theta <- seq(0,315,45)*pi/180
loc_x <- 7*cos(-loc_theta)
loc_y <- 7*sin(-loc_theta)

fix_radius <- 2.5

#####seqno############
search_fixations <- search_fixations %>%  #look for which position was fixed
  mutate(which_fix = ifelse(fix_1==1,1,ifelse(fix_2==1,2,ifelse(fix_3==1,3,ifelse(fix_4==1,4,ifelse(fix_5==1,5,ifelse(fix_6==1,6,
                                                                                                                      ifelse(fix_7==1,7,ifelse(fix_8==1,8,NA))))))))) %>%
  group_by(participant, tno) %>% mutate(fix_change = dists[which_fix - lag(which_fix)+8])

search_fixations %>% filter(!is.na(fix_change)) %>% mutate(direction=dyndirections[group]) %>% group_by(direction) %>% 
  ggplot(aes(x=fix_change,y=..density..)) + geom_histogram() + facet_wrap(direction~.) + theme_classic()                 

search_saccades <- all_saccades %>% filter(!is.na(tar_pos)) 

fix_patterns <- group_by(search_fixations, participant, tno) %>%
  summarize(fix_pattern = fix_pattern(fix_type)) 

search_fixations <- full_join(search_fixations, fix_patterns, by=c('participant','tno'))


# ---- Find first target fixations and saccades ----

# ---- Count number of saccades until reaching target ----
tar_sacc_count1 <- first_tar_sacc %>% group_by(tar_pos_type, block_type,awareness, participant) %>%
  filter(correct == 1, !outlier,awareness == "Aware") %>%summarize(sc=mean(sacc_count)) #filter(SearchResponse.corr==1)

tar_sacc_count1$block_type<-factor(tar_sacc_count1$block_type,levels = c("fix","mix"),labels = c("Fixed","Mixed"))

tar_sacc_count_fig <- tar_sacc_count1 %>% summarize(msc=mean(sc), se_sc=sd(sc)/sqrt(n()-1)) %>%
  ggplot(aes(x=block_type, y=msc, color=tar_pos_type, group=tar_pos_type)) +  
  geom_point(position=pj, size=3.6) + 
  geom_errorbar(aes(ymin=msc-se_sc, ymax=msc+se_sc), width=0.2, position=pj,size=0.8) +
  geom_line(position=pj,size=0.8) + 
  theme_bw() +                                                                                        # facet_wrap(~awareness) + 
  labs(x="Target constancy", y="Number of saccades until reaching target", color="Target transition condition")+
  theme(legend.position = "top",legend.title=element_blank())+
  scale_color_manual(values = c("#000000", "#696969","#A9A9A9"),breaks=c("Random", "Infreq. neighbour","Freq. neighbour"),labels = c("Random", "Infrequent", "Frequent"))+
  coord_cartesian(xlim=c(1.2,2))+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) 


# ##########---- Plot first target fixation proportion ---- 

search_fixations6<-search_fixations %>% 
  group_by(tar_pos_type, block_type, participant, awareness,tno) %>% 
  filter(fix_type=="Target") %>% filter(correct==1, !outlier,awareness =="Aware" ) %>% # "Aware""Unaware" 
  summarize(dur=sum(duration)*1000) %>%  group_by(tar_pos_type, block_type,awareness, participant) %>%  #within a trial
  summarize(dur = mean(dur)) %>%    #across traisl
  summarize(mdur = mean(dur), se_dur=sd(dur)/sqrt(n()-1))  

fix_dur_fig <-search_fixations6  %>% 
  ggplot(aes(x=block_type, y=mdur, color=tar_pos_type, group=tar_pos_type)) + 
  geom_point( size=4) + 
  geom_line(size=0.8) + 
  theme_bw() + # facet_wrap(~awareness) + 
  geom_errorbar(aes(ymin=mdur - se_dur, ymax=mdur + se_dur), width=0.1, size=0.8) +
  scale_color_manual(values = c("#000000", "#696969","#A9A9A9"),breaks=c("Random", "Infreq. neighbour","Freq. neighbour"),labels = c("Random", "Infrequent", "Frequent"))+
  theme(legend.position = "top",legend.title=element_blank())+
  labs(x="Target constancy", y="Total target fixation duration",color="Target transition condition")+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) 

# Time until first saccade on target - all trials

tar_lat_dat <- first_tar_sacc %>% group_by(tar_pos_type,block_type,awareness, participant) %>% #the last one is participant,then it is summarise trials, when it is trials, then it is summarize participants
  filter( correct==1, !outlier,awareness == "Aware") %>% # correct==1,
  summarize(latency=median(latency)) 
pj = position_dodge(width = 0.4)

first_sacc_tar_fig_aware <- tar_lat_dat %>% 
  summarize(mlatency = mean(latency), se_lat=sd(latency)/(sqrt(n()-1))) %>%
  ggplot(aes(x=block_type, y=mlatency, color=tar_pos_type, group=tar_pos_type)) + 
  geom_point(position=pj, size=3.6) + 
  geom_errorbar(aes(ymin=mlatency-se_lat, ymax=mlatency+se_lat),width=0.2, position=pj,size=0.8) +
  geom_line(position=pj,size=0.8) + 
  theme_classic() + 
  scale_color_manual(values = c("#000000", "#696969","#A9A9A9"),breaks=c("Random", "Infreq. neighbour","Freq. neighbour"),labels = c("Random", "Infrequent", "Frequent"))+
  theme(legend.position = "top",legend.title=element_blank())+
  labs(x="Target constancy", y="Average time before first saccade to target (ms)")+ 
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) 


######### the latencies of the first saccade to the frequent, random, infrequent locations

# Find time of first fixation on target
first_fix_time <- search_fixations %>% mutate(dyn_dir=dyndirections[group], 
                                                  f_t_dist = dists[dyn_dir*(which_fix - p_tar_pos)+8],
                                                  fix_cond = ifelse(is.na(which_fix), "Non-search", 
                                                                    ifelse(f_t_dist==0, "Repeat",
                                                                           ifelse(f_t_dist==1, "Freq.",
                                                                                  ifelse(f_t_dist ==-1, "Infreq.", "Random"))))) %>%
  group_by(awareness, block_type,  fix_cond , participant) %>%   # it is the proportion of looking at the target directly
  filter(fix_type!="Fixation") %>% filter(fix_count==min(fix_count)) %>% 
  select(tno, participant, first_fix_time = start_time)

# Find first saccade to target
  search_saccades_new <- full_join(search_saccades,first_fix_time , by=c("tno", "participant","awareness","block_type"))

 first_fix_sacc <-  search_saccades_new %>% filter(start_time < first_fix_time) %>% 
  group_by(tno,participant) %>% filter(start_time==max(start_time)) %>% 
  mutate(latency=(start_time - fixcross_shape.started)*1000)


 first_fix_lat_dat <- first_fix_sacc %>% group_by( fix_cond , block_type,awareness,participant) %>%
   filter( correct == 1, !outlier,!fix_cond=="NA", !fix_cond=="Non-search", awareness=="Aware",latency>50,fix_cond!="Random") %>%  #SearchResponse.corre==1,
   summarize(latency=mean(latency)) 
 
 pj = position_dodge(width = 0.4)  #（saccdic latency)
 Fig_first_tar_lat_dat_fig <-  first_fix_lat_dat  %>% 
   summarize(mlatency = mean(latency), se_lat=sd(latency)/(sqrt(n()-1))) %>%
   ggplot(aes(x= block_type, y=mlatency, color= fix_cond , group= fix_cond )) + geom_point(position=pj, size=2) + 
   geom_errorbar(aes(ymin=mlatency-se_lat, ymax=mlatency+se_lat),width=0.8, position=pj) +
   geom_line(position=pj) +
   theme_classic() +  
   #facet_wrap(~awareness) + 
   labs(x="Block", y=" Average latencies of the first saccade to different conditions (ms)", color="Target transition condition")
 
##### average fixation duration before fixating the target ####

search_fixations_ave<-search_fixations %>% 
 group_by(tar_pos_type, block_type, participant, awareness,tno) %>% filter(correct==1, !outlier,awareness =="Aware", fix_count< first_tar_fix) %>% 
   summarize(dur=mean(duration)*1000) %>%  group_by(tar_pos_type, block_type,awareness, participant) %>%  #within a trial
   summarize(dur = mean(dur))
   
 search_fixations_ave$block_type<-factor( search_fixations_ave$block_type,levels = c("fix","mix"),labels = c("Fixed","Mixed"))
 
 average_fix_dur <-  search_fixations_ave  %>% 
   summarize(mdur = mean(dur), se_dur=sd(dur)/sqrt(n()-1)) %>% 
   ggplot(aes(x= block_type, y=mdur, color= tar_pos_type, group=tar_pos_type )) + 
   geom_point(position=pj, size=3.5) + 
   geom_errorbar(aes(ymin=mdur-se_dur, ymax=mdur+se_dur),width=0.2, position=pj,size=0.8) +
   geom_line(position=pj,size = 0.8) +
   theme_bw() + 
   theme(legend.position = "top",legend.title=element_blank())+
   scale_color_manual(values = c("#000000", "#696969","#A9A9A9"),breaks=c("Random", "Infreq. neighbour","Freq. neighbour"),labels = c("Random", "Infrequent", "Frequent"))+
   labs(x="Target constancy", y="Average pre-target fixation duration (ms)", color="Target transition condition")+
   theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
       legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
   theme(strip.text = element_text(face="bold", size=9)) 


# ---- Plot first fixation proportion (line plot) ---- 

dyndirections <- as.numeric((seqno+1)%%8 < 4)*2 - 1
first_fix_prop <-search_fixations %>% mutate(dyn_dir=dyndirections[group], 
                                             f_t_dist = dists[dyn_dir*(which_fix - p_tar_pos)+8],
                                             fix_cond = ifelse(is.na(which_fix), "Non-search", 
                                                               ifelse(f_t_dist==0, "Repeat",
                                                                      ifelse(f_t_dist==1, "Freq.",
                                                                             ifelse(f_t_dist ==-1, "Infreq.", "Random"))))) %>%
  group_by(awareness, block_type, tar_pos_type, participant,tno) %>%  
  filter(fix_type!="Fixation") %>% filter(fix_count==min(fix_count)) %>% 
  filter(correct==1, tno %% 60 > 0, fix_cond!="Non-search") %>%
  summarize(repeated_prop = mean(fix_cond=="Repeat"), freq_prop = mean(fix_cond=="Freq."), infreq_prop=mean(fix_cond=="Infreq."),
            rand_prop=mean(fix_cond=="Random")) %>% 
  gather("type", "ffp", -c("awareness", "block_type", "tar_pos_type", "participant")) # !outlier 


first_fix_prop$type <- factor(first_fix_prop$type, levels=c("infreq_prop", "repeated_prop", "freq_prop", "rand_prop"), 
                              labels=c("Inreq.", "Repeated", "Freq.", "Random"), ordered = TRUE)

sfirst_fix_prop <- first_fix_prop %>% group_by(type, awareness,  block_type, tar_pos_type) %>%
  summarize(mffp = mean(ffp), se_ffp = sd(ffp)/sqrt(n()-1))

sfirst_fix_prop$block_type<-factor(sfirst_fix_prop$block_type,levels = c("fix","mix"),labels = c("Fixed","Mixed"))

Fix_Pro_Fig_heat <-  sfirst_fix_prop %>% filter(type!="Random",awareness=="Aware") %>%                   #,awareness=="Aware"
  ggplot(aes(x=type, y=mffp,  color= tar_pos_type, group= tar_pos_type)) +
  geom_point(position=pj, size=4) +  
  geom_line(position=pj) + 
  theme_bw() + 
  facet_wrap(~block_type) + 
  geom_errorbar(aes(ymin = mffp - se_ffp, ymax = mffp + se_ffp), width = 0.5,size=0.8,position=pj) + 
  labs(x="Fixated location", y="Proportion first fixations")+
  scale_color_manual(name ="Target position",values = c("#1B1B1B", "#696969","#A9A9A9"),breaks=c("Random", "Infreq. neighbour","Freq. neighbour"),labels = c("Random", "Infrequent", "Frequent"))+ #"#C77CFF", "#00BFC4","#F8766D"
  theme(legend.position = c(0.85,0.2),legend.title=element_blank())+ #scale_x_discrete(labels=c("Random"="Random","Infreq. neighbour"="Infrequent","Freq. neighbour"= "Frequent"))+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) 
   #'top'c("#000000", "#696969","#A9A9A9")


################################################## plot second fixation. ###############################

ssecond_fix_prop$block_type<-factor(ssecond_fix_prop$block_type,levels = c("fix","mix"),labels = c("fixed","mixed"))

second_Pro_Fig_heat <-  ssecond_fix_prop  %>% filter(type!="Random",awareness=="Aware") %>%                   #,awareness=="Aware"
  ggplot(aes(x=type, y=mffp,  color= tar_pos_type, group= tar_pos_type)) +
  geom_point(size=4) +  
  geom_line() + 
  theme_classic() + 
  facet_wrap(~block_type) + 
  geom_errorbar(aes(ymin = mffp - se_ffp, ymax = mffp + se_ffp), width = 0.2,size=0.8) + 
  labs(x="Fixated location", y="Proportion second fixations")+
  scale_color_manual(name ="Target position",values = c("#C77CFF", "#00BFC4","#F8766D"),breaks=c("Random", "Infreq. neighbour","Freq. neighbour"),labels = c("Random", "Infrequent", "Frequent"))+ #"#C77CFF", "#00BFC4","#F8766D"
  theme(legend.position = c(0.85,0.7),legend.title=element_blank()) #scale_x_discrete(labels=c("Random"="Random","Infreq. neighbour"="Infrequent","Freq. neighbour"= "Frequent"))

################ plot third fixation################

third_fix_prop$block_type<-factor(third_fix_prop$block_type,levels = c("fix","mix"),labels = c("fixed","mixed"))

third_Pro_Fig_heat <-  third_fix_prop  %>% filter(type!="Random",awareness=="Aware") %>%                   #,awareness=="Aware"
  ggplot(aes(x=type, y=mffp,  color= tar_pos_type, group= tar_pos_type)) +
  geom_point(size=4) +  
  geom_line() + 
  theme_classic() + 
  facet_wrap(~block_type) + 
  geom_errorbar(aes(ymin = mffp - se_ffp, ymax = mffp + se_ffp), width = 0.2,size=0.8) + 
  labs(x="Fixated location", y="Proportion third fixations")+
  scale_color_manual(name ="Target position",values = c("#C77CFF", "#00BFC4","#F8766D"),breaks=c("Random", "Infreq. neighbour","Freq. neighbour"),labels = c("Random", "Infrequent", "Frequent"))+ #"#C77CFF", "#00BFC4","#F8766D"
  theme(legend.position = c(0.85,0.7),legend.title=element_blank()) #scale_x_discrete(labels=c("Random"="Random","Infreq. neighbour"="Infrequent","Freq. neighbour"= "Frequent"))
################ plot cumulative fixations ################

cum_fix_prop <- cum_fix_prop %>% group_by(awareness,fix_loc,block_type, tar_pos_type,rep) %>%
  summarize(mffp = mean(ffp), se_ffp = sd(ffp)/sqrt(n()-1))

cum_fix_prop$block_type<-factor(cum_fix_prop$block_type,levels = c("fix","mix"),labels = c("Fixed","Mixed"))
cum_fix_prop$fix_loc<-factor(cum_fix_prop$fix_loc,levels = c("freq","infreq","repeat"),labels = c("Frequent","Infrequent","Repeated"))
cum_fix_prop$rep<-factor(cum_fix_prop$rep,levels = c("1.prop","2.prop","3.prop"),labels = c("1","2","3"))

cum_fix_prop_figure <- cum_fix_prop %>% filter(awareness=="Aware") %>%             
  ggplot(aes(x=rep, y=mffp, color= tar_pos_type, group= tar_pos_type)) +
  geom_point(size= 2.8) + 
  geom_line(size=0.8) + 
  theme_bw() + 
  scale_color_manual(name ="Target position",values = c("#000000", "#696969","#A9A9A9"),breaks=c("Random", "Infreq. neighbour","Freq. neighbour"),labels = c("Random", "Infrequent", "Frequent"))+ 
  geom_errorbar(aes(ymin = mffp - se_ffp, ymax = mffp + se_ffp), width = 0.2,size=0.8) + #,position=pj
  facet_nested(block_type~fix_loc) +
  labs(x="Fixation number", y="Cumulative proportions")+
  theme(legend.position = "top",legend.title=element_blank())+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) 


####PC effect based proportions of first fixations between frequent and infrequent condition, and explore
### correlation and confidence level(questionnaire)

awareness_data2 <-awareness_data2  %>% filter(awareness=="Aware") 

pc_Q1_eye<-awareness_data2%>% ggplot(aes(x=q1_rating, y=pc, color= awareness)) + 
  geom_point(size=1.8,alpha = 0.5 )+ 
  stat_poly_line(aes(fill= awareness))+
  stat_poly_eq() +
  labs(x="Q1 confidence rating", y="Probability-cueing effect")+
  theme_bw() +
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) 
pc_Q1_eye

pc_Q3_eye<-awareness_data2 %>% ggplot(aes(x=q3_resp, y=pc,color= awareness)) + 
  geom_point(size=1.8,alpha = 0.5 )+
  stat_poly_line(aes(fill= awareness))+
  stat_poly_eq()+ 
  labs(x="Q3 frequency rating", y="Probability-cueing effect")+
  theme_bw() + 
  theme(legend.position = "none")+ 
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(-0.5,0.7))+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9))
 
pc_Q3_eye


########intertrial priming  #########
first_fix_prop2 <-search_fixations %>% group_by(participant) %>% 
  mutate(previous_tar_pos_type = lag(tar_pos_type),
         dyn_dir=dyndirections[group], 
         f_t_dist = dists[dyn_dir*(which_fix - p_tar_pos)+8],
         fix_cond = ifelse(is.na(which_fix), "Non-search", 
                           ifelse(f_t_dist==0, "Repeat",
                                  ifelse(f_t_dist==1, "Freq.",
                                         ifelse(f_t_dist ==-1, "Infreq.", "Random"))))) %>%
  group_by(awareness, block_type, previous_tar_pos_type, participant) %>%   #looking at the target directly
  filter(fix_type!="Fixation") %>% filter(fix_count==min(fix_count)) %>% 
  filter(correct==1, tno %% 60 > 0, fix_cond!="Non-search") %>%
  summarize(repeated_prop = mean(fix_cond=="Repeat"), freq_prop = mean(fix_cond=="Freq."), infreq_prop=mean(fix_cond=="Infreq."),
            rand_prop=mean(fix_cond=="Random")) %>% 
  gather("type", "ffp", -c("awareness", "block_type", "previous_tar_pos_type", "participant")) 


first_fix_prop2$type <- factor(first_fix_prop2$type, levels=c("infreq_prop", "repeated_prop", "freq_prop", "rand_prop"), 
                               labels=c("Inreq.", "Repeated", "Freq.", "Random"), ordered = TRUE)

first_fix_prop3 <- group_by(first_fix_prop2,previous_tar_pos_type,type,awareness,block_type,participant) %>%  filter(previous_tar_pos_type!= 'Random',awareness=="Unaware") %>% summarize(mffp = mean(ffp)) %>% spread(type,mffp) %>% mutate( pc = `Freq.` - `Inreq.`)  

first_fix_prop4 <- group_by(first_fix_prop3,block_type,awareness,previous_tar_pos_type) %>%  filter(awareness== "Unaware" ,previous_tar_pos_type!= 'Random') %>%   
  summarize(mpc=mean(pc), sert=sd(pc)/sqrt(n()))

sfirst_fix_prop_new <-ggplot(first_fix_prop4, aes(x= previous_tar_pos_type, y=mpc,fill= previous_tar_pos_type)) + 
  geom_bar( stat="identity") + 
  scale_fill_grey(start = 0.4,end = 0.8)+ 
  geom_errorbar(aes(ymin = mpc - sert, ymax = mpc + sert), width=0.2,size=0.8)+ 
  theme_bw() + 
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) +
  theme(legend.title = element_blank())+
  facet_grid(~block_type)+
  labs(x="Previous target position", y="Mean probability-cueing effect") 

sfirst_fix_prop_new<-sfirst_fix_prop_new+scale_x_discrete(labels=c("Infreq. neighbour"="Infrequent","Freq. neighbour"= "Frequent"))
sfirst_fix_prop_new


######Proportion first fixations######

first_tar_sacc_latency<-first_tar_sacc_SE %>% filter( fix_cond!="Non-search",fix_cond!="Random",awareness=="Aware") %>%             
  ggplot(aes(x=fix_cond, y=mlatency, group= block_type)) + # color=block_type
  geom_point(aes(shape = block_type),size=3.4) +  
  geom_line(aes(linetype = block_type),size=0.8) + 
  theme_bw() + 
  geom_errorbar(aes(ymin = mlatency - se_latency, ymax = mlatency + se_latency), width = 0.2,size=0.8) + 
  labs(x="Fixated location", y="Latencies of initial saccades")+
  theme(legend.position = 'top',legend.direction = "horizontal")+
  theme(legend.title = element_blank())+
  theme(legend.position = c(0.61,0.1))+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9)) 


####exclude random
first_tar_sacc_new_unit<- first_tar_sacc_new %>% unite("fix_block",c(fix_cond,block_type))  
first_tar_sacc_new_unit<- spread(first_tar_sacc_new_unit,key= "fix_block",value = latency) 
first_tar_sacc_noblo <- first_tar_sacc  %>% group_by(fix_cond, awareness,  participant) %>% summarize(latency=mean(latency))

###distribution analyses of latency###

first_tar_distr <- first_tar_sacc %>% group_by(fix_cond, participant) %>%
  mutate(bin = bin(latency, nbins=9, labels = c("1","2","3","4","5","6","7","8","9"), method="content")) %>%
  filter(fix_cond %in% c("Freq.", "Infreq."),awareness=="Aware")  %>%  group_by(fix_cond, bin, participant) %>% 
  summarize(latency = mean(latency)) 

sfirst_tar_distr<- first_tar_distr %>% summarize(mlatency = mean(latency))
sfirst_tar_distr$fix_cond<-factor(sfirst_tar_distr$fix_cond,levels = c("Freq.","Infreq."),labels = c("Frequent","Infrequent"))

sfirst_tar_distr_plot<-sfirst_tar_distr %>% 
  ggplot(aes(x=bin, y=mlatency, color = fix_cond, group=fix_cond)) + 
  geom_point() + 
  geom_line() + 
  geom_hline(aes(yintercept=150),linetype=2) + 
  theme_classic()+
  theme(axis.title = element_text(face = "bold"))+
  theme(legend.position = "top",legend.title=element_blank())+
  labs(x="Vincentitles",y="First saccade latency (ms)")
sfirst_tar_distr_plot


#################### pattern #######

#fix_pattern both, one and neither

first_fix_prop_ppot <-  first_fix_prop_pattern2 %>%              
  ggplot(aes(x=type, y=mfp_prop,fill=type)) +
  geom_bar(stat="identity")+
  scale_fill_grey(start = 0.4,end = 0.8)+
  facet_wrap(~block_type) + 
  geom_errorbar(aes(ymin = mfp_prop - se_ffp, ymax = mfp_prop + se_ffp), width = 0.2,size=0.8) + 
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        legend.title=element_blank(),legend.text=element_text(face="bold",size=9)) +
  theme(strip.text = element_text(face="bold", size=9))+
  theme(legend.position = "none")+
  labs(x="Post-target saccades", y="Proportion of saccades")
first_fix_prop_ppot


