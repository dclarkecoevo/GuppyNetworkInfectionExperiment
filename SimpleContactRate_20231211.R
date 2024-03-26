#Clear the working environment to be thorough
rm(list=ls())

#Read in the RDS file that has the processed and organized information for the videos of interest.
#For this video, I'm reading in GNIE_videos->all_XYdata->LAinf_Socinf_transmission.RDS, which is in the main GNIE_videos folder.
guppy_positions_uncorrected=readRDS(file.choose())

View(guppy_positions_uncorrected)
#Note that video 164 is a bit broken bc fish 8 is missing almost all of the time
#and fish 7 is also missing much of the time. There was no frame with all
#fish present. So let's drop that video.

View(guppy_positions_uncorrected)

#removing problematic videos that are not useful in contact rate calculation
#guppy_positions_uncorrected[1]<-NULL
# guppy_positions_uncorrected[[(1,2,3,4,5,6,7,8,9,10,11,12,13,41,35,238,239)]]<-NULL
#Removing all the videos specified below but the naming is messed up so have to index
guppy_positions_uncorrected[[1]]<-NULL
guppy_positions_uncorrected[[2]]<-NULL
guppy_positions_uncorrected[[3]]<-NULL
guppy_positions_uncorrected[[4]]<-NULL
guppy_positions_uncorrected[[5]]<-NULL
guppy_positions_uncorrected[[6]]<-NULL
guppy_positions_uncorrected[[7]]<-NULL
guppy_positions_uncorrected[[8]]<-NULL
guppy_positions_uncorrected[[9]]<-NULL
guppy_positions_uncorrected[[10]]<-NULL
guppy_positions_uncorrected[[11]]<-NULL
guppy_positions_uncorrected[[12]]<-NULL
guppy_positions_uncorrected[[13]]<-NULL
guppy_positions_uncorrected[[41]]<-NULL
guppy_positions_uncorrected[[35]]<-NULL
guppy_positions_uncorrected[[238]]<-NULL
guppy_positions_uncorrected[[239]]<-NULL

#Adding in detailed description of which videos to remove
#guppy_positions_uncorrected[["LAunf1_Day3_1200"]]<-NULL
#guppy_positions_uncorrected[["LAunf1_Day3_0930"]]<-NULL
#guppy_positions_uncorrected[["LAunf1_Day1_0930"]]<-NULL
#guppy_positions_uncorrected[["RCinf1_Day2_1200"]]<-NULL
#guppy_positions_uncorrected[["UYinf1_Day9_0930"]]<-NULL
#guppy_positions_uncorrected[["ADinf1_Day1_1200"]]<-NULL
#guppy_positions_uncorrected[["RWinf1_Day3_12_g"]]<-NULL
#guppy_positions_uncorrected[["RWinf1_Day3_4_gu"]]<-NULL
#guppy_positions_uncorrected[["RWinf1_Day3_930_"]]<-NULL
#guppy_positions_uncorrected[["RWinf1_Day4_12_g"]]<-NULL
#guppy_positions_uncorrected[["RWinf1_Day4_3_gu"]]<-NULL
#guppy_positions_uncorrected[["RWinf1_Day4_930_"]]<-NULL
#guppy_positions_uncorrected[["UYinf1_Day9_1200"]]<-NULL
#guppy_positions_uncorrected[["UYinf1_Day9_1500"]]<-NULL




#Read in the file that has info on all the videos about pop, day, num_fish,
#I usually call this vid_summaries or something similar.
vid_info=readRDS(file.choose())
View(vid_info)
#Converting between actual internal width of the cage, 59.5 cm, to what Trex thinks it is, 24.6 cm
conversion_factor=59.5/24.6

#Get some parameters like the number of videos and the number of fish. For now, I'm assuming the number of fish is the same in
#all videos
num_videos=length(guppy_positions_uncorrected)

#Removing Inf and replacing with NA
#for(h in 1:num_videos){
  #guppy_positions_uncorrected[[h]][which(!is.finite(guppy_positions_uncorrected[[h]]))]<-NA
  
#}

#View(guppy_positions)
#View(guppy_positions_uncorrected)
#Correct all guppy positions by the conversion_factor
guppy_positions=list()

for(v in 1:num_videos){
  num_fish=(dim(guppy_positions_uncorrected[[v]])[2]-1)/2
  #Copy guppy_positions_uncorrected
  guppy_positions[[v]]=guppy_positions_uncorrected[[v]]
  
  #Then multiply all positions by the conversion factor
  guppy_positions[[v]][,(2:(1+2*num_fish))]=guppy_positions_uncorrected[[v]][,(2:(1+2*num_fish))]*conversion_factor
}

#Get the distance between all pairs of guppies. This takes a while, e.g, 8 sec
#per video.
guppy_distances=list()
for(v in 1:num_videos){
  num_fish=(dim(guppy_positions_uncorrected[[v]])[2]-1)/2
  #There are a small percentage of examples where guppy_positions are Inf bc the position was lost
  #for a frame. Overwrite these with NAs to make sure they don't factor into the computation.
  guppy_positions[[v]][which(is.infinite(guppy_positions[[v]]))]<-NA
  
  #Calculate distance between any two fish
  guppy_distances[[v]]=array(NA,dim=c(dim(guppy_positions[[v]])[1],num_fish,num_fish))
  for(row_number in 1:dim(guppy_positions[[v]])[1]){
    for(guppya_number in 1:num_fish){
      for(guppyb_number in 1:num_fish){
        #Calculate distance for each guppy pair and time as sqrt((xa-xb)^2+(ya-yb)^2)
        guppy_distances[[v]][row_number,guppya_number,guppyb_number]=sqrt(sum((guppy_positions[[v]][row_number,(2:3)+2*(guppya_number-1)]-guppy_positions[[v]][row_number,(2:3)+2*(guppyb_number-1)])^2))
      }
    }
  }
  print(v)
}

#Check that the names make sense and then keep the names across the lists.
names(guppy_positions_uncorrected)
names(guppy_positions)<-names(guppy_positions_uncorrected)
names(guppy_distances)<-names(guppy_positions)

#Inputting fisha index, fishb index, the chosen contact distance, and the
#matrix of fish distances, calculate the number of contacts between fish a
#and fish b.
frequency_times=function(fisha,fishb,contact_distance,input_list){
  num_videos=length(input_list)
  frequency_output=list()
  temp_diffs=list()
  for(v in 1:num_videos){
    #Count 0, NA x times, 0 y times -1 as 1
    #For a given -1, is there a NA more recently than a 1?
    #If so, was there a 0 immediately before that set of NAs?
    temp_diffs[[v]]=NA
    temp_diffs[[v]]=diff(input_list[[v]][,fisha,fishb]<contact_distance)
    
    #When did the fish leave move apart outside the contact distance?
    indices_neg1=which(temp_diffs[[v]]==(-1))
    
    #Need to count times when a distance was at least one zero before an NA,
    #any number of NAs and zeros, then a -1. This means that the fish entered
    #contact distance without being noticed and has now left. Technically,
    #we could miss an NA then fish entering contact distance then staying there
    #until the end but this should be negligible as long as the frequency of NAs is small
    neg1_contribution=0
    if(length(indices_neg1)>0){
      for(j in 1:length(indices_neg1)){
        #Default neg1_status is 0
        neg1_status=0
        index_neg1=indices_neg1[j]
        
        #If it's at the very beginning, count it as a contact
        if(index_neg1==1){
          neg1_status=1
        }else{
          #Starting from right before the -1 index and moving backward.
          #With the break statements, this most nested loop scans 
          #backwards until it either finds a 1 leading to the -1 or finds
          #a 0, NA sequence. In the latter case, this means that there was no
          #1 before the -1 but after a 0, NA sequence.
          for(i in (index_neg1-1):1){
            if(!is.na(temp_diffs[[v]][i])&&temp_diffs[[v]][i]==1){
              #If the previous index had a 1, then stop this most nested loop. 
              #This index_neg1 doesn't contribute anything.
              break
            }
            if(!is.na(temp_diffs[[v]][i])&&temp_diffs[[v]][i]==0&&is.na(temp_diffs[[v]])[i+1]){
              #If the previous index was a zero right before an NA, then this
              #index_neg1 does contribute to the number of contacts. 
              neg1_status=1
              break
            }
          }
        }
        #Keep a running total of neg1_contribution.
        neg1_contribution=neg1_contribution+neg1_status
      }
    }
    
    #Count starting TRUE as 1. Count 1 as 1. And count 0, X NAs, Y 0s, -1 as 1
    frequency_output[[v]]=length(which(input_list[[v]][1,fisha,fishb]<contact_distance))+length(which(temp_diffs[[v]]==1))+neg1_contribution
  }
  return(frequency_output)
}

#Contacts per minute, averaged across videos
con_per_min=function(fisha,fishb,contact_distance,input_list){
  num_frames=0
  for(v in 1:length(input_list)){
    num_frames=num_frames+length(which(!is.na(input_list[[v]][,fisha,fishb])))
  }
  #This assumes a frame rate of 20 fps for all videos.
  return(sum(unlist(frequency_times(fisha,fishb,contact_distance,input_list))))
}


#Read in files of worm counts and ID_trajectories
worm_counts=read.csv(file.choose())
ID_trajectories=read.csv(file.choose())
GZID<-subset(ID_trajectories, Population == "GZinf1")
#UYID<-subset(ID_trajectories, Population == "UYinf1")
UGID<-subset(ID_trajectories, Population == "UGinf1")
#For now, work with a distance of 2 cm as this seems to be in the ballpark of
#correct.
distance=3.5

#Get rate of contact between overallID X and overallID Y on day(s) Z for Chosen_population
#Note that males should always have overallIDs 1-3 while the index gets 10
#By default, this returns one value for each day put into the days argument
#Note that this function will mess up if one of the fish is the index (10)
#and you include days before 3. Because the index isn't present then.
dyad_contact_rate=function(overallIDX,overallIDY,days,Chosen_population){
  num_days=length(days)
  contact_rate=rep(NA,length(days))
  for(d in 1:length(days)){
    #You have to intelligently grab the Trex ID for the chosen fish and increment
    #it by 1 to be an index since Trex IDs start at 1.
    X_vid_index=unique(ID_trajectories[which(ID_trajectories$Day==days[d]&ID_trajectories$Population==Chosen_population&ID_trajectories$IDtype=="Trex"),overallIDX+4])+1
    Y_vid_index=unique(ID_trajectories[which(ID_trajectories$Day==days[d]&ID_trajectories$Population==Chosen_population&ID_trajectories$IDtype=="Trex"),overallIDY+4])+1
    contact_rate[d]=con_per_min(#Input the correct IDs for the two fish on those days
      X_vid_index,Y_vid_index,distance,
      #Input the part of the guppy_distances object that corresponds
      #to the videos you want based on population and day
      guppy_distances[which(vid_info[,1]==Chosen_population&vid_info[,2]==days[d])])
  }
  return(contact_rate)
}

#A vector that helps count up which dyads go with which fish
index_list=list(1:9,10:17,18:24,25:30,31:35,36:39,40:42,43:44,45)

#This automatically incorporates how many fish are present on each day.
population_contact_rate=function(days,Chosen_population){
  #Make an array to store the results for the all the different dyads
  #These columns are populated with the 9 contact rates for the dyads
  #including fish with overall ID 1, then the 8 other dyads including fish
  #with overall ID 2 but not the dyad already counted of fish 1 and 2, etc.
  #Contact rate is NA for dyads involving the index when the index
  #is not yet present (i.e., days before 3)
  
  #Then the 46 column is the average across all dyads.
  all_dyads=array(NA,dim=c(length(days),46))
  for(d in 1:length(days)){
    #Do the 9 fish version if it's before day 3
    if(days[d]<3){
      #Loop through fish1 identities. There is no need to go up to 10
      #because all dyads including the 10th fish will already be accounted for.
      for(f1 in 1:8){
        #Loop through fish2 identities. Only consider fish2 having a higher
        #ID than fish 1.
        for(f2 in (f1+1):9){
          #Populate the correct columns for each dyad
          all_dyads[d,index_list[[f1]][f2-f1]]=dyad_contact_rate(f1,f2,days[d],Chosen_population)
        }
      }
      
    }else{
      #Else do the 10 fish version
      #Loop through fish1 identities. There is no need to go up to 10
      #because all dyads including the 10th fish will already be accounted for.
      for(f1 in 1:9){
        #Loop through fish2 identities. Only consider fish2 having a higher
        #ID than fish 1.
        for(f2 in (f1+1):10){
          #Populate the correct columns for each dyad
          all_dyads[d,index_list[[f1]][f2-f1]]=dyad_contact_rate(f1,f2,days[d],Chosen_population)
        }
      }
    }
    #Then the last column is the average of all dyads for that population and that day
    all_dyads[d,46]=sum(all_dyads[d,which(!is.na(all_dyads[d,]))])/length(which(!is.na(all_dyads[d,])))
  }
  return(all_dyads)
}



GZinf1=data.frame(population_contact_rate(1:3,"GZinf1"))
GZinf1$pop<-"GZinf1"
UYinf1=data.frame(population_contact_rate(1:3,"UYinf1"))
UYinf1$pop<-"UYinf1"
UGinf1=data.frame(population_contact_rate(1:3,"UGinf1"))
UGinf1$pop<-"UGinf1"
GZunf1=data.frame(population_contact_rate(1:3,"GZunf1"))
GZunf1$pop<-"GZunf1"
GDinf1=data.frame(population_contact_rate(1:3,"GDinf1"))#Problem with this one
GDinf1$pop<-"GDinf1"#Problem with this one
LAunf1=data.frame(population_contact_rate(1:3,"LAunf1"))
LAunf1$pop<-"LAunf1"
RCinf1=data.frame(population_contact_rate(2:3,"RCinf1"))
RCinf1$pop<-"RCinf1"
TRinf1=data.frame(population_contact_rate(1:3,"TRinf1"))#Problem
TRinf1$pop<-"TRinf1"#Problem
SBunf1=data.frame(population_contact_rate(1:3,"SBunf1"))
SBunf1$pop<-"SBunf1"
ADinf1=data.frame(population_contact_rate(1:3,"ADinf1"))
ADinf1$pop<-"ADinf1"

Popcontacts<-rbind(ADinf1,SBunf1,LAunf1,RCinf1,UGinf1,UYinf1,GZinf1,GZunf1)

write.csv(Popcontacts, "20230221_Popcontactsonly.csv")
write.csv(LAunf1, "20231213_LAunfcontacts.csv")
write.csv(ADinf1, "20231206_LAunfcontacts.csv")

#f1 populates the first 9, f2 populates the next 8, f3 populates the next 7, etc.
#This returns all dyads in the first 45 columns and the average in the 46 column
test_colors=colorRampPalette(c("purple","red"))
test_colors(45)

test=GDinf1


plot(1:3,1:3+NA,type="l",lwd=3,col="black",ylim=range(test,na.rm=T),xlab="Day",ylab="Contacts per minute per dyad",main="GZinf1")
for(i in 1:45){
  points(1:3,test[,i],type="l",lwd=1,lty=2,col=test_colors(45)[i])
}
points(1:3,test[,46],type="l",lwd=3,col="black",ylim=range(test,na.rm=T),xlab="Day",ylab="Contacts per minute per dyad",main="GZinf1")

ID_trajectories3<-subset(ID_trajectories, Day=="3")
View(ID_trajectories3)
PopulationContacts<-read.csv(file.choose())

library(ggplot2)
library(visreg)
library(lme4)
library(glmmTMB)


ggplot(PopulationContacts, aes(y=Day3PopAvg, x=Treatment))+geom_boxplot()


ggplot(PopulationContacts, aes(y=PopulationContacts$Inf., x=Day2Popavg))+geom_point()+geom_smooth()+facet_wrap(~Sex)


lmPopcon<-lm(Day3PopAvg~Treatment, PopulationContacts)
summary(lmPopcon)
visreg(lmPopcon)

lmPopconprev<-glm(prev~Day3PopAvg, family=binomial, PopulationContacts)
summary(lmPopconprev)
visreg(lmPopconprev)

lmIndcon<-lm(Day3ContactRate~Treatment*Sex, PopulationContacts)
summary(lmIndcon)
visreg(lmIndcon)

lmIndcon2<-lm(Day2ContactRate~Treatment*Sex, PopulationContacts)
summary(lmIndcon2)
visreg(lmIndcon2)

lmIndconInf<-glmmTMB(Inf.~Day3ContactRate*Sex+(1|Pop), family=betabinomial(), PopulationContacts)
summary(lmIndconInf)
visreg(lmIndconInf, scale='response')

PopConIndex<-subset(PopulationContacts, Index=="1")
lmIndexCon<-glm(Day3ContactRate~Treatment, family=Gamma(link="identity"), PopConIndex)
summary(lmIndexCon)
visreg(lmIndexCon, scale="response")
