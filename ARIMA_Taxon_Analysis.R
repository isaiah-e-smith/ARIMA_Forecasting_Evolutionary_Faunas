# Isaiah E. Smith
# For Analytical Paleobiology Class Research Project
# Summer 2021

# This script follows ARIMA_Taxon_Preprocessing.R


##############################################################

library(mgcv)


# We'll use SQS data from "ARIMA_Taxon_Preprocessing.R"
tmp_div<-sqsDD$divRT

# Just in case any NA's remain, replace with 0's
tmp_div[is.na(tmp_div)] <- 0

# Let's take another look at the plot
tp(ylim=c(0,200), ylab="richness")	
lines(stages$mid, tmp_div, lwd=2)



# create empty data frame to store predictions for each interval
fo_storage<-data.frame(matrix(ncol =20, nrow = 10))

# and one to fill it with each loop
pdq<-data.frame(matrix(ncol =3, nrow = 10))

colnames(pdq)<-c("p", "d", "q")



# for-loop to give prediction at each interval 
# Also plots prediction for each loop, which makes a nice animated visual if saved as GIF

for (i in 1:85){ # only go to 85 to account for the 10 bin burn-in
  
  # In first loop, analyze first 10 intervals. Second loop: 11 intervals. Third loop: 12 intervals...
  a<-tmp_div[1:(i+9)] # currently using diversity
  b<-stages$mid[1:(i+9)]
  
  print(i)
  
  # construct generalized additive model
  ga<-gam(a ~ s(b)) 
  
  # fit ARIMA model to ga
  fit2 <- auto.arima(predict(ga))
  
  # forecast from the ARIMA model, with 20 forecasting periods
  fo <- forecast(fit2, h=20)
  
  plot(fo, ylim=c(-125, 125))
  
  # record average values 
  fo_storage<-rbind(fo_storage, fo$mean)
  
  # transpose the ARIMA order values
  ao<-t(data.frame(arimaorder(fit2)))
  
  # and save them to object pdq
  pdq<-rbind(pdq, ao)
  
}


# create data frame to record correlation values
cor_val<-data.frame(matrix(ncol =2, nrow = 10))

# and name the columns 
colnames(cor_val)<-c("estimate", "p.value")

for (i in 11:95){
  
  reg<-(95-i)
  
  if (reg>=19){ # compare 20 forecast values to corresponding 20 actual values
    
    cor<-cor.test((as.numeric(fo_storage[i,])), (tmp_div[i:(i+19)]))  
    
  } else{ # if less than 20 actual values remain in time series, compare the remaining actual values to the corresponding forecast values
    
    cor<-cor.test((as.numeric(fo_storage[i,1:reg])), (tmp_div[(i+1):(i+reg)]))
  }
  
  cor2<-as.numeric(t(cor[c("estimate", "p.value")]))
  
  # record results
  cor_val<-rbind(cor_val, cor2 )
}

# goes until 92, after which point there are not enough observations to assess correlation
# fill it out to maintain length
cor_val[93:95,]<-NA

# plot prediction for interval of interest, if desired
ex_int <- 34 # example interval selected
plot(as.numeric(fo_storage[ex_int,])) 


# create data frame to store values of interest
KPg_df_fo<-data.frame(matrix(ncol=6, nrow = nrow(fo_storage)))

colnames(KPg_df_fo)<-c("initial", "min", "max", "last", "set_mean", "stage_mean")

# use for loop to extract all values of interest
for (i in 1:nrow(fo_storage)){
  
  KPg_df_fo$initial[i]<-fo_storage[i,1] # first predicted interval
  
  KPg_df_fo$last[i]<-fo_storage[i,20] # last predicted interval
  
  KPg_df_fo$min[i]<-min(fo_storage[i,]) # minimum predicted interval
  
  KPg_df_fo$max[i]<-max(fo_storage[i,]) # maximum predicted interval
  
  KPg_df_fo$set_mean[i]<-mean(as.numeric(fo_storage[i,])) # average
}

#Calculate average predicted values for each stage

for (i in 11:95){
  
  reg<-(95-i)
  
  if (reg>=19){
    
    for (j in 1:19){
      
      tmp<-list()
      
      tmp<-append(tmp, as.numeric(fo_storage[(i-1)+j,j]))
      
      KPg_df_fo$stage_mean[i]<-mean(as.numeric(tmp))
    }
    
  } else {
    
    for (j in 1:reg){
      
      tmp<-list()
      
      tmp<-append(tmp, as.numeric(fo_storage[(i-1)+j,j]))
      
      KPg_df_fo$stage_mean[i]<-mean(as.numeric(tmp))
    }
  }
}



# normalize relative to zero, and adjust everything accordingly

low_val<-min(as.numeric(KPg_df_fo$min), na.rm = TRUE)*(-1)

KPg_df_fo<-KPg_df_fo+low_val


# add column to data frame which calculates proportion of change from
# first predicted point to lowest point, and first predicted point
# to last predicted point. Also calculate the log of these values.

KPg_df_fo$max_change_prop<-KPg_df_fo$min/KPg_df_fo$initial

KPg_df_fo$net_change_prop<-KPg_df_fo$last/KPg_df_fo$initial

KPg_df_fo$max_change_lgrt<-log(KPg_df_fo$min/KPg_df_fo$initial)

KPg_df_fo$net_change_lgrt<-log(KPg_df_fo$last/KPg_df_fo$initial)


# create empty column to fill in the next step
ngv_prop_remaining<-data.frame(matrix(ncol = 1, nrow = nrow(KPg_df_fo)))

KPg_df_fo$ngv_prop_remaining<-as.numeric(unlist(ngv_prop_remaining))


# What proportion of the remaining predictions predict a decrease in species abundance?
for (i in 1:nrow(KPg_df_fo)){
  
  y<-length(which(KPg_df_fo$net_change_lgrt[i:nrow(KPg_df_fo)]>=0))
  
  z<-length(which(KPg_df_fo$net_change_lgrt[i:nrow(KPg_df_fo)]<0))
  
  KPg_df_fo$ngv_prop_remaining[i]<-(z/(z+y))
  
}

# replace NA's with zeros
KPg_df_fo$net_change_lgrt[is.na(KPg_df_fo$net_change_lgrt)] <- 0

# combine columns into one data frame:
KPg_df_fo<-cbind(KPg_df_fo, cor_val)
KPg_df_fo<-cbind(KPg_df_fo, pdq)

# save to new object 
plot_dat<-KPg_df_fo




############ Exploratory Plotting #################

# Explore the data in different ways by experimenting with plotting 
# Nothing too serious, just messing around

plot(stages$mid, plot_dat$net_change_lgrt, ylim=c(-2,2), xlim=c(500,0))
plot(stages$mid, plot_dat$p, ylim=c(-2,2), xlim=c(500,0))
plot(stages$mid, plot_dat$d, ylim=c(-2,2), xlim=c(500,0))
plot(stages$mid, plot_dat$q, ylim=c(-2,2), xlim=c(500,0))
plot(stages$mid, plot_dat$estimate, ylim=c(0,1), xlim=c(500,0))
lines(stages$mid, plot_dat$p.value, ylim=c(1,0), col="red")
lines(stages$mid, plot_dat$ngv_prop_remaining, ylim=c(1,0), col="green")
points(stages$mid, plot_dat$q, ylim=c(1,0), col="blue")
plot(stages$mid, KPg_df_fo$ngv_prop_remaining, ylim=c(0,1), xlim=c(500,0))
tp(ylim=c(0,1.2), ylab="Proportion of remaining predictions showing decline")	
lines(stages$mid[14:95], plot_dat$ngv_prop_remaining[14:95], ylim=c(0,1.1), )


#x11()
par(mfrow=c(1,1)) 
tp(ylim=c(0,150), ylab="Number of Genera")	
lines(stages$mid, tmp_div, lwd=2)
tp(ylim=c(0,500), ylab="By-Stage Predicted Number of Genera")	
lines(stages$mid, plot_dat$set_mean)
h<-plot_dat$estimate

# Replace correlation values at or below 0.1 with NA
# is.na(h) <- h <= 0.1

desired_range <- (1:57) # change this based on taxon, target intervals, etc....

tp(ylim=c(0,200), ylab="By-Stage Average Predicted Number of Genera")	
lines(stages$mid[desired_range], plot_dat$stage_mean[desired_range], lwd = 2)
par(new=T)
plot(stages$mid, h, 
     type="p",axes=F, xlab=NA, ylab=NA,col="blue",
     ylim=c(0,1), pch=16, cex=1, legend = c("A", "B")
     , text.col = c("red", "blue")
     , pt.bg = c("red","blue")) # Removes axis and labels so they don't overlap
axis(side = 4)


# more plotting
tp(ylim=c(0,1.1), ylab="Prop. of Remaining Predictions that Show Declines")	
lines(stages$mid, as.numeric(unlist(KPg_df_fo$ngv_prop_remaining)))


par(new=T)
plot(stages$mid, plot_dat$estimate, type="p",axes=F, xlab=NA, ylab=NA,col="blue", ylim=c(0.,1)) 
axis(side = 4)
lines(stages$mid, plot_dat$p.value, ylim=c(1,0), col="red")


# Allow a second plot on the same graph
par(new=TRUE)
x11()
plot(stages$mid, tmp_div, ylim=c(0,100), xlim=c(500,0))
lines(stages$mid, tmp_div, col="black",  ylim=c(0,100), xlim=c(500,0))





##############################################################
# predict() with totally cumulative data

# start value should be earliest interval with diversity data for target taxon
# The first loop will then examine this bin and the following 10 bins to make its first prediction. 
# With each loop, it will have one more historic interval to inform the prediction of 200 future values

fore <- 200 # how many future intervals to predict
start <- 14 # earliest datum to be considered 

# add the appropriate number of empty rows ahead to keep things consistent
pre_storage<-data.frame(matrix(ncol = fore, nrow = (start+9)))



for (i in (start+10):81){  

  a<-sqsDD$divRT[start:i] # start will remain constant, but i will increase each loop
  
  b<-stages$mid[start:i]
  
  print(i)
  
  ga<-gam(a ~ s(b))
  
  pre <- predict(ga, data.frame(b=seq(1, fore)))
  
  pre_storage<-rbind(pre_storage, pre)
  
}

ex_int <- 70 # example interval selected

plot((1:fore), pre_storage[ex_int,], xlim=c(fore,1), ylim=c(0,100))




# create data frame to store values of interest

KPg_df_pre<-data.frame(matrix(ncol=3, nrow = nrow(pre_storage)))

colnames(KPg_df_pre)<-c("initial", "min", "max")


# use for loop to extract all values of interest
for (i in 1:nrow(pre_storage)){
  
  KPg_df_pre$initial[i]<-pre_storage[i,1]
  
  KPg_df_pre$min[i]<-min(pre_storage[i,])
  
  KPg_df_pre$max[i]<-max(pre_storage[i,])
}


# add column to data frame which calculates proportion of change from
# first predicted point to lowest point

KPg_df_pre$prediction<-KPg_df_pre$min/KPg_df_pre$initial





##############################################################
# predict() with moving window

# start value should be earliest interval with diversity data for target taxon
# The first loop will then examine this bin and the following 10 bins to make its first prediction. 
# With each loop, it will progress to the next 10 bin set to inform the prediction of 200 future values

fore <- 200 # how many future intervals to predict
start <- 14 # earliest datum to be considered 

# add the appropriate number of empty rows ahead to keep things consistent
pre_storage<-data.frame(matrix(ncol = fore, nrow = (start+9)))


for (i in (start+10):81){  
  
  a<-sqsDD$divRT[(i):(i+9)]
  
  b<-stages$mid[(i):(i+9)]
  
  print(i)
  
  ga<-gam(a ~ s(b))
  
  stages$mid
  
  pre <- predict(ga, data.frame(b=seq(stages$mid[i+10], stages$mid[i])))
  
  pre_storage<-rbind(pre_storage, pre)
  
}


# plot prediction for interval of interest, if desired

ex_int <- 70 # example interval selected

plot((1:fore), pre_storage[ex_int,], xlim=c(fore,1), ylim=c(0,100))



# create data frame to store values of interest
KPg_df_pre_win<-data.frame(matrix(ncol=3, nrow = nrow(pre_storage)))
colnames(KPg_df_pre_win)<-c("initial", "min", "max")


# use for loop to extract all values of interest
for (i in 1:nrow(pre_storage)){
  
  KPg_df_pre_win$initial[i]<-pre_storage[i,1]
  
  KPg_df_pre_win$min[i]<-min(pre_storage[i,])
  
  KPg_df_pre_win$max[i]<-max(pre_storage[i,])
  
}


# add column to data frame which calculates proportion of change from
# first predicted point to lowest point
KPg_df_pre_win$prediction<-KPg_df_pre_win$min/KPg_df_pre_win$initial











