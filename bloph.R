#Written by Garima Gupta May 2017.


#remove anything from global environment
rm(list=ls())

#create a list of packages that we might want to install
packages<-c("ggplot2","plyr","reshape2","dplyr","pander","captioner","bbmle","broom","tidyr","gridExtra","rgdal","sp","maptools","raster","scales","ape","rgeos",
            "maps","mapdata","spatstat","GISTools","gdalUtils","MODISTools","MODIS","rts","RCurl","data.table")

#install any CRAN packages if not already installed
inst <- packages %in% installed.packages()
if(length(packages[!inst]) > 0) install.packages(packages[!inst])

#use lapply to loop through all packages and load them
lapply(packages, require, character.only = TRUE)

#set options for dplyr so can see tables properly
options(dplyr.width = Inf,dplyr.length=Inf)

#set the ggplot2 theme to be theme_bw() globally
theme_set(theme_bw())

#create a special version of length that handles NAs
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)
}

#set wd - change to wherever your shapefiles are stored
setwd("H:/Jon_material/R_GIS")


######################################################### MAIN CODE STARTS HERE#######################################

#Step 1: ctrl + Enter to read these functions I've written

#function for counterexamples
fn.counterexample<-function(shapefile,ni){
  
  #store the time we started running the function
  started.at=proc.time()
  
  #create a new filename depending on species
  bloph<-unique(levels(factor(shapefile$SpcsC)))
  filename.year <- paste(bloph, "_year_c.csv", sep="") 
  filename.count <- paste(bloph, "_count_c.csv", sep="") 
  filename.graphs <- paste(bloph, "_graphs_c.png", sep="") 
  
  #remove rows with na in year column as it messes things up if not removed
  shapefile<-subset(shapefile, !is.na(YrFrA))
  
  #get each unique year
  n.years<-sort(unique((shapefile$YrFrA)))
  
  #initialise empty vector of area
  empty.area<-NULL
  
  #initialise empty vector of count
  empty.count<-NULL
  
  #project shapefile
  shapefile.p<-spTransform(shapefile,CRS("+proj=moll"))
  
  #run a loop to repeat ni times
  for (i in 1:ni){
    
    #create a random year column
    shapefile.p$year.ran<-sample(shapefile.p$YrFrA,replace=F)
    #store those years
    #years.random <- sort(unique(shapefile$year.ran))
    
    #split the shapefile
    dp.list<-lapply(sort(unique(shapefile.p$year.ran)), function(i) shapefile.p[shapefile.p$year.ran <= i, ])
    
    #calculate the mcp area for each part of the list
    dp.list.mcp<-lapply(dp.list,function(x) gConvexHull(x,byid=F,id=NULL))
    
    #calculate area of mcp for each part of list
    dp.list.mcp.area<-sapply(dp.list.mcp,function(x) gArea(x,byid=F))
    
    #convert area list into dataframe
    dp.list.mcp.area.df<-as.data.frame(dp.list.mcp.area)
    
    #add each area dataframe to a list
    empty.area[i]<-dp.list.mcp.area.df
    
    #calculate the nrow for each part of the list
    dp.list.count<-sapply(dp.list,function(x) nrow(x))
    
    #convert list to a dataframe
    dp.list.count.df<-as.data.frame(dp.list.count)
    
    #add each row count to a list
    empty.count[i]<-dp.list.count.df
    
    #inform us about progress of loop
    cat('Processed iteration', i, 'of', ni,'at',paste(Sys.time()),'\n')
  }
  
  #once loop finished
  
  #convert every item in the mcp area list to a dataframe
  df.area<-data.frame(matrix(unlist(empty.area), ncol=1, byrow=T))
  
  #create a vector called years where each sequence of unique years is repeated ni times
  years<-rep(n.years,ni)
  
  #convert the vector to a df
  df.years<-as.data.frame(years)
  
  #convert every item in the count list to a dataframe
  df.count<-as.data.frame(matrix(unlist(empty.count),ncol=1,byrow=T))
  
  #df.count<-as.data.frame(dp.list.count)
  
  #combine the three dataframe columns
  df.final<-cbind(df.years,df.area,df.count)
  
  #rename the columns
  names(df.final)[1]<-"year"
  names(df.final)[2]<-"mcp_area"
  names(df.final)[3]<-"count"
  
  #get some summary statistics for area
  df.summary.area<-df.final%>%
    group_by(year)%>%
    summarise(mean.area=mean(mcp_area),se.area=sd(mcp_area)/sqrt(length2(mcp_area)))
  
  #get some summary statistics for count
  df.summary.count<-df.final%>%
    group_by(count)%>%
    summarise(mean.area=mean(mcp_area),se.area=sd(mcp_area)/sqrt(length2(mcp_area)))
  
  #export the summary data
  write.csv(df.summary.area,file=filename.year,row.names=F)
  write.csv(df.summary.count,file=filename.count,row.names=F)
  
  #create two graphs and export
  pd<-position_dodge(.1)
  
  g.area<-ggplot(df.summary.area,aes(x=as.numeric(year),y=mean.area/1000))+geom_point()+geom_path()+
    geom_errorbar(aes(ymin=mean.area/1000-se.area/1000, ymax=mean.area/1000+se.area/1000), width=0.1, position=pd,linetype="solid")+
    xlab("Year")+ylab("log10(EOO area sq-km)")+stat_smooth(method="loess",se=F,colour="red")+
    scale_y_log10()+ggtitle(paste(bloph,"counterexample year: iterations = ",ni,"\n"))
  
  g.count<-ggplot(df.summary.count,aes(x=as.numeric(count),y=mean.area/1000))+geom_point()+geom_path()+
    geom_errorbar(aes(ymin=mean.area/1000-se.area/1000, ymax=mean.area/1000+se.area/1000), width=0.1, position=pd,linetype="solid")+
    xlab("Record count")+ylab("log10(EOO area sq-km)")+stat_smooth(method="loess",se=F,colour="red")+
    scale_y_log10()+ggtitle(paste(bloph,"counterexample count: iterations = ",ni,"\n"))
  
  g.combine<-grid.arrange(g.area,g.count,ncol=1)
  
  ggsave(file=filename.graphs,g.combine,dpi=300,height=6,width=6)
  
  #inform us the process has finished
  cat('Finished processing',ni, 'iterations in', timetaken(started.at),"\n")
  
}

#function for examples
fn.example<-function(shapefile){
  
  #store the time we started running the function
  started.at=proc.time()
  
  #create a new filename depending on species
  bloph<-unique(levels(factor(shapefile$SpcsC)))
  filename <- paste(bloph, "_example.csv", sep="") 
  filename.graphs <- paste(bloph, "_example_graphs.png", sep="") 
  
  #remove rows with na in year column as it messes things up if not removed
  shapefile<-subset(shapefile, !is.na(YrFrA))
  
  #get each unique year
  n.years<-sort(unique((shapefile$YrFrA)))
  
  #initialise empty vector of area
  empty.area<-NULL
  
  #initialise empty vector of count
  empty.count<-NULL
  
  #project shapefile
  shapefile.p<-spTransform(shapefile,CRS("+proj=moll"))
  
  #split the shapefile
  dp.list<-lapply(sort(unique(shapefile.p$YrFrA)), function(i) shapefile.p[shapefile.p$YrFrA <= i, ])
  
  #calculate the mcp area for each part of the list
  dp.list.mcp<-lapply(dp.list,function(x) gConvexHull(x,byid=F,id=NULL))
  
  #calculate area of mcp for each part of list
  dp.list.mcp.area<-sapply(dp.list.mcp,function(x) gArea(x,byid=F))
  
  #convert area list into dataframe
  dp.list.mcp.area.df<-as.data.frame(dp.list.mcp.area)
  
  #calculate the nrow for each part of the list
  dp.list.count<-sapply(dp.list,function(x) nrow(x))
  
  #convert list to a dataframe
  dp.list.count.df<-as.data.frame(dp.list.count)
  
  #convert the vector to a df
  df.years<-as.data.frame(n.years)
  
  #combine the three dataframe columns
  df.final<-cbind(df.years,dp.list.mcp.area.df,dp.list.count.df)
  
  #rename the columns
  names(df.final)[1]<-"year"
  names(df.final)[2]<-"mcp_area"
  names(df.final)[3]<-"count"
  
  
  #export the summary data
  write.csv(df.final,file=filename,row.names=F)
  
  #create two graphs and export
  g.area<-ggplot(df.final,aes(x=as.numeric(year),y=mcp_area/1000))+geom_point()+
    geom_path()+
    xlab("Year")+ylab("log10(EOO area sq-km)")+stat_smooth(method="loess",se=F,colour="red")+
    scale_y_log10()+ggtitle(paste(bloph,"year","\n"))
  
  g.count<-ggplot(df.final,aes(x=as.numeric(count),y=mcp_area/1000))+geom_point()+
    geom_path()+
    xlab("Record count")+ylab("log10(EOO area sq-km)")+stat_smooth(method="loess",se=F,colour="red")+
    scale_y_log10()+ggtitle(paste(bloph,"count","\n"))
  
  g.combine<-grid.arrange(g.area,g.count,ncol=1)
  
  ggsave(file=filename.graphs,g.combine,dpi=300,height=6,width=6)
  
  #inform us the process has finished
  cat('Finished processing in', timetaken(started.at),"\n")
  
}
#Step 2: read in the shapefile for the species you want (make sure you have exported it by splitting our shapefile in R as before)
bloph<-readOGR(dsn=".",layer="bloph") #note that the number of features may be higher than the record count we eventually see
#because some of these locality records will not have a year associated with them or lat/lon information associated with them

#Step 3: give it a geographic projection if it doesn't already have one - our function will do the projected coordinate
#system stuff, so don't worry about it for now
proj4string(bloph) <- CRS("+proj=longlat +ellps=WGS84")

#Step 4: plot records
plot(bloph)

#Step 5: check there are no records in obviously wrong countries

xtabs(~Cntry,scmph)

#Step 6: remove records in wrong country

scmph<-subset(scmph,Cntry=="India"|Cntry=="Nepal"|Cntry=="China"|Cntry=="Bhutan")


xtabs(~Cntry,scmph)


#Step 7: use the functions.  In fn.example, just write the name of your species' shapefile in the brackts
#In fn.counterexample, write the species shapefile AND the number of iterations you want.

#NOTE: the processing time of the counterexample function increases by approximately 10x for every 10x iteration.  i.e. 
#for himph, if you have it at 10, it takes 6.56s. If you have it at 100, it takes 40s.  I predict it will take 400s for 1000
#and 4000 for 10000.  Thus ~7 minutes for 1000 and 66 minutes for 10,000.  Obviously with species with larger ranges
#and more records, the longer it will take.  I'd run the fn.counterexample with just 10 iterations at first, to get an idea
#of how long we would need to run everything. e.g. quail of chuka might be the species that take the longest to run.

fn.example(bloph)
fn.counterexample(bloph,1000)

#######finding a plateau#####

###example year###
data1<-read.csv("bloph_example.csv")
head(data1)
(graph1<-ggplot(data1,aes(x=year,y=mcp_area))+geom_point()+stat_smooth(method="loess")+geom_line())+xlab("Year")+ylab("Area")
fn<-function(x){
  
  #find nrow
  number_row<-nrow(x)
  
  #find mcp area at max year
  max.mcp_area<-x[number_row,2]
  
  #find 5% of mcp area of max year
  max.mcp_area_5<-(max.mcp_area/100)*0.05
  
  #reverse order of rows so when we calculate our cumulative average we start with the latest year not the earliest
  xr<-x[rev(rownames(x)),]
  
  #add a column with the cumulative average
  xr$cum_ave<-cumsum(xr$mcp_area)/seq_along(xr$mcp_area) #if using counterexamples, etc, change mcp_area to whatever you are looking at
  
  #add a column with difference between max.mcp_area and cum_ave
  xr$diff<-max.mcp_area-xr$cum_ave
  
  #now subset to ensure we only have rows where the difference is less than
  #say 0.05% of the maximum MCP area
  
  xr.diff<-subset(xr,diff<max.mcp_area_5)
  
  #find the number of years that are 'stable'
  xr.diff.nyear<-max(xr.diff$year)-min(xr.diff$year)
  
  #now find the minimum year to say when plateau began
  xr.diff.min<-min(xr.diff$year)
  
  #print year where we say plateau starts and how many years are stable
  cat('Plateau starts at', xr.diff.min,"and lasts for",xr.diff.nyear,"years")
  
  #create graph with vertical line showing start of plateau
  ggplot(x,aes(x=year,y=mcp_area))+geom_point()+geom_line()+
    geom_vline(xintercept=xr.diff.min,col="red")
  
}

#use function with dataframe
fn(data1)

###example_count####
data2<-read.csv("bloph_example.csv")
head(data2)
(graph2<-ggplot(data2,aes(x=count,y=mcp_area))+geom_point()+stat_smooth(method="loess")+geom_line())+xlab("Count")+ylab("Area")

fn.count<-function(x){
  
  #find nrow
  number_row<-nrow(x)
  
  #find mcp area at max year
  max.mcp_area<-x[number_row,2]
  
  #find 5% of mcp area of max year
  max.mcp_area_5<-(max.mcp_area/100)*0.05
  
  #reverse order of rows so when we calculate our cumulative average we start with the latest year not the earliest
  xr<-x[rev(rownames(x)),]
  
  #add a column with the cumulative average
  xr$cum_ave<-cumsum(xr$mcp_area)/seq_along(xr$mcp_area) #if using counterexamples, etc, change mcp_area to whatever you are looking at
  
  #add a column with difference between max.mcp_area and cum_ave
  xr$diff<-max.mcp_area-xr$cum_ave
  
  #now subset to ensure we only have rows where the difference is less than
  #say 0.05% of the maximum MCP area
  
  xr.diff<-subset(xr,diff<max.mcp_area_5)
  
  #find the number of records that are 'stable'
  xr.diff.ncount<-max(xr.diff$count)-min(xr.diff$count)
  
  #now find the minimum record count to say when plateau began
  xr.diff.min<-min(xr.diff$count)
  
  #print record where we say plateau starts and how many records are stable
  cat('Plateau starts at', xr.diff.min,"records and lasts for",xr.diff.ncount,"records")
  
  #create graph with vertical line showing start of plateau
  ggplot(x,aes(x=count,y=mcp_area))+geom_point()+geom_line()+
    geom_vline(xintercept=xr.diff.min,col="red")
  
}

fn.count(data2)

######for counter_example_year#####
data3<-read.csv("bloph_year_c.csv")
head(data3)
(graph3<-ggplot(data3,aes(x=year,y=mean.area))+geom_point()+stat_smooth(method="loess")+geom_line())+xlab("Year")+ylab("Area")

fn<-function(x){
  
  #find nrow
  number_row<-nrow(x)
  
  #find mcp area at max year
  max.mean.area<-x[number_row,2]
  
  #find 5% of mcp area of max year
  max.mean.area_5<-(max.mean.area/100)*0.05
  
  #reverse order of rows so when we calculate our cumulative average we start with the latest year not the earliest
  xr<-x[rev(rownames(x)),]
  
  #add a column with the cumulative average
  xr$cum_ave<-cumsum(xr$mean.area)/seq_along(xr$mean.area) #if using counterexamples, etc, change mcp_area to whatever you are looking at
  
  #add a column with difference between max.mcp_area and cum_ave
  xr$diff<-max.mean.area-xr$cum_ave
  
  #now subset to ensure we only have rows where the difference is less than
  #say 0.05% of the maximum MCP area
  
  xr.diff<-subset(xr,diff<max.mean.area_5)
  
  #find the number of years that are 'stable'
  xr.diff.nyear<-max(xr.diff$year)-min(xr.diff$year)
  
  #now find the minimum year to say when plateau began
  xr.diff.min<-min(xr.diff$year)
  
  #print year where we say plateau starts and how many years are stable
  cat('Plateau starts at', xr.diff.min,"and lasts for",xr.diff.nyear,"years")
  
  #create graph with vertical line showing start of plateau
  ggplot(x,aes(x=year,y=mean.area))+geom_point()+geom_line()+
    geom_vline(xintercept=xr.diff.min,col="red")
  
}

#use function with dataframe
fn(data3)


#########counter_example_count#####
data4<- read.csv("bloph_count_c.csv")
head(data4)
(graph4<-ggplot(data4,aes(x=count,y=mean.area))+geom_point()+stat_smooth(method="loess")+geom_line())+xlab("Count")+ylab("Area")


fn.count<-function(x){
  
  #find nrow
  number_row<-nrow(x)
  
  #find mcp area at max year
  max.mean.area<-x[number_row,2]
  
  #find 5% of mcp area of max year
  max.mean.area_5<-(max.mean.area/100)*0.05
  
  #reverse order of rows so when we calculate our cumulative average we start with the latest year not the earliest
  xr<-x[rev(rownames(x)),]
  
  #add a column with the cumulative average
  xr$cum_ave<-cumsum(xr$mean.area)/seq_along(xr$mean.area) #if using counterexamples, etc, change mcp_area to whatever you are looking at
  
  #add a column with difference between max.mcp_area and cum_ave
  xr$diff<-max.mean.area-xr$cum_ave
  
  #now subset to ensure we only have rows where the difference is less than
  #say 0.05% of the maximum MCP area
  
  xr.diff<-subset(xr,diff<max.mean.area_5)
  
  #find the number of records that are 'stable'
  xr.diff.ncount<-max(xr.diff$count)-min(xr.diff$count)
  
  #now find the minimum record count to say when plateau began
  xr.diff.min<-min(xr.diff$count)
  
  #print record where we say plateau starts and how many records are stable
  cat('Plateau starts at', xr.diff.min,"records and lasts for",xr.diff.ncount,"records")
  
  #create graph with vertical line showing start of plateau
  ggplot(x,aes(x=count,y=mean.area))+geom_point()+geom_line()+
    geom_vline(xintercept=xr.diff.min,col="red")
  
}

fn.count(data4)


######merging two plots- year###
#read the dataframes you need
d1<-read.csv("bloph_example.csv")
d2<-read.csv("bloph_year_c.csv")
d3<-read.csv("bloph_count_c.csv")
head(d1)
head(d2)
head(d3)

#let's add some extra columns so we can differentiate between dataframes
#when we merge them later
d1$type<-"example"
d2$type<-"counterexample"
d3$type<-"counterexample"

#let's also call the area columns by the same thing
names(d1)[2]<-"area"
names(d2)[2]<-"area"
names(d3)[2]<-"area"


################One way of doing it###############
#I would probably do the following
d1.sub<-subset(d1,select=c("year","area","type"))
d2.sub<-subset(d2,select=c("year","area","type"))

d1.subc<-subset(d1,select=c("count","area","type"))
d3.sub<-subset(d3,select=c("count","area","type"))

#now use rbind
dy.combine<-rbind(d1.sub,d2.sub)
head(dy.combine)

dc.combine<-rbind(d1.subc,d3.sub)
head(dc.combine)

#then add in the se column from d2
d2.se<-subset(d2,select=c("year","se.area","type"))
head(d2.se)

#do same for d3 for the count stuff
d3.se<-subset(d3,select=c("count","se.area","type"))
head(d3.se)
#now merge
dy<-merge(dy.combine,d2.se,by=c("year","type"),all.x=T)
head(dy)

dc<-merge(dc.combine,d3.se,by=c("count","type"),all.x=T)
head(dc)
#plot (I would probably get rid of loess stat_smooth bit to make it clearer) and have the y axis in logs
(gy<-ggplot(dy,aes(x=year,y=area,col=type))+
    geom_point()+
    geom_line()+
    xlab("Year")+
    ylab("log10(area),sq-km")+
    geom_errorbar(aes(ymin=area-se.area,ymax=area+se.area),position=position_dodge(width=0.9))+
    scale_y_continuous(trans='log10')+ggtitle("A"))

(gc<-ggplot(dc,aes(x=count,y=area,col=type))+
    geom_point()+
    geom_line()+
    xlab("Number of records")+
    ylab("log10(area),sq-km")+
    geom_errorbar(aes(ymin=area-se.area,ymax=area+se.area),position=position_dodge(width=0.9))+
    scale_y_continuous(trans='log10')+ggtitle("B"))

#now combine plots
g.combine<-grid.arrange(gy,gc,ncol=1)

ggsave("merge_combine.png",g.combine,dpi=300,height=6,width=6)


