require("doParallel")
require("foreach")

# Note this example was intended for a situation with many different time
# series that each need stl interpolation performed.  For example:
# many different district's agricultural time series, each with NA values 
# to be estimated. Contact me if you need any help getting this to work.


####Run STLinterp for Parallelization####
#first break up into (n_cores) even sets the districts with missing values
#to improve efficiency of parallel CPUs. 

breaker<-floor(length(missing_names)/3)
b1 <- breaker
b2 <- breaker*2+1
part1 <- missing_names[1:b1]
part2 <- missing_names[(b1+1):(b2)]
part3 <- missing_names[(b2+1):length(missing_names)]

joblist <- list(part1, part2, part3)

partition <- vector(mode="list",length=3)
for (i in 1:3){
  partition[[i]]<-master %>% filter(state.dist %in% joblist[[i]])
}

#init results storage (requires grid to have been created)
colnam<-  c("s.window", "s.degree", "t.window", "t.degree",
            "fc.window", "fc.degree","CVscore")
r1 <- data.frame(matrix(nrow=length(part1), ncol=8,
                  dimnames=list(c(),c("state.dist",colnam))),
           stringsAsFactors=F)
r2 <- data.frame(matrix(nrow=length(part2), ncol=8,
                           dimnames=list(c(),c("state.dist",colnam))),
                    stringsAsFactors=F)
r3 <- data.frame(matrix(nrow=length(part3), ncol=8,
                           dimnames=list(c(),c("state.dist",colnam))),
                    stringsAsFactors=F)

results <- list(r1,r2,r3) 

#### Run over all in parallel####
n_c <- detectCores()-1
cl<- makeCluster(n_c,type = "PSOCK")
registerDoParallel(cl)

# dopar over the joblsit (takes about 15 hours)
st2 <- system.time(
  res <- foreach(i=1:3,
                      .packages =c("tidyverse","stlplus"))
   %dopar% {
     for (j in 1:length(joblist[[i]])){
      ####need to do more to make temp a correct time series object
      temp <- partition[[i]] %>% filter(state.dist == joblist[[i]][j]) %>%
        arrange(time)
      temp.ts <- ts(temp$med.welldepth,
         start= c(temp$year[1],time_func(temp$season[1])),
         frequency = 4)

      #store name, store results
      results[[i]][j,1]<- joblist[[i]][j]
      results[[i]][j,2:8] <- STLinterp(temp.ts, s.window=c(4,10,20,35,55,70), t.window = c(NA,4,10,20,35,55,70),
                                       t.degree = c(1,2), s.degree=c(1,2), fc.window = c(NA, 4, 10, 50), fc.degree = NA,
                                       n =40, type="kfold")
     }
    results[[i]]
  }
)
stopCluster(cl)
