## HIGH NU ##
#!/usr/bin/env Rscript
#options(repos=structure(c(CRAN="http://cran.mirror.garr.it/mirrors/CRAN/")))
#library(meanShiftR)

### Model selection effects from a Gamma distributions, instead of fixed s.
############ Model WITHOUT ECOLOGICAL SELECTION ###################################

### This algorithm models adaptation across different mutation rates (Us) and keeps the other parameters constant.
## It will output a folder containing summarized data files for each condition
# and it does NOT accumulate memory during the process.

args<-commandArgs(trailingOnly = TRUE) ## Enables to read parameters from command line

##########################################  Set the fixed parameters ##########################################
##### Seeding parameter ######
### change if you want two run further simulations with the same set set of parameters
seed_x=0

######## Algorithm Parameters #######
t_tot=10000    ## Total time/generations
Pop=10^7       ## Population size

#### Evolution Parameters ####
Pb=1            ## Probability of beneficial mutations. (1-Pb) is deleterious
n0=1            ## initial density of the mutants
#########

### Output parameters ###
n_sample=100 ## sample size 
p_sample=100 ## sampling period
t_sample=c(1,seq(p_sample,t_tot,p_sample)) ## sample at given generations


###
model="DirSel_Gamma" ##suffix of the output file

##########################################  Set and check the custom parameters ##########################################
### parameters are given by the user through command line as:
##  Rscript **.R arg1 arg2 arg3 arg4 \n
##  The 1st argument (arg1) is the beneficial mutational load NU_b (must be >0) and assumes Pb and N as defined above.
##  The 2nd argument (arg2) is the expected mutation step (must be >0).
##  The 3th argument (arg3) is the initial fitness of the population (must be >=0).
##  The 4th argument (arg4) is the number of replicas/populations (must be integer and >1).
##  The 5th argument (arg5) is the shape parameter of the Gamma distribution (must be >0).

## double-check the correct number of parameters
if(length(args)!=5){stop("Incorrect number of parameters!
       Try again as follows:
       Rscript **.R arg1 arg2 arg3 arg4 \n
       The 1st argument (arg1) is the beneficial mutational load NU_b (must be >0).
       The 2nd argument (arg2) is the expected mutation step (must be >0).
       The 3th argument (arg3) is the initial fitness of the population (must be >=0).
       The 4th argument (arg4) is the number of replicas/populations (must be integer and >1).
       The 5th argument (arg5) is the shape parameter of the Gamma distribution (must be >0).")}

NU_b=as.numeric(args[1])        ## Beneficial mutational load (get U by assuming Pb and N as defined above)
Delta=as.numeric(args[2])       ## Mutation step (fixed or standard deviation to be set at line X)
inFit=as.numeric(args[3])       ## Initial fitness
simulations=as.numeric(args[4]) ## Number of simulations
shape_G=as.numeric(args[5])     ## Shape parameter of the Gamma distribution

## double-check the correct format of the parameters
if(NU_b>0 && Delta>=0 && inFit>=0 && simulations>1 && simulations%%1==0 && shape_G>0){
  print("Input given by the user:")
  print(paste("Beneficial mutation load = ",NU_b,sep=""))
  print(paste("Mean beneficial mutation step = ",Delta,sep=""))
  print(paste("Initial fitness of the population = ",inFit,sep=""))
  print(paste("Number of replicas = ",simulations,sep=""))
  print(paste("Shape of Gamma = ",shape_G,sep=""))
}else{stop("Input not valid!
       Try again as follows:
       Rscript **.R arg1 arg2 arg3 arg4 \n
       The 1st argument (arg1) is the beneficial mutational load NU_b (must be >0).
       The 2nd argument (arg2) is the expected mutation step (must be >0).
       The 3th argument (arg3) is the initial fitness of the population (must be >=0).
       The 4th argument (arg4) is the number of replicas/populations (must be integer and >1).
       The 5th argument (arg5) is the shape parameter of the Gamma distribution (must be >0).\n
       Note: Math expressions (e.g. 10^-3) are NOT allowed!
       Use full format instead (e.g. 0.001).")}

## Ancestor fitness
ancestor=inFit 

## Mutation rate
rate_G=shape_G/Delta ## so that the mean mutation effect is Delta
Us=NU_b/(Pop*Pb)

##### Output file name ###
OutFile_suffix=paste("_N",Pop,"_Delta",Delta,"_Pb",Pb,"_W0_",inFit,"_",t_tot,"gen_",simulations,"sim_seed",seed_x,"_",model,"_shape",shape_G,sep="")
#########

print("###########################")
print("Fixed parameters:")
print(paste("Population size = ",Pop,sep=""))
print(paste("Proportion of ben. mutations = ",Pb,sep=""))
print(paste("Number of generations = ",t_tot,sep=""))
print(paste("Sampling size = ",n_sample,sep=""))
print(paste("Simulation seed = ",seed_x,sep=""))
print("###########################")
if(NU_b>5){
  print("Warning:")
  print("Simulations with mutational load >5 are slow and might take long time.")
}

##########################################  End of parameters ##########################################


##########################################  Built-in functions ##########################################
## convert a number variable into string
ID2string <- function(x) {
     return(format(x,scientific = FALSE,trim=TRUE)) 
}

## Mutation function
## takes as input the #of resources or dimensions (D), the standard deviation (sd) and the pleiotropic effect (p)
## and outputs a vector of dimension D of mutations drawn from a normal distribution of mutation effects (DME)
Normal_DME <- function(D,sd,p) {
  Deltas=rnorm(D,mean=0,sd=sd*p)        ##if pleiotropy=0 this is equivalent to Deltas=rep(0,D)
  gene=sample(1:D,1)
  Deltas[gene]=rnorm(1,mean=0,sd=sd)
  return(Deltas) 
}

## takes as input the #of resources or dimensions (D), the mutation effect (s) and the pleiotropic effect (p)
## and outputs a vector of dimension D where D-1 traits change by (+/- s)*p and one trait changes by +/- s.
Fixed_s <- function(D,s,p) {
  Signs=sample(c(1,-1),D,replace = T,prob = c(Pb,1-Pb))       ## probability of beneficial (Pb) and deleterious (1-Pb)
  Deltas=Signs*s*p                                            ## if pleiotropy=0 this is equivalent to Deltas=rep(0,D)
  gene=sample(1:D,1)
  Deltas[gene]=Signs[gene]*s
  return(Deltas) 
}

##########################################  End of built-in functions ##########################################

################### 
print("Model of directional selection (No Epistasis)")
print(paste("#simulations:",simulations,", Population size:",Pop,", Mutation rate (U):",Us,", Mean mutation step:",shape_G/rate_G))

############ Start the adaptation simulations ###############

##Initialize Data structures
counter=0 ## parameter iteration

### Computational time table
CPtime_table=data.frame(matrix(rep(0,length(Us)),nrow =1,ncol = length(Us)))
colnames(CPtime_table)=Us
  
### iterations over different mutation rates
for (u in 1:length(Us)){
    U=Us[u]
    pt0=proc.time()[3] ## reset the CP time
    counter=counter+1 ## update the count of parameter combinations
      
    ######### Initialize matrices for summary statistics ##########
    Fit_Var=matrix(nrow=simulations,ncol=length(t_sample))
    colnames(Fit_Var)=t_sample
    Fit_Mean=matrix(nrow=simulations,ncol=length(t_sample))
    colnames(Fit_Mean)=t_sample
    
    ######## Create Folder where to save the output ###########
    folder=paste("Data_U",Us[u],OutFile_suffix,"/",sep="")
    dir.create(folder,showWarnings = FALSE)
    
    ### iterate independent populations
    for(replica in 1:simulations){
      set.seed(replica+simulations*seed_x) ##set the seeds
      print(paste("simulation",replica, " has started!"))
      
      #### Initial conditions 
      N=1  ### number of initial species
      Ndym=c(N)  ## strain diversity over time
      
      #### Fitness ####
      ## build an empty data frame for the fitness
      fitness=ancestor
        
      history=list()
      history[["1"]]=c(1) ## history stores the phylogenetic info
        
      Ns=c(1) ## count all the genotypes
      trackID=c(1) ## to keep track of the genotypes' identity (haplotype)
      mutID=c(1) ## to keep track of the mutations' identity (single mutation)
      surv_pos=c(1)
        
      State=c(Pop,rep(0,N-1)) ## Initial density
      
      ######### Initialize vectors for summary statistics ##########
      fit_Var=0
      fit_Mean=sum(State*fitness)/sum(State)
      Mut_count=matrix(1,nrow = 1,ncol = 1+length(t_sample))
      colnames(Mut_count)=c("mut_ID",t_sample)
      mut_parent=c(1)
      mut_fit=fitness
      
      ########### Iterate over time/generations ############
      ######################################################
      k=1
      for (generation in 2:t_tot){
        t_fitness=fitness[surv_pos] ##get the temporary fitness vector
        Nt=Ns 
        parents=vector()
        
        ## mutation iteration for each genotype
        for(x in 1:N){
          parental=t_fitness[x] ## parental fitness
          ###### Poisson process ######
          exp_mut=rpois(1,State[x]*U) ## number of mutations on genotype x
          
          ## iteration of mutations on genotype x  
          if(exp_mut>0){
            parents=c(parents,rep(mutID[surv_pos[x]],exp_mut))
            for (trial in 1:exp_mut){
              ## Mutation effect
              mut_effect=rgamma(1, shape=shape_G, rate = rate_G) ## Gamma distributed ben. mutations
              mutant=parental+mut_effect
              # update genotypes counts
              Nt=Nt+1
              trackID=c(trackID,Nt)
              State=c(State,n0)
              t_fitness=c(t_fitness,mutant)
            }
            
             State[x]=max(0,State[x]-exp_mut) ## update densities
          }
            
        }
        tN=N+(Nt-Ns) ## temporary number of strains
        
        ########  Directional Selection: The growth depends only on fitness 
        sum_Fit=sum(State*exp(t_fitness))
        X_m=log(sum_Fit/sum(State))
        Probs=State*(1+t_fitness-X_m)/sum(State) ## n(t+1)=n(t)(1+x-x_m)
        
        ## multinomial sampling from the expected densities
        Tout=rmultinom(1,sum(State),prob = Probs)          ## DirSelection + Drift
        
        densities=as.numeric(Tout)
        surv_id=which(densities>=n0) ## eliminate the extinct
        New=setdiff(surv_id,1:N)     ## get the new genotypes
        old=setdiff(surv_id,New)
          
        ## update the phylogenetic info and store only the mutations that are still present
        trackID=trackID[old]
        
        ## get all the mutations within the surv genotypes
        type=ID2string(trackID)
        surv_mut=unique(unlist(history[type]))
        type_mut=unique(ID2string(c(surv_mut,parents))) #save also the parent history because it get get lost by extinction 
        history=history[type_mut]    ##update mutation
       
        old_mut=mutID%in%surv_mut
        mutID=mutID[old_mut]
        fitness=fitness[old_mut]
        surv_pos=which(mutID %in% trackID)
        
        ### add new mutation 
        if(length(New)>0){
          for(new in 1:length(New)){
            mutID=c(mutID,Ns+new)
            trackID=c(trackID,Ns+new)
            surv_pos=c(surv_pos,length(mutID))
            parent=ID2string(parents[New[new]-N]) 
            history[[ID2string(Ns+new)]]=c(Ns+new,history[[parent]])
          }
        }
        
        ## update fitness
        fitness=c(fitness,t_fitness[New])
        
        ## update genotypes count and densities
        N=length(surv_id)
        Ndym=c(Ndym,N)
        Ns=Ns+length(New)
        State=densities[surv_id]
       
        ################################ sample and output stats ##########################
        
        if(generation %in% t_sample){
          densities=State
          surv=surv_pos
          #set.seed(replica)
          samples=sample(rep(surv,round(densities)),n_sample) ##sampled population
          types=ID2string(mutID[samples])
          fit=exp(fitness[samples]) ##sampled fitnesses
          
          # Mean and Variance fitness
          fit_Mean=c(fit_Mean,log(mean(fit)))
          fit_Var=c(fit_Var,var(log(fit)))
          
          ## mutation trajectories
          k=k+1
          mutations=unlist(history[types])
          mutations=mutations[mutations!=1]
          
          for (j in mutations){
              if(!any(Mut_count[,1]==j)){
                Mut_count=rbind(Mut_count,c(j,rep(0,length(t_sample))))
                mut_parent=c(mut_parent,history[[ID2string(j)]][2])
                mut_pos=which(mutID==j)
                mut_fit=c(mut_fit,fitness[mut_pos])
              }
              x=which(Mut_count[,1]==j)
              Mut_count[x,k+1]=Mut_count[x,k+1]+1/n_sample ##add frequency
          }
          
        } ######### end of sampling #########
        
      } ### end of 1 single simulation #######
      
      ## update output matrices
      Fit_Mean[replica,]=fit_Mean
      Fit_Var[replica,]=fit_Var
      
      x=which(Mut_count[,1]!=1)
      mut_parent=mut_parent[x]            ##remove ancestor
      mut_fit=mut_fit[x]                  ##remove ancestor but keep matrix form
      pop_ID=rep(replica,length(x))
      
      #update mutation counts, remove ancestor id and keep matrix form
      if(length(x)==1){
        Mut_count=t(as.matrix(Mut_count[x,])) 
        Mut=c(pop_ID,Mut_count)
      }else{Mut_count=as.matrix(Mut_count[x,])
      Mut=cbind(pop_ID,Mut_count)}
      
      mut_ID=Mut_count[,1] 
      Mut_Trajectory=rbind(get0("Mut_Trajectory"),Mut)
      
      ##compute the total derived allele frequency (DAF)
      M_t=colSums(Mut_count)[2:dim(Mut_count)[2]]   ##Sums the frequencies of all mutations at each time point
      Mut_sum=c(Pop,U,Delta,inFit,Pb,replica,M_t)   ##concatenate parameters, population ID and M_t
      Mut_Accumulation=rbind(get0("Mut_Accumulation"),Mut_sum)
      
      ##add mutations details
      Mut_d=cbind(pop_ID,mut_ID,mut_parent,mut_fit)
      Mut_Details=rbind(get0("Mut_Details"),Mut_d)
      
      ######## Output summary data for the simulated condition ##########
      pop_ID=c(1:replica,rep(NA,(simulations-replica)))
      
      data=Fit_Var
      data=cbind(pop_ID,data)
      #write.csv(data,paste(folder,"Fit_Var_U",Us[u],OutFile_suffix,".csv",sep=""),row.names = F)
      
      data=Fit_Mean
      data=cbind(pop_ID,data)
      #write.csv(data,paste(folder,"Fit_Mean_U",Us[u],OutFile_suffix,".csv",sep=""),row.names = F)
      
      data=Mut_Trajectory
      write.csv(data,paste(folder,"MutationTrajectories_U",Us[u],OutFile_suffix,".csv",sep=""),row.names = F)
      
      data=Mut_Accumulation
      colnames(data)=c("N","U","delta","w0","Pb","replica",t_sample)
      #write.csv(data,paste(folder,"Mut_Accumulation_U",Us[u],OutFile_suffix,".csv",sep=""),row.names = F)
      
      data=Mut_Details
      #write.csv(data,paste(folder,"Mut_Details_U",Us[u],OutFile_suffix,".csv",sep=""),row.names = F)
      
      ################## end of output #############################
      
    } ### end of all simulations for a given u ###
  
  CPtime_table[1,u]=(proc.time()[3]-pt0)/simulations      
  print(paste("It took:",round(CPtime_table[1,u],digits=3), " seconds" , "(= ",round(CPtime_table[1,u]/3600,digits=3)," hours)  per simulation."))
  print("#####################################################################")
  
  rm(list=c("Mut_Trajectory","Mut_Details","Mut_Accumulation"))
  
  ## save every iteration
  #save.image(file=OutFile)
}
 
#save.image(file=OutFile)

