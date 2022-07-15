## HIGH NU ##
#!/usr/bin/env Rscript
#options(repos=structure(c(CRAN="http://cran.mirror.garr.it/mirrors/CRAN/")))
#library(meanShiftR)

### Adapted from Amicone & Gordo 2021 to output mutation trajectories over time.
### This algorithm models adaptation across different mutation rates (Us) and keeps the other parameters constant.
## It will output a folder containing summarized data files for each condition
# and it does NOT accumulate memory during the process.

args<-commandArgs(trailingOnly = TRUE) ## Enables to read parameters from command line

##########################################  Set the fixed parameters ##########################################
##### Seeding parameter ######
### change if you want to run further simulations with the same set of parameters
seed_x=0

######## Algorithm Parameters #######
t_tot=10000    ## Total time/generations
d=1           ## death rate
Pop=10^7      ## Population size

#### Evolution Parameters ####
Pb=0.5          ## Probability of beneficial mutations. (1-Pb) is deleterious
pleiotropy=0    ## Pleiotropic effect
n0=1            ## initial density of the mutants
E_constraint=1  ## Energy constraint (has to be < Inf)
#########

### Output parameters ###
n_sample=100 ## sample size 
p_sample=100 ## sampling period
t_sample=c(1,seq(p_sample,t_tot,p_sample)) ## sample at given generations

##########################################  Set and check the custom parameters ##########################################
### parameters are given by the user through command line as:
##  Rscript **.R arg1 arg2 arg3 arg4 arg5 \n
##  The 1st argument (arg1) is the number of resources/traits (must be integer and >0).
##  The 2nd argument (arg2) is the beneficial mutational load NU_b (must be >0) and assumes Pb=0.5, N=10^7.
##  The 3rd argument (arg3) is the expected mutation step (must be >0).
##  The 4th argument (arg4) is the initial fitness of the population (must be 0<arg4<=E_constraint).
##  The 5th argument (arg5) is the number of replicas/populations (must be integer and >1).

## double-check the correct number of parameters
if(length(args)!=5){stop("Incorrect number of parameters!
       Try again as follows:
       Rscript **.R arg1 arg2 arg3 arg4 arg5 \n
       The 1st argument (arg1) is the number of resources/traits (must be integer and >0).
       The 2nd argument (arg2) is the beneficial mutational load NU_b (must be >0).
       The 3rd argument (arg3) is the expected mutation step (must be >0).
       The 4th argument (arg4) is the initial fitness of the population (must be 0<arg4<=1).
       The 5th argument (arg5) is the number of replicas/populations (must be integer and >1).")}

R=as.numeric(args[1])
NU_b=as.numeric(args[2])        ## Beneficial mutational load (get U by assuming Pb=0.5, N=10^7)
Delta=as.numeric(args[3])       ## Mutation step (fixed or standard deviation to be set at line X)
inFit=as.numeric(args[4])       ## Sum of the alphas (initial fitness)
simulations=as.numeric(args[5]) ## Number of simulations

## double-check the correct format of the parameters
if(R>1 && R%%1==0 && NU_b>0 && Delta>=0 && inFit>0 && inFit<=E_constraint && simulations>1 && simulations%%1==0){
  print("Input given by the user:")
  print(paste("Number of resources/traits = ",R,sep=""))
  print(paste("Beneficial mutation load = ",NU_b,sep=""))
  print(paste("Mean beneficial mutation step = ",Delta,sep=""))
  print(paste("Initial fitness of the population = ",inFit,sep=""))
  print(paste("Number of replicas = ",simulations,sep=""))
}else{stop("Input not valid!
       Try again as follows:
       Rscript **.R arg1 arg2 arg3 arg4 arg5 \n
       The 1st argument (arg1) is the number of resources/traits (must be integer and >0).
       The 2nd argument (arg2) is the beneficial mutational load NU_b (must be >0).
       The 3rd argument (arg3) is the expected mutation step (must be >0).
       The 4th argument (arg4) is the initial fitness of the population (must be 0<arg4<=1).
       The 5th argument (arg5) is the number of replicas/populations (must be integer and >1). \n
       Note: Math expressions (e.g. 10^-3) are NOT allowed!
       Use full format instead (e.g. 0.001).")}

## Ancestor phenotypes
ancestor=rep(inFit/R,R) ## here assigned equal preference

## Mutations
Sigma_s=Delta*sqrt(pi)/sqrt(2) ## so that the expected beneficial delta=Delta
Us=NU_b/(Pop*Pb)

############ resource supply ##########
sS=unlist(lapply(1:R, function(x)(Pop/R))) ##Here, resources are equally distributed

##### Output file name ###
model=paste("EcoEvo_",R,"Res",sep="") ##suffix of the output file

OutFile_suffix=paste("_N",Pop,"_Delta",Delta,"_Pb",Pb,"_W0_",inFit,"_",t_tot,"gen_",simulations,"sim_seed",seed_x,"_",model,sep="")
#########

print("###########################")
print("Fixed parameters:")
print(paste("Population size = ",Pop,sep=""))
print(paste("Proportion of ben. mutations = ",Pb,sep=""))
print(paste("Number of generations = ",t_tot,sep=""))
print(paste("Energetic constraint = ",E_constraint,sep=""))
print(paste("Sampling size = ",n_sample,sep=""))
print(paste("Simulation seed = ",seed_x,sep=""))
print("###########################")
if(NU_b>5){
  print("Warning:")
  print("Simulations with mutational load >5 are slow and might take very long time.")
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
Gamma_DME <- function(D,shape,rate) {
  Deltas=rep(0,D)        ##initialize vector
  gene=sample(1:D,1)     ##sample 1 random trait
  Deltas[gene]=rgamma(1, shape=shape, rate = rate)
  return(Deltas) 
}

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


print(paste("#simulations:",simulations,", Population size:",Pop,", Mutation rate (U):",Us,", Mutation step:",Delta,", Pleiotropy: ",pleiotropy))
print(paste("Number of resources:",R))


############ Start the adaptation simulations ###############

##Initialize Data structures
counter=0 ## parameter iteration

### Computational time table
CPtime_table=data.frame(matrix(rep(0,length(Us)),nrow =1,ncol = length(Us)))
rownames(CPtime_table)=pleiotropy
colnames(CPtime_table)=Us
  
print(paste("##########################################",R," Resources ###########################################"))
  
### iterations over different mutation rates
for (u in 1:length(Us)){
    U=Us[u]
    pt0=proc.time()[3] ## reset the CP time
    counter=counter+1 ## update the count of parameter combinations
      
    print(paste("With U:",Us[u],", mutation step ",Delta," and pleiotropy:",pleiotropy,""))
    
    ######### Initialize matrices for summary statistics ##########

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
      
      #### alphas (phenotypes) ####
      ## build an empty data frame for the phenotypes
      alphas=data.frame(matrix(rep(0,N*R),nrow=N,ncol=R))
      ## assign th ancestor phenotype
      alphas[1,]=ancestor
        
      history=list()
      history[["1"]]=c(1) ## history stores the phylogenetic info
        
      Ns=c(1) ## count all the genotypes
      trackID=c(1) ## to keep track of the genotypes' identity (haplotype)
      mutID=c(1) ## to keep track of the mutations' identity (single mutation)
      surv_pos=c(1)
        
      State=c(Pop,rep(0,N-1)) ## Initial density
      enzymes=unlist(lapply(1:R,function(x)(ancestor[x]*State[1])))
      
      ######### Initialize vectors for summary statistics ##########
      fit_Mean=inFit
      Mut_count=matrix(1,nrow = 1,ncol = 1+length(t_sample))
      colnames(Mut_count)=c("mut_ID",t_sample)
      mut_parent=c(1)
      mut_pheno=t(matrix(ancestor))
      colnames(mut_pheno)=colnames(alphas)
      
      ########### Iterate over time/generations ############
      ######################################################
      k=1
      for (generation in 2:t_tot){
        t_alphas=alphas[surv_pos,] ## get the alphas only of the alive genotypes
        Nt=Ns 
        parents=vector()
        
        ## mutation iteration for each genotype
        for(x in 1:N){
          parental=t_alphas[x,] ## parental phenotype
          ###### Poisson process ######
          exp_mut=rpois(1,State[x]*U) ## number of mutations on genotype x
          
          ## iteration of mutations on genotype x  
          if(exp_mut>0){
            parents=c(parents,rep(mutID[surv_pos[x]],exp_mut))
            for (trial in 1:exp_mut){
              mutant=rep(Inf,R)
              ## draw the mutation effect within the allowed region (sum(alphas)<=E)
              while(sum(mutant)>E_constraint){
                
                deltas=Normal_DME(R,Sigma_s,pleiotropy) ## Normal distribution
                #deltas=Fixed_s(R,Delta,pleiotropy)    ## Fixed effect of size +/- Delta
                mutant=pmax(rep(0,R),parental+deltas)
              }
              # update genotypes counts
              Nt=Nt+1
              trackID=c(trackID,Nt)
              State=c(State,n0)
              t_alphas=rbind(t_alphas,mutant)
            }
            
             State[x]=max(0,State[x]-exp_mut) ## update densities
          }
            
        }
        tN=N+(Nt-Ns) ## temporary number of strains
        enzymes=unlist(lapply(1:R,function(x)(sum(State*t_alphas[,x])))) ## compute the state of the ecosystem
        growths=unlist(lapply(1:tN, function(s) sum(unlist(lapply(1:R,function(x) t_alphas[s,x]*sS[x]/enzymes[x]))))) # compute the expected growth
        exTout=pmax(rep(0,tN),State*(1+growths-d)) #expected update of the density
        
        ## multinomial sampling from the expected densities
        Tout=rmultinom(1,sum(exTout),exTout/sum(exTout))  ## Selection + Drift
        #Tout=rmultinom(1,sum(State),State/sum(State))    ## Pure drift
        #Tout=exTout                                      ## Pure selection
        
        densities=as.numeric(Tout)
        surv_id=which(densities>=n0) ##eliminate the extinct
        New=setdiff(surv_id,1:N) ## get the new genotypes
        old=setdiff(surv_id,New)
          
        ## update the phylogenetic info and store only the mutations that are still present
        trackID=trackID[old]
        
        ## get all the mutations within the surv genotypes
        type=ID2string(trackID)
        surv_mut=unique(unlist(history[type]))
        #type_mut=ID2string(surv_mut)
        type_mut=unique(ID2string(c(surv_mut,parents))) #save also the parent history because it get get lost by extinction 
        history=history[type_mut] ##update mutation
       
        old_mut=mutID%in%surv_mut
        mutID=mutID[old_mut]
        alphas=alphas[old_mut,]
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
        
        ## update alphas
        alphas=rbind(alphas,t_alphas[New,])
        rownames(alphas)=mutID
        
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
          s_alphas=alphas[samples,] ##sampled alphas
          
          ##phenotypic
          fit=rowSums(s_alphas)
          fit_Mean=c(fit_Mean,mean(fit))
          
          ## mutation trajectories
          k=k+1
          mutations=unlist(history[types])
          mutations=mutations[mutations!=1]
          
          for (j in mutations){
              if(!any(Mut_count[,1]==j)){
                Mut_count=rbind(Mut_count,c(j,rep(0,length(t_sample))))
                mut_parent=c(mut_parent,history[[ID2string(j)]][2])
                mut_pos=which(mutID==j)
                mut_pheno=rbind(mut_pheno,alphas[mut_pos,])
              }
              x=which(Mut_count[,1]==j)
              Mut_count[x,k+1]=Mut_count[x,k+1]+1/n_sample ##add frequency
          }
          
        } ######### end of sampling #########
        
      } ### end of 1 single simulation #######
      
      ## update output matrices
      Fit_Mean[replica,]=fit_Mean 
      
      x=which(Mut_count[,1]!=1)
      mut_parent=mut_parent[x]            ##remove ancestor
      mut_pheno=as.matrix(mut_pheno[x,])  ##remove ancestor but keep matrix form
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
      Mut_sum=c(R,Pop,U,Delta,inFit,Pb,replica,M_t) ##concatenate parameters, population ID and M_t
      Mut_Accumulation=rbind(get0("Mut_Accumulation"),Mut_sum)
      
      ##add mutations details
      Mut_d=cbind(pop_ID,mut_ID,mut_parent,mut_pheno)
      Mut_Details=rbind(get0("Mut_Details"),Mut_d)
      
      ######## Output summary data for the simulated condition ##########
      pop_ID=c(1:replica,rep(NA,(simulations-replica)))
      
      data=Fit_Mean
      data=cbind(pop_ID,data)
      #write.csv(data,paste(folder,"Fit_Mean_U",Us[u],OutFile_suffix,".csv",sep=""),row.names = F)
      
      data=Mut_Trajectory
      write.csv(data,paste(folder,"MutationTrajectories_U",Us[u],OutFile_suffix,".csv",sep=""),row.names = F)
      
      data=Mut_Accumulation
      colnames(data)=c("Resources","N","U","delta","w0","Pb","replica",t_sample)
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

