
attach.nimble <- function(mcmc.output=mcmc.output){
  library(data.table)
  nsaved = dim(mcmc.output$chain1)[1]
  nchains = length(mcmc.output)
  
  for(i in 1:nchains){
    if(i==1){chain=rep(i,nsaved)}else{chain=c(chain,rep(i,nsaved))}
  }#i
  
  samples <- rbindlist(lapply(mcmc.output, as.data.frame)) %>% add_column(chain=chain)
  
  #Remove brackets indicating matrix structure
  refs <- unique(gsub("\\[[^\\]]*\\]", "", names(samples), perl=TRUE))
  
  for(j in 1:length(refs)){
    
    tempobj <-select(samples,contains(refs[j]))
    
    if(grepl(",",names(tempobj)[1])){
      maxarray <- tail(names(tempobj),1)
      arraydims <- as.numeric(unlist(strsplit(gsub(" ", "",sub(".*\\[([^][]+)].*","\\1",maxarray)),",")))

      output.array <- array(dim=c(dim(samples)[1],arraydims))
      
      for(k in 1:ncol(tempobj)){

        eval(parse(text=paste0("output.array[,",gsub(".*\\[([^.]+)\\].*","\\1",names(tempobj)[k]),"] <- as.data.frame(tempobj)[,",k,"]")))
        
      }#k
      
      assign(refs[j],output.array,envir = parent.frame())
      
    }else{
      assign(refs[j],as.matrix(tempobj),envir = parent.frame())
    }#ifelse
    
  }#j
}#func

