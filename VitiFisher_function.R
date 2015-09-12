# pathway enrichemnts from 
# a. a two-column list of genes and their pathways
# b. a single column list of all possible gene names
# c. a file, or files with lists of genes in a group, one per line

## all three files must match up (genenames)

# will output a file for every input file with categories present in the gene list
# including fisher p-value, number of genes in category, number of genes in genome
# for that category, and adjusted p-value. 

### now lets discuss the input in more detail

# this next file is file a. It never changes unless pathway annotations are updated
genes.annotations <- read.table("vitisnet_funcCats_clean.txt",sep="\t",colClasses=rep("character",2),header=F)
colnames(genes.annotations) <- c("genename","category")
# this next file is file b. It never changes unless gene names change. If they do, then
# you must also update file a
vitis.genenames <- as.character(read.table("vitis.genenames.txt")[,2])
# this is the variable character in the list of files, could just be one file
# if you don't have more than one gene list
groups <- c("C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11",
            "C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23")
# change the "paste" command inside the loop to match the filenames, 
# change the "filename" structure at the end of the loop to match what you want your output to look like
# and you're set... run this bad boy
####################

library('plyr')
for(i in 1:length(groups)){
    which.genes <- read.table(paste("clusters",groups[i],"genelist.txt",sep="."),header=F,colClasses="character")[,1]
    which.genes <- as.factor(as.numeric(vitis.genenames %in% which.genes))
    names(which.genes) <- vitis.genenames
    functional.category.counts <- rep(0,length(table(genes.annotations$category)))
    names(functional.category.counts) <- names(table(genes.annotations$category))
    locations <- which(genes.annotations$genename %in% 
                           names(which.genes[which.genes == 1]))
    functions <- table(genes.annotations$category[locations]) #get functional category counts
    print(length(functions))
    functional.category.counts[which(names(functional.category.counts) %in% 
                                         names(functions))] <- functions #save to table
    print(length(which(names(functional.category.counts) %in% 
                           names(functions))))
    
    #now to do the enrichment
    genome.totals <- table(genes.annotations$category)
    sigtallysall <- c()
    total.gene.count <- length(functional.category.counts[functional.category.counts > 0])
    testlocs <- functional.category.counts[which(functional.category.counts > 0)];
    if(length(testlocs) > 0){
        for(l in 1:length(testlocs)){
            thistest <- names(testlocs[l])
            in.group <- testlocs[l] #count subset in
            out.group <- total.gene.count - in.group #count subset out
            in.genome <- genome.totals[which(names(genome.totals) %in% 
                                                 thistest)] - in.group #count genome in
            out.genome <- 29971 - total.gene.count - in.genome #count genome out  
            fisherset <- matrix(c(in.group,in.genome,out.group,out.genome),
                                nrow=2,ncol=2)
            testRes <- fisher.test(fisherset,alternative="greater")
            if(in.group > 0 & length(sigtallysall) < 1){
                sigtallysall <- c(names(in.genome),testRes$p.value,in.group,
                                  genome.totals[which(names(genome.totals) %in% thistest)])
            }else if(in.group > 0){
                sigtallysall <- rbind(sigtallysall,c(names(in.genome),testRes$p.value,
                                                     in.group,genome.totals[which(names(genome.totals) %in% thistest)]))
            }
        }
    }
    rownames(sigtallysall) <-c(1:length(sigtallysall[,1]))
    colnames(sigtallysall) <-c("category","p.value",
                               "count.in.category","count.in.comparison")
    filename = paste("clusters",groups[i],"genelist.enrichments.txt",sep=".")
    sigtallysall <- cbind(sigtallysall,"adjusted.p"=p.adjust(sigtallysall[,2]))
    write.table(file=filename,sigtallysall)
}