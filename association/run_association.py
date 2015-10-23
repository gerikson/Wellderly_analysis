import os, sys, datetime
"""
Run Association
"""
def create_job_file():


    #infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/chr."+str(sample)
    
    vcffile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_AFmore0.01/final_vcf_allChrom_snps_AF0.01.vcf.gz"
    infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_vcf_allChrom_snps_AF0.01"
    
    
    #command = "/gpfs/home/nwineing/plink --vcf "+vcffile+" --double-id --vcf-half-call m --make-bed --out "+infile + "\n"


    jobfile = jobs_folder  + ".job"         
    outjob = open(jobfile, 'w')
    outjob.write("#!/bin/csh\n")                    
    outjob.write("#PBS -S /bin/bash\n")
    outjob.write("#PBS -l nodes=1:ppn=1\n")
    outjob.write("#PBS -l mem=24gb\n")
    outjob.write("#PBS -l walltime=500:00:00\n")
    outjob.write("#PBS -l cput=9600:00:00\n")
    outjob.write("#PBS -m n\n")

    # name markers according to position #
    '''
    sign = '"_"'
    command += "awk '{print $1, $1"+sign+"$4"+sign+"$5"+sign+"$6, $3, $4, $5, $6}' < " +infile +".bim > "+infile+".temp\n"
    command += "rm "+infile+".bim\n"
    command += "mv "+infile+".temp "+ infile+".bim\n"

    ## ANY OTHER FILTERS (biallelic) ##
    ## ASSSOCIATION TESTING ON WHATEVER FILES CREATED AFTER FILTERS HERE ##

    ##### PCA #####

    # filter by MAF #
    command +="/gpfs/home/nwineing/plink --bfile "+infile+" --maf 0.05 --make-bed --out "+infile+"-0.05_MAF\n"

    # filter by LD #
    command += "/gpfs/home/nwineing/plink --bfile "+infile+"-0.05_MAF --indep-pairwise 50 5 0.5 --out "+infile+"-LD\n"
    command += "/gpfs/home/nwineing/plink --bfile "+infile+"-0.05_MAF --extract "+infile+"-LD.prune.in --make-bed --out "+infile+"-MAF_LD_pruned\n"

    # calculate PCs (first 10) #
    command +="/gpfs/home/nwineing/plink --bfile "+infile+"-MAF_LD_pruned --pca 10 --out "+infile+"-PCA\n"

    
    #Do this later on the combined dataset
    ### plot PCs in R ###
    
    library(ggplot2)
    dat <- read.table("test-PCA.eigenvec", F)
    # first 519 are wellderly, remaining are inova #
    COLOR <- c(rep("red", times=519), rep("blue", times=dim(dat)[1]-519))
    p <- ggplot() +
        geom_point(aes(dat$V3, dat$V4), col=COLOR) +
        xlab("pc1") +
        ylab("pc2")
    tiff("pca.tiff", width=2000, height=2000, res=300, compression="lzw")
    p
    dev.off()
    q()
    


    ##### ASSOCIATION TESTING #####

    # create phenotype file in R #
    # wellderly = "affected" #
    dat <- read.table("test-MAF_LD_pruned.fam", F)
    pheno <- c(rep(2, times=519), rep(1, times=dim(dat)[1]-519))
    out <- data.frame(dat$V1, dat$V2, pheno)
    write.table(out, "test.pheno", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    q()

    '''
    # perform association testing (logistic) using pcs as covariate #
    command = "plink --bfile "+infile+" --logistic --pheno "+infile+".pheno --covar "+infile+"-PCA.eigenvec --allow-no-sex --hide-covar --out "+infile+"-results\n"


    ##### POST-ASSOCIATION PRUNING (for QQ plot) #####
    ##### THIS SCRIPT ASSUMES test-0.05_MAF IS THE FINAL FILTERED DATASET SET #####
    ##### OTHER FILTERS WILL BE USED (and a different dataset) IN THE FINAL ANALYSIS #####
    ##### filter by LD for QQ plot, but not for results #####

    # reduce asssociation testing output #
    command += "awk '{ print $2, $9}' < "+infile+"-results.assoc.logistic > "+infile+"-results.assoc.logistic.short\n"

    # reduce marker names output #
    command += "awk '{ print $2}' < "+ infile+"-0.05_MAF.bim > "+infile+"-0.05_MAF.bim.short\n"         # this is the filtered marker list

    '''
    # prune in R #
    results <- read.table("test-results.assoc.logistic.short", T)
    snp_names <- read.table("test-0.05_MAF.bim.short", F)
    pruned_results <- results[results$SNP %in% snp_names$V1,]
    final_results <- pruned_results[!is.na(pruned_results$P),]

    # create QQ plot (same R session) #
    obs_P <- -log10(final_results$P[order(final_results$P, decreasing=FALSE)])
    exp_P <- -log10(c(1:length(obs_P))/length(obs_P))

    library(ggplot2)
    p <- ggplot() +
        geom_point(aes(exp_P, obs_P)) +
        geom_line(aes(c(0,max(exp_P)), c(0,max(exp_P))), color="red") +
        ylim(c(0,max(exp_P)))

    tiff("test-qq_plot.tiff", width=2000, height=2000, res=300, compression="lzw")
    p
    dev.off()
    
    '''

    outjob.write(command)
    outjob.close()  
    execute = QSUB + ' -e '+ jobs_folder + '.job.err -o ' + jobs_folder   + '.job.out ' + jobfile
    print execute 

    sys.stdout.flush()
    clustnum = os.popen(execute, 'r')
    jobnum = clustnum.readline().strip()
    clustnum.close()


print 'start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q stsi -M gerikson@scripps.edu -l mem=24G -l cput=9600:00:00 -l walltime=500:00:00 "
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/association/run_association."

create_job_file()



