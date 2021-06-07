######## Diversity Calcs ##########
library(vegan)
library(phyloseq)

leaf.div<-estimate_richness(leaf, measures=c('InvSimpson', "ACE"))

root.div<-estimate_richness(root.4, measures=c('InvSimpson', "ACE"))

soil.div<-estimate_richness(soil.4, measures=c('InvSimpson', "ACE"))

write.csv(leaf.div, '~/Documents/canola/leaf.div.csv')

leaf.div.1<-read.csv('~/Documents/canola/leaf.div.csv')

root.sam<-data.frame(sample_data(root.4))

root.div.1<-na.omit(merge(pi.root, root.div, by='row.names'), row.names=1)
rownames(root.div.1)<-root.div.1[,1]
root.div.1<-root.div.1[,-1]
root.div.2<-merge(root.sam, root.div.1, by='row.names')
write.csv(root.div.2, '~/Documents/root.div.csv')         

soil.div.1<-na.omit(merge(pi.soil, soil.div, by='row.names'), row.names=1)
rownames(soil.div.1)<-soil.div.1[,1]
soil.div.1<-soil.div.1[,-1]
soil.div.2<-merge(soil.sam, soil.div.1, by='row.names')
write.csv(soil.div.2, '~/Documents/soil.div.csv')         

leaf.div.1<-na.omit(merge(pi, leaf.div, by='row.names'), row.names=1)
rownames(leaf.div.1)<-leaf.div.1[,1]
leaf.div.1<-leaf.div.1[,-1]
leaf.div.2<-merge(leaf.sam, leaf.div.1, by='row.names')
write.csv(leaf.div.1, '~/Documents/leaf.div.csv')         



names(root.div.1)

r.simp.mean<- ddply(root.div.1, "Week", summarise,
                    N    = length(InvSimpson),
                    mean = mean(InvSimpson),
                    sd   = sd(InvSimpson),
                    se   = sd / sqrt(N))
write.csv(r.simp.mean, '~/Documents/r.simp.mean.csv')

r.ACE.mean<- ddply(root.div.1, "Week", summarise,
                   N    = length(ACE),
                   mean = mean(ACE),
                   sd   = sd(ACE),
                   se   = sd / sqrt(N))
write.csv(r.ACE.mean, '~/Documents/r.ACE.mean.csv')

s.simp.mean<- ddply(soil.div.1, "Week", summarise,
                    N    = length(InvSimpson),
                    mean = mean(InvSimpson),
                    sd   = sd(InvSimpson),
                    se   = sd / sqrt(N))
write.csv(s.simp.mean, '~/Documents/s.simp.mean.csv')

s.ACE.mean<- ddply(soil.div.1, "Week", summarise,
                   N    = length(ACE),
                   mean = mean(ACE),
                   sd   = sd(ACE),
                   se   = sd / sqrt(N))


write.csv(s.ACE.mean, '~/Documents/s.ACE.mean.csv')


#Pielou's evenness

H.leaf<- diversity(leaf.abun)
pi.leaf <- as.data.frame(H.leaf/log(specnumber(leaf.abun)))
names(pi.leaf)[1] <- "pi"

soi.abun<-data.frame(t(otu_table(soil.4)))

H.soil<- diversity(soi.abun)
pi.soil <- as.data.frame(H.soil/log(specnumber(soi.abun)))
names(pi.soil)[1] <- "pi"

pi.soil.1<-na.omit(merge(soil.sam, pi.soil, by='row.names'))
names(pi.soil.1)

s.pi.mean<-ddply(pi.soil.1, "Week", summarise,
                 N    = length(pi),
                 mean = mean(pi),
                 sd   = sd(pi),
                 se   = sd / sqrt(N))


write.csv(s.pi.mean, '~/Documents/s.pi.mean.csv')

root.abun<-data.frame(t(otu_table(root.4)))

H.root<- diversity(root.abun)
pi.root <- as.data.frame(H.root/log(specnumber(root.abun)))
names(pi.root)[1] <- "pi"

pi.root.1<-na.omit(merge(root.sam, pi.root, by='row.names'))
names(pi.root.1)

r.pi.mean<-ddply(pi.root.1, "Week", summarise,
                 N    = length(pi),
                 mean = mean(pi),
                 sd   = sd(pi),
                 se   = sd / sqrt(N))

write.csv(r.pi.mean, '~/Documents/r.pi.mean.csv')

l.simp.mean<- ddply(leaf.div.1, "Week", summarise,
                    N    = length(InvSimpson),
                    mean = mean(InvSimpson),
                    sd   = sd(InvSimpson),
                    se   = sd / sqrt(N))
write.csv(l.simp.mean, '~/Documents/l.simp.mean.csv')

l.ACE.mean<- ddply(leaf.div.1, "Week", summarise,
                   N    = length(ACE),
                   mean = mean(ACE),
                   sd   = sd(ACE),
                   se   = sd / sqrt(N))

write.csv(l.ACE.mean, '~/Documents/l.ACE.mean.csv')

leaf.abun<-data.frame(t(otu_table(leaf)))

H.leaf<- diversity(leaf.abun)
pi.leaf <- as.data.frame(H.leaf/log(specnumber(leaf.abun)))
names(pi.leaf)[1] <- "pi"

write.csv(l.pi.mean, '~/Documents/leaf.pi.csv')

pi.leaf.1<-read.csv('~/Documents/canola/leaf.pi.csv')

l.pi.mean<-ddply(pi.leaf.1, "week", summarise,
                 N    = length(pi),
                 mean = mean(pi),
                 sd   = sd(pi),
                 se   = sd / sqrt(N))  



################## NRI/NTI #############################
library(picante)

### Soil ###

phydist.soil <- cophenetic(soil.tree.sub)

#week 1
w1.s<-subset_samples(soil.4, Week=='1')
w1.s.abun<-data.frame(otu_table(w1.s))
w1.s.abun.t<-as.data.frame(t(w1.s.abun))

w1.mntd <- ses.mntd(w1.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

#week 2
w2.s<-subset_samples(soil.4, Week=='2')
w2.s.abun<-data.frame(otu_table(w2.s))
w2.s.abun.t<-as.data.frame(t(w2.s.abun))

w2.mntd <- ses.mntd(w2.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

#week 3
w3.s<-subset_samples(soil.4, Week=='3')
w3.s.abun<-data.frame(otu_table(w3.s))
w3.s.abun.t<-as.data.frame(t(w3.s.abun))

w3.mntd <- ses.mntd(w3.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

#week 4
w4.s<-subset_samples(soil.4, Week=='4')
w4.s.abun<-data.frame(otu_table(w4.s))
w4.s.abun.t<-as.data.frame(t(w4.s.abun))

w4.mntd <- ses.mntd(w4.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)
#week 5
w5.s<-subset_samples(soil.4, Week=='5')
w5.s.abun<-data.frame(otu_table(w5.s))
w5.s.abun.t<-as.data.frame(t(w5.s.abun))

w5.mntd <- ses.mntd(w5.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

#week 6
w6.s<-subset_samples(soil.4, Week=='6')
w6.s.abun<-data.frame(otu_table(w6.s))
w6.s.abun.t<-as.data.frame(t(w6.s.abun))

w6.mntd <- ses.mntd(w6.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

#week 7
w7.s<-subset_samples(soil.4, Week=='7')
w7.s.abun<-data.frame(otu_table(w7.s))
w7.s.abun.t<-as.data.frame(t(w7.s.abun))

w7.mntd <- ses.mntd(w7.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

#week 8
w8.s<-subset_samples(soil.4, Week=='8')
w8.s.abun<-data.frame(otu_table(w8.s))
w8.s.abun.t<-as.data.frame(t(w8.s.abun))

w8.mntd <- ses.mntd(w8.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

#week 9
w9.s<-subset_samples(soil.4, Week=='9')
w9.s.abun<-data.frame(otu_table(w9.s))
w9.s.abun.t<-as.data.frame(t(w9.s.abun))

w9.mntd <- ses.mntd(w9.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

#week 10
w10.s<-subset_samples(soil.4, Week=='10')
w10.s.abun<-data.frame(otu_table(w10.s))
w10.s.abun.t<-as.data.frame(t(w10.s.abun))

w10.mntd <- ses.mntd(w10.s.abun.t, phydist.soil, null.model="taxa.labels", 
                     abundance.weighted = TRUE, runs=999)




w1.mpd.s <- ses.mpd(w1.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w2.mpd.s <- ses.mpd(w2.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w3.mpd.s <- ses.mpd(w3.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w4.mpd.s <- ses.mpd(w4.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w5.mpd.s <- ses.mpd(w5.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w6.mpd.s <- ses.mpd(w6.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w7.mpd.s <- ses.mpd(w7.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w8.mpd.s <- ses.mpd(w8.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w9.mpd.s <- ses.mpd(w9.s.abun.t, phydist.soil, null.model="taxa.labels", 
                    abundance.weighted = TRUE, runs=999)

w10.mpd.s <- ses.mpd(w10.s.abun.t, phydist.soil, null.model="taxa.labels", 
                     abundance.weighted = TRUE, runs=999)



####### Root and Leaf NTI/NRI was calcalculated in an indetical manner with the appropriate bacterial abundances

########### BMNTD #############
phylo<-root.tree.sub

#soil w1
otu<-w1.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.1.weighted_bNTI.csv",quote=F);

#week 2
otu<-w2.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.2.weighted_bNTI.csv",quote=F);

### week 3

otu<-w3.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.3.weighted_bNTI.csv",quote=F);

### week 4

otu<-w4.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.4.weighted_bNTI.csv",quote=F);

#week 5

otu<-w5.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.5.weighted_bNTI.csv",quote=F);

#week 6
otu<-w6.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.6.weighted_bNTI.csv",quote=F);


### week 7

otu<-w7.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.7.weighted_bNTI.csv",quote=F);


### week 8

otu<-w8.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.8.weighted_bNTI.csv",quote=F);

#week 9 ####################3
otu<-w9.s.abun.t

match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.9.weighted_bNTI.csv",quote=F);

###### week 10 ###################3

otu<-w10.s.abun.t

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"~/Documents/canola/root/s.10.weighted_bNTI.csv",quote=F);


####### Root and Leaf BMNTD was calcalculated in an indetical manner with the appropriate bacterial abundances


######################### RCBray ################################
library(iCAMP)


#### Leaf ####

l.w1.bray<-as.data.frame(RC.pc(w1.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w2.bray<-as.data.frame(RC.pc(w2.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w3.bray<-as.data.frame(RC.pc(w3.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w4.bray<-as.data.frame(RC.pc(w4.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w5.bray<-as.data.frame(RC.pc(w5.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w6.bray<-as.data.frame(RC.pc(w6.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w7.bray<-as.data.frame(RC.pc(w7.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w8.bray<-as.data.frame(RC.pc(w8.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w9.bray<-as.data.frame(RC.pc(w9.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                               meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                               detail.null=FALSE,output.bray=FALSE,silent=FALSE))

l.w10.bray<-as.data.frame(RC.pc(w10.abun.t, rand = 1000, na.zero = TRUE, weighted = TRUE,
                                meta.ab = NULL,sig.index=c("RC","Confidence","SES"),
                                detail.null=FALSE,output.bray=FALSE,silent=FALSE))




####### Root and Rhizopshere RCBray was calculated in the same manner using the approprite abundances 

