Paras =  readLines("./MRmethod/allParas.txt")
path = Paras[1]
savePath = Paras[2]

# ??Folders?еõ???ͬpower???ļ?????
i = 0
fileCup = data.frame()

Folders = list.files(path)
for(filename in Folders){
    i = i + 1
    fileCup[i,"file"] = filename
    
    S = stringr::str_split(filename,"_")[[1]]
    fileCup[i,"H"] = S[1]
    fileCup[i,"gamma"] = S[2]
    fileCup[i,"tau"] = S[3]
    fileCup[i,"q"] = S[4]
    fileCup[i,"eta"] = S[5]
    fileCup[i,"omega"] = S[6]
    fileCup[i,"sample1"] = S[7]
    fileCup[i,"sample2"] = S[8]
    fileCup[i,"n1"] = S[9]
    fileCup[i,"n2"] = S[10]
    fileCup[i,"h1"] = S[11]
    fileCup[i,"h2"] = S[12]
}

fileCup = fileCup[order(fileCup$q),]
row.names(fileCup) = seq(dim(fileCup)[1])



# filter 
S = sort(fileCup[(fileCup$sample1 == 40000 )&(fileCup$sample2 == 40000 ) ,1])

write.table(S,"./MRmethod/all_runFiles.txt",quote = F,row.names = F,col.names = F)