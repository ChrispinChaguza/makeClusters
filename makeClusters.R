#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))
 
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Input similarity file", metavar="character"),
    make_option(c("-t", "--thres"), type="character", default=NULL,
              help="Similarity threshold for clustering", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.lin.tsv", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if("file" %in% names(opt)){
    if("file" %in% names(opt)){
        if(file.exists(opt$file)){
            invisible("OK")
        }else{
            print("Input fastANI file not found")
            quit()
        } 
    }else{
        print("Specify output file name with --out <file>")
        quit()
    }
}else{
    print("Specify input fastANI file with --file <file>")
    quit()
}

strain.df.final.tmp<-as_tibble(fread(opt$file,header=TRUE,sep="\t")) %>%
    dplyr::select(1:3)
colnames(strain.df.final.tmp)<-c("seq1","seq2","pid")

strain.df.final<-strain.df.final.tmp %>%
    dplyr::rename(reference="seq1",query="seq2",ani="pid") %>%
    rowwise() %>%
    mutate(query=gsub("\\.fasta","",query),
         reference=gsub("\\.fasta","",reference)) #%>% dplyr::select(-V4,-V5)

set.seed(1)

edges<-strain.df.final %>% mutate(weight1=ani) %>%
    dplyr::rename(from=reference,to=query) %>% dplyr::select(from,to,weight1) %>%
    mutate(weight=ifelse(weight1>=as.numeric(opt$thres),1,0)) %>% 
    dplyr::select(-weight1)

edges1<-edges %>% dplyr::filter(weight==1) %>% 
    mutate(XX=ifelse(from>to,paste0(from,to),paste0(to,from))) %>% 
    group_by(XX) %>% dplyr::filter(row_number()==1) %>% 
    ungroup() %>% dplyr::select(-XX)

nodes1<-data.frame(name=unique(c(edges1$from,edges1$to))) %>% as_tibble() 

net1<-graph_from_data_frame(d=edges1,vertices=nodes1,directed=FALSE) 

g<-as_undirected(igraph::simplify(net1,remove.loops=TRUE,remove.multiple=TRUE))

strains.cls<-igraph::groups(components(g))

tmp.df.final<-as_tibble(data.frame(cluster.id=paste0("CLS-",as.numeric(names(strains.cls))),
                    cluster.size=sapply(names(strains.cls), function(i){ length(strains.cls[[i]]) }),
                    cluster.strains=sapply(names(strains.cls), function(i){ paste0(strains.cls[[i]],collapse=",") }),
                    cluster.rep.strain=sapply(strains.cls, function(i){ i[[1]][1] }))) %>%
  dplyr::arrange(desc(cluster.size)) %>%
  ungroup() %>% mutate(cluster.id=paste0("CLS-",row_number())) 

write.table(tmp.df.final,file=opt$out,sep="\t",quote=FALSE,row.names=FALSE)

write_graph(graph=g,file=paste0(opt$out,".graphml"),format="graphml")

