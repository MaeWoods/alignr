#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Code to pull the full length TCR sequence from 10X genomics all_contig_annotations.json file
#'
#' This function will import a set of single cell barcodes, extract the full length TCR transcript details
#' and output the details ranked by the cell with the greatest functional score.
#'
#' @param barcodes Barcodes used to identify cells in full contig annotations JSON file.
#' @param contig.annotations contig_annotations.csv file.
#' @param json.path Directory where JSON file is saved.
#' @param save.dir Directory where full length TCR transcripts will be saved.
#' @param score Orders the cells be classification score. 
#' @param verbose Print progress bars and output
#' @return Void
#' @concept Genomics
#'
#' @export
GetTCRs <- function(
    barcodes,
    contig.annotations,
    json.path=".",
    save.dir=".",
    score="",
    ranked.sheet=FALSE,
    TCRs=c(1),
    verbose = TRUE
){
  
  dir.create(file.path(mainDir, "Cloning"), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  result <- fromJSON(file = json.path)
  for(s in 1:length(barcodes)){
    CHR1 <- lapply(result, function(x) { x$barcode == barcodes[s] })
    
    len.cells=dim(contig.annotations)[1]
    indicies=0
    for(k in 1:len.cells){
      if(length(CHR1[[k]])==0){
        
      }
      else{
        
        if(CHR1[k]==TRUE){
          indicies=append(indicies,k)
        }
      }
    }
    indicies=indicies[-1]
  length.fragments=length(indicies)
    
  frag.number = vector(mode = "list", length = length.fragments)
    for(k in 1:length.fragments){
      frag.number[[k]]=result[[indicies[k]]]
      capture.output(frag.number[[k]], file = paste(save.dir,paste(paste(paste("Cell",s,sep=""),"_Chain",sep=""),k,sep=""),".txt",sep=""), append = TRUE)
    }
  
  }

}

#' This function will import a set of single cell barcodes, extract the full length TCR transcript details
#' and output the details ranked by the cell with the greatest functional score.
#'
#' @param barcodes Barcodes used to identify cells in full contig annotations JSON file.
#' @param contig.annotations contig_annotations.csv file.
#' @param json.path Directory where JSON file is saved.
#' @param save.dir Directory where full length TCR transcripts will be saved.
#' @param score Orders the cells be classification score. 
#' @param verbose Print progress bars and output
#' @return Void
#' @concept Genomics
#'
#' @export
RankTCRs <- function(
    cell.data,
    gene="IFNG",
    meta.data=FALSE,
    data.name="Mdist",
    save.dir=".",
    verbose = TRUE
){
  
  if(packageVersion("Seurat")<'5.0.0'){
    RnaStoreUMO=cell.data[['RNA']]$counts
    gene_id=match(gene,row.names(cell.data[['RNA']]@counts))
    clones=levels(factor(cell.data@meta.data$cdr3_na))
    out_df_pre=data.frame(barcodes=rep(0,length(clones)),
                          clonotypes=clones,
                          cdr3=rep(0,length(clones)),
                          max_ele=rep(0,length(clones)),
                          mean_ele=rep(0,length(clones)))
    cell.meta=cell.data@meta.data
    subject.meta.array=rep(0,length(cell.data@meta.data$orig.ident))
    if(meta.data==TRUE){
      pos=match(,names(cell.data@meta.data))
      subject.meta.array=cell.data@meta.data[pos,]
    }
    else{
      subject.meta.array=RnaStoreUMO[gene_id,]
    }
    for(j in 1:length(clones)){
      out_df_pre$cdr3=subset(cell.meta,cdr3_na==clones[j])$cdr3[1]
      
    }
    orderbyfactor=order(out_df_pre$max_ele)
    out_df=data.frame(barcodes=out_df_pre$barcodes[orderbyfactor],
                      clonotypes=out_df_pre$clones[orderbyfactor],
                      cdr3=out_df_pre$cdr3[orderbyfactor],
                      max_ele=out_df_pre$max_ele[orderbyfactor],
                      mean_ele=out_df_pre$mean_ele[orderbyfactor])
    write.csv(out_df,paste(save.dir,"/tcrs_max.csv",sep=""))
    orderbyfactor=order(out_df_pre$mean_ele)
    out_df=data.frame(barcodes=out_df_pre$barcodes[orderbyfactor],
                      clonotypes=out_df_pre$clones[orderbyfactor],
                      cdr3=out_df_pre$cdr3[orderbyfactor],
                      max_ele=out_df_pre$max_ele[orderbyfactor],
                      mean_ele=out_df_pre$mean_ele[orderbyfactor])
    write.csv(out_df,paste(save.dir,"/tcrs_avg.csv",sep=""))
  }
  else{
  RnaStoreUMO=cell.data[['RNA']]$counts
  gene_id=match(gene,row.names(cell.data[['RNA']]$counts))
  clones=levels(factor(cell.data@meta.data$cdr3_na))
  out_df_pre=data.frame(barcodes=rep(0,length(clones)),
                        clonotypes=clones,
                        cdr3=rep(0,length(clones)),
                        max_ele=rep(0,length(clones)),
                        mean_ele=rep(0,length(clones)))
  cell.meta=cell.data@meta.data
  subject.meta.array=rep(0,length(cell.data@meta.data$orig.ident))
  if(meta.data==TRUE){
    
    subject.meta.array=cell.data@meta.data[pos,]
  }
  else{
    subject.meta.array=RnaStoreUMO[gene_id,]
  }
  for(j in 1:length(clones)){
    out_df_pre$cdr3=subset(cell.meta,cdr3_na==clones[j])$cdr3[1]
    
  }
  
  orderbyfactor=order(out_df_pre$max_ele)
  out_df=data.frame(barcodes=out_df_pre$barcodes[orderbyfactor],
                        clonotypes=out_df_pre$clones[orderbyfactor],
                        cdr3=out_df_pre$cdr3[orderbyfactor],
                        max_ele=out_df_pre$max_ele[orderbyfactor],
                        mean_ele=out_df_pre$mean_ele[orderbyfactor])
  write.csv(out_df,paste(save.dir,"/tcrs_max.csv",sep=""))
  orderbyfactor=order(out_df_pre$mean_ele)
  out_df=data.frame(barcodes=out_df_pre$barcodes[orderbyfactor],
                    clonotypes=out_df_pre$clones[orderbyfactor],
                    cdr3=out_df_pre$cdr3[orderbyfactor],
                    max_ele=out_df_pre$max_ele[orderbyfactor],
                    mean_ele=out_df_pre$mean_ele[orderbyfactor])
  write.csv(out_df,paste(save.dir,"/tcrs_avg.csv",sep=""))
  }
  
}

#' This function will import a set of single cell barcodes, extract the full length TCR transcript details
#' and output the details ranked by the cell with the greatest functional score.
#'
#' @param barcodes Barcodes used to identify cells in full contig annotations JSON file.
#' @param contig.annotations contig_annotations.csv file.
#' @param json.path Directory where JSON file is saved.
#' @param save.dir Directory where full length TCR transcripts will be saved.
#' @param score Orders the cells be classification score. 
#' @param verbose Print progress bars and output
#' @return Void
#' @concept Genomics
#'
#' @export
GetFASTQ <- function(
    gene.name,
    n.cells,
    path
){
  
  for(){
  
  TRAC="XIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSN 
SAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGF 
RILLLKVAGFNLLMTLRLWSS"
  
  TRBC2="XDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDP 
QPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQI 
VSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDSRG"
  
  TRBC1="XDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDP 
QPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQI 
VSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDF"
  
  if((JSON_output$X.barcode[141]=="[1] TRUE")&
     (JSON_output$X.barcode[143]=="[1] TRUE")&
     (JSON_output$X.barcode[145]=="[1] TRUE")&
     (JSON_output$X.barcode[147]=="[1] TRUE")&
     (JSON_output$X.barcode[149]=="[1] TRUE")&
     (JSON_output$X.barcode[151]=="[1] TRUE")&
     (JSON_output$X.barcode[153]=="[1] TRUE")){
    print("yes")
  }
  
  if((JSON_output_TRB$X.barcode[141]=="[1] TRUE")&
     (JSON_output_TRB$X.barcode[143]=="[1] TRUE")&
     (JSON_output_TRB$X.barcode[145]=="[1] TRUE")&
     (JSON_output_TRB$X.barcode[147]=="[1] TRUE")&
     (JSON_output_TRB$X.barcode[149]=="[1] TRUE")&
     (JSON_output_TRB$X.barcode[151]=="[1] TRUE")&
     (JSON_output_TRB$X.barcode[153]=="[1] TRUE")){
    print("yes")
  }
  
  if(JSON_output$X.barcode[98]=="TRA"){
  gene.name="TRAC"
  path="/Users/maewoodsphd/mVSTManuscript/TCRCloning/Cloning/EBVFullLengthTCRs/TruePositive/Cell1Cloning/"
  curl_download("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz","/Users/maewoodsphd/mVSTManuscript/TCRCloning/Cloning/EBVFullLengthTCRs/TruePositive/Cell1Cloning/HUMAN_9606_idmapping.dat.gz",mode="w",quiet=FALSE)
  locusID = match("TRAC",read.table(gzfile("/Users/maewoodsphd/mVSTManuscript/TCRCloning/Cloning/EBVFullLengthTCRs/TruePositive/Cell1Cloning/HUMAN_9606_idmapping.dat.gz"), fill = TRUE)$V3)
  PiD = read.table(gzfile("/Users/maewoodsphd/mVSTManuscript/TCRCloning/Cloning/EBVFullLengthTCRs/TruePositive/Cell1Cloning/HUMAN_9606_idmapping.dat.gz"), fill = TRUE)$V1[locusID]
  curl_download(paste(paste("https://rest.uniprot.org/uniprotkb/",PiD,sep=""),
                      ".fasta",sep=""),"bb.txt",mode="w",quiet=FALSE)
  }
  
  }
}
  
  

