# 04-28-2026
# functions for calling ID30 from Koh et al 2023 (Nature Genetics)
# template for ID30 channels is taken from "Supplementary_Table_S2_Signature" from the same manuscript
# one final function written by K.T.S to pad out the matrix (i.e., include channels with a frequency of zero)

# the following functions (and indel template) were written by members of the Nik Zainal lab

call_legacy_indels <- function(df,
                               callers=list(microhomology_mediated_del,repeat_mediated_indel),
                               decision=current_decision){
  df$slice3_1bp_pyr <- substr(df$slice3_pyr,1,1)
  calls <- purrr::map(callers,~ purrr::map2(df[['change']],df[['slice3']],.))
  calls <- purrr::pmap_dfr(c(calls,list(type=df[['indel.type']])),current_decision)
  colnames(calls) <- c('legacy_classification','legacy_score','legacy_seq')
  df <- dplyr::bind_cols(df,calls)
  df <- assign_legacy_channels(df)
  
  return(df)
}

assign_legacy_channels <- function(df){
  ## check colnames
  required_cols <- c('indel.type','indel.length','change_pyr','slice3_1bp_pyr','legacy_score','legacy_classification')
  if(any(!(required_cols %in% colnames(df)))){
    message("required columns missing")
    stop()
  }
  
  #indel.classified_details2$length_type <- indel.classified_details2$indel.length
  #indel.classified_details2[indel.classified_details2$indel.length>1,"length_type"] <- ">1bp"
  #indel.classified_details2$Subtype <- paste0(indel.classified_details2$Type,"_",indel.classified_details2$legacy_classification,"_",indel.classified_details2$length_type)
  #indel.classified_details2[indel.classified_details2$indel.length==1,"Subtype"] <- paste0(indel.classified_details2[indel.classified_details2$indel.length==1,"Subtype"], "_",indel.classified_details2[indel.classified_details2$indel.length==1,"change_pyr"])
  
  df$legacy_channel <- NULL
  df[df$indel.type=="I" & df$change_pyr=="C" & df$legacy_score==0 & df$slice3_1bp_pyr=="A","legacy_channel"] <- "[+C]A"
  df[df$indel.type=="I" & df$change_pyr=="C" & df$legacy_score==0 & df$slice3_1bp_pyr=="G","legacy_channel"] <- "[+C]G"
  df[df$indel.type=="I" & df$change_pyr=="C" & df$legacy_score==0 & df$slice3_1bp_pyr=="T","legacy_channel"] <- "[+C]T"
  df[df$indel.type=="I" & df$change_pyr=="C" & df$legacy_score==1,"legacy_channel"] <- "[+C]C"
  df[df$indel.type=="I" & df$change_pyr=="C" & df$legacy_score==2,"legacy_channel"] <- "[+C]CC"
  df[df$indel.type=="I" & df$change_pyr=="C" & df$legacy_score>2,"legacy_channel"] <- "[+C]LongRep"
  
  df[df$indel.type=="I" & df$change_pyr=="T" & df$legacy_score==0 & df$slice3_1bp_pyr=="A","legacy_channel"] <- "[+T]A"
  df[df$indel.type=="I" & df$change_pyr=="T" & df$legacy_score==0 & df$slice3_1bp_pyr=="C","legacy_channel"] <- "[+T]C"
  df[df$indel.type=="I" & df$change_pyr=="T" & df$legacy_score==0 & df$slice3_1bp_pyr=="G","legacy_channel"] <- "[+T]G"
  df[df$indel.type=="I" & df$change_pyr=="T" & df$legacy_score==1,"legacy_channel"] <- "[+T]T"
  df[df$indel.type=="I" & df$change_pyr=="T" & df$legacy_score==2,"legacy_channel"] <- "[+T]TT"
  df[df$indel.type=="I" & df$change_pyr=="T" & df$legacy_score>2,"legacy_channel"] <- "[+T]LongRep"
  
  df[df$indel.type=="I" & df$indel.length>1 & df$legacy_score<1,"legacy_channel"] <- "[+>1]NonRep"
  df[df$indel.type=="I" & df$indel.length>1 & df$legacy_score>=1,"legacy_channel"] <- "[+>1]Rep"
  
  df[df$indel.type=="D" & df$change_pyr=="C" & df$legacy_score==0 & df$slice3_1bp_pyr=="A","legacy_channel"] <- "[-C]A"
  df[df$indel.type=="D" & df$change_pyr=="C" & df$legacy_score==0 & df$slice3_1bp_pyr=="G","legacy_channel"] <- "[-C]G"
  df[df$indel.type=="D" & df$change_pyr=="C" & df$legacy_score==0 & df$slice3_1bp_pyr=="T","legacy_channel"] <- "[-C]T"
  df[df$indel.type=="D" & df$change_pyr=="C" & df$legacy_score==1,"legacy_channel"] <- "[-C]C"
  df[df$indel.type=="D" & df$change_pyr=="C" & df$legacy_score==2,"legacy_channel"] <- "[-C]CC"
  df[df$indel.type=="D" & df$change_pyr=="C" & df$legacy_score>2,"legacy_channel"] <- "[-C]LongRep"
  
  df[df$indel.type=="D" & df$change_pyr=="T" & df$legacy_score==0 & df$slice3_1bp_pyr=="A","legacy_channel"] <- "[-T]A"
  df[df$indel.type=="D" & df$change_pyr=="T" & df$legacy_score==0 & df$slice3_1bp_pyr=="C","legacy_channel"] <- "[-T]C"
  df[df$indel.type=="D" & df$change_pyr=="T" & df$legacy_score==0 & df$slice3_1bp_pyr=="G","legacy_channel"] <- "[-T]G"
  df[df$indel.type=="D" & df$change_pyr=="T" & df$legacy_score==1,"legacy_channel"] <- "[-T]T"
  df[df$indel.type=="D" & df$change_pyr=="T" & df$legacy_score==2,"legacy_channel"] <- "[-T]TT"
  df[df$indel.type=="D" & df$change_pyr=="T" & df$legacy_score>2,"legacy_channel"] <- "[-T]LongRep"
  
  df[df$indel.type=="D" & df$legacy_classification !="Microhomology-mediated" & df$indel.length>1 & df$legacy_score<1,"legacy_channel"] <- "[->1]Others"
  df[df$indel.type=="D" & df$legacy_classification=="Repeat-mediated" & df$indel.length>1 & df$legacy_score>=1,"legacy_channel"] <- "[->1]Rep"
  df[df$indel.type=="D" & df$legacy_classification =="Microhomology-mediated","legacy_channel"] <- "[->1]Mh"
  df[df$indel.type=="DI","legacy_channel"] <- "[Com]"
  
  return(df) 
}

current_decision <- function(type,...){
  res <- list(...)
  calls <- unlist(map(res,magrittr::extract,'call'))
  scores <- unlist(map(res,magrittr::extract,'score'))
  seqs <- unlist(map(res,magrittr::extract,'seq'))
  names(scores)<-names(seqs)<-names(res)<-calls
  
  if(type=='DI'){
    return(list(call='Unassigned',score=-1,seq='Unassigned'))
  }
  else if(type=='D'){
    if(all(scores!=0)){
      if(scores['Repeat-mediated'] * nchar(seqs['Repeat-mediated']) >= scores['Microhomology-mediated']){
        return(res[['Repeat-mediated']])
      }else{
        return(res[['Microhomology-mediated']])
      }
    } else if(sum(scores) > 0){
      return(res[[which(scores!=0)]])
    } else{
      return(list(call='None',score=0,seq=''))
    }
  }
  else if(type=='I'){
    if(scores['Repeat-mediated'] > 0){
      return(res[['Repeat-mediated']])
    } else{
      return(list(call='None',score=0,seq=''))
    }
  }
  else{
    return(list(call='Failed',score=0,seq=''))
  }
}

left_shared_subsequence <- function(seq1,seq2){
  for(i in 1L:min(nchar(seq1),nchar(seq2))){
    if(str_sub(seq1,1L,i) != str_sub(seq2,1L,i)){
      i <- i - 1L
      break()
    }
  }
  return(substr(seq1,1L,i))
}

left_shared_tandem_subsequence <- function(seq1,seq2){
  pat <- str_c('^(',seq1,')+')
  res <-str_match(seq2,pat)[1]
  res <- ifelse(is.na(res),'',res)
  return(res)
}

microhomology_mediated_del <- function(change,slice3){
  ## Use  helper function to find homology between the deletion and 3' context from the left
  mh_seq <- left_shared_subsequence(as.character(change),as.character(slice3))
  result <- list(call='Microhomology-mediated',score=nchar(mh_seq),seq=mh_seq)
  return(result)
}

repeat_mediated_indel <- function(change,slice3){
  smallest_rep_unit <- smallest_exact_repetitive_subsequence(change)
  if(nchar(change) == 1){
    rpt_seq <- left_shared_tandem_subsequence(change,slice3)
    result <- list(score=nchar(rpt_seq),seq=rpt_seq)
  }
  else if(smallest_rep_unit == str_sub(slice3,1,nchar(smallest_rep_unit))){ ## catches homopolymers.
    match <- left_shared_tandem_subsequence(smallest_rep_unit,slice3)
    result <- list(score=nchar(match)/nchar(smallest_rep_unit),seq=match)
  }
  else {
    srs <- smallest_repetitive_subsequence(seq = change,min_size = 1)
    
    if(srs[2] == change){
      result <- list(score=0,seq='')
    }
    else if(nchar(srs[1])/nchar(change) > 2/3){ ## the repetitive sequence is 2/3rds the indel.
      match <- left_shared_tandem_subsequence(srs[2],slice3)
      
      tot_matches <- (nchar(srs[1]) + nchar(match))/nchar(srs[2])
      
      if(tot_matches >= 3 & nchar(match) > 1){
        result <- list(score=nchar(match)/nchar(srs[2]),seq=match)
      }
      else{
        result <- list(score=0,seq='')
      }
    }
    else{
      result <- list(score=0,seq='')
    }
  }
  return(c(call='Repeat-mediated',result))
}

smallest_exact_repetitive_subsequence <- function(seq){
  len <- as.integer(nchar(seq))
  factors <- seq_len(len)[ len %% seq_len(len) == 0L]
  for(i in seq_along(factors)){
    if(str_dup(str_sub(seq,1L,factors[i]),len/factors[i]) == seq){
      break()
    }
  }
  return(str_sub(seq,1L,factors[i]))
}

smallest_repetitive_subsequence <- function(seq,min_size=1){
  for( i in nchar(seq):min_size){
    res <- str_match(seq,str_c('^(.{',i,'})\\1+'))
    if(!is.na(res[1])){
      pr2 <- smallest_repetitive_subsequence(res[2],min_size) #unique(str_split(res[2],'')[[1]])
      res[2] <- ifelse(pr2[2]==res[2],res[2],pr2[2])
      return(res)
    }
  }
  return(c(seq,seq))
}

# plotting function written by members of the Nik-Zainal lab, but edited by K.T.S 
  # input to this function is output of: get_matrix_ID30() (as below)
gen_plot_catalouge_legacy_single <- function(muts_basis,text_size,plot_title){
  
  # define the indel template
  indel_template_legacy <- data.frame("IndelType"=c("[+C]A", "[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LongRep",
                                                    "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LongRep",
                                                    "[+>1]NonRep","[+>1]Rep",
                                                    "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LongRep",
                                                    "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LongRep",
                                                    "[->1]Others","[->1]Rep",
                                                    "[->1]Mh",
                                                    "[Com]"),
                                      "Indel"=c(
                                        rep("Ins(C)",6),
                                        rep("Ins(T)",6),
                                        rep('Ins(2,)',2),
                                        rep("Del(C)",6),
                                        rep("Del(T)",6),
                                        rep('Del(2,)',2),'Del(Mh)',
                                        'Complex')
  )
  
  # prepare the df for input
  # first, add rownames back as a column
  muts_basis$IndelType <- rownames(muts_basis)
  
  # then, reorder rows as in the indel_template_legacy df
  muts_basis <- muts_basis[indel_template_legacy$IndelType,]
  
  # then, add indel column from the indel_template_legacy df
  muts_basis$Indel <- indel_template_legacy$Indel
  
  # add a sample column
  muts_basis$Sample <- rep(colnames(muts_basis)[1], 30)
  
  # rename the frequency column
  colnames(muts_basis)[1] <- "freq"
  
  # remove rownames
  rownames(muts_basis) <- NULL
  
  # re-order the columns
  muts_basis <- muts_basis %>% dplyr::select(IndelType, Indel, Sample, freq)
  
  # rename variable so we don't need to edit the rest of the function
  muts_basis_melt = muts_basis
  
  # define colour palettes
  indel_mypalette_fill2 <- c("#000000", # FEABB9 Complex
                             "#EE3377", # Del(2,)
                             "#004488", # Del(C)
                             "#762A83", # Mh
                             "#997700", # Del(T)
                             "#EE99AA", #"#668D3C"  Ins(2,)
                             "#6699CC", #"#007996"  Ins(C)
                             "#EECC66") # Ins(T)
  indel_mypalette_fill <- c("#000000", # FEABB9 Complex
                            "#009E73", # Del(2,)
                            "#56B4E9", # Del(C)
                            "#762A83", # Mh
                            "#E69F00", # Del(T)
                            "#CC79A7", #"#668D3C"  Ins(2,)
                            "#0072B2", #"#007996"  Ins(C)
                            "#D55E00") # Ins(T)
  
  indel_positions <- indel_template_legacy$IndelType
  entry <- table(indel_template_legacy$Indel)
  order_entry <- c("Ins(C)", "Ins(T)", "Ins(2,)", "Del(C)", "Del(T)", "Del(2,)", "Del(Mh)", "Complex")
  
  
  entry <- entry[order_entry]
  blocks <- data.frame(Type=unique(indel_template_legacy$Indel),
                       fill=indel_mypalette_fill,
                       xmin=c(0,cumsum(entry)[-length(entry)])+0.5,
                       xmax=cumsum(entry)+0.5)
  blocks$ymin <- max(muts_basis_melt$freq)*1.08
  blocks$ymax <- max(muts_basis_melt$freq)*1.2
  #  blocks$labels <-c("+C", "+T", "+M", "+N", "-C", "-T", "-M", "-N", "-Mh", "X")
  #  blocks$cl <-c("white", "white", "white", "white", "black", "black", "black", "black", "black", "white")
  
  blocks$labels <-c("1bp C", "1bp T", ">=2bp", "1bp C", "1bp T", ">=2bp", "Mh", "X")
  blocks$cl <-c("white", "white", "white","black", "black", "black", "white",  "white")
  
  
  p <- ggplot2::ggplot(data=muts_basis_melt, ggplot2::aes(x=IndelType, y=freq,fill=Indel))+ 
    ggplot2::geom_bar(stat="identity",position="dodge", width=.8)+
    ggplot2::xlab("Indel Types")+ggplot2::ylab("Count")
  p <- p+ggplot2::scale_x_discrete(limits = indel_positions)+ ggplot2::ggtitle(plot_title)
  p <- p+ggplot2::scale_fill_manual(values=indel_mypalette_fill)+ggplot2::coord_cartesian(ylim=c(0,unique(blocks$ymax)), expand = FALSE)
  p <- p+ggplot2::theme_classic()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
                                                 axis.text.y=ggplot2::element_text(size=10,colour = "black"),
                                                 legend.position = "none",
                                                 axis.title.x = ggplot2::element_text(size=15),
                                                 axis.title.y = ggplot2::element_text(size=15))
  
  ## Add the overhead blocks
  p <- p+ggplot2::geom_rect(data = blocks, ggplot2::aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,fill=Type),inherit.aes = F)+
    ggplot2::geom_text(data=blocks,ggplot2::aes(x=(xmax+xmin)/2,y=(ymax+ymin)/2,label=labels, colour=cl),size=text_size,fontface="bold",inherit.aes = F)+ggplot2::scale_colour_manual(values=c("black", "white"))
  
  # return plot
  return(p)

}

# the following function is written by K.T.S, as a complement to the functions above and to generate the full matrix
  # df here is the output from "call_legacy_indels" function
  # function supports multiple sample input
get_matrix_ID30 <- function(df) { 
  
  # channels:
  channels <- data.frame(MutationType = c("[->1]Mh", "[->1]Others", "[->1]Rep", "[-C]A", "[-C]C", "[-C]CC", "[-C]G", 
                                          "[-C]LongRep", "[-C]T", "[-T]A", "[-T]C", "[-T]G", "[-T]LongRep", "[-T]T", 
                                          "[-T]TT", "[+>1]NonRep", "[+>1]Rep", "[+C]A", "[+C]C", "[+C]CC", "[+C]G", 
                                          "[+C]LongRep", "[+C]T", "[+T]A", "[+T]C", "[+T]G", "[+T]LongRep", "[+T]T", 
                                          "[+T]TT", "[Com]"))
  
  # 
  sample_list <- df$Sample |> unique()
  
  # 
  for (s in 1:length(sample_list)) {
    
    # subset for sample
    indels <- df %>% dplyr::filter(Sample %in% sample_list[s])
    
    # 
    matrix <- indels$legacy_channel |> table() |> as.data.frame()
    
    # rename columns
    colnames(matrix) <- c("MutationType", sample_list[s])
    
    # 
    matrix_final <- dplyr::left_join(x = channels, y = matrix, by = "MutationType") %>%
      mutate(across(where(is.numeric), ~replace_na(.x, 0)))
    
    # assign mutation type column as rownames
    rownames(matrix_final) <- matrix_final$MutationType
    
    # remove mutation type column
    matrix_final$MutationType <- NULL
    
    # 
    if (exists("final_df") == TRUE) {
      final_df <- cbind(final_df, matrix_final)
    } else {
      final_df <- matrix_final
    }
  }
  
  return(final_df)
}






                                    