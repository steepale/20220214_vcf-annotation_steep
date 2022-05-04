#' ExtractNamedAbundance1A
#' A function to create a matrix with viallabels (or bids) and metabolite names for named metabolites
#' @param export_df the named metabolite data from the UMICH metabolomics core
#' @param manifest_df the sample info from the UMICH metabolomics core
#' @param transpose logical: whether to transpose the dataframe. If transpose = FALSE, then metabolites and arrange by row and samples by column
#' @param sample_id which sample identifier to use for labeling experimental samples (excludes references). Choices include: "viallabel" or "bid"
#' @param format_out the object output format. One of "df" or "mat"
#' @return an object with viallabels and metabolite names for named metabolites
# Function to load vep annotated germline vcfs
################################################################################
load_vep_vcf <- function(x) {
        # read in dependencies
        library("ff")
        # Ensure 'x' is a string
        if (!"character" %in% class(x)  ) {
                stop("'x' must be a string", class.= FALSE)
        }
        
        # Gunzip the file if necessary
        if (endsWith(x, '.gz')) {
                compressed <- TRUE
        }else{
                compressed <- FALSE 
        }
        
        if (compressed) {
                # Determine the number of lines to skip (with ##)
                system_cmd <- paste0('bgzip -d -c ',x,' | grep "^##" | wc -l')
                n_skip <- system(system_cmd, intern = TRUE) %>% as.numeric()
                # Determine the columns
                system_cmd <- paste0('bgzip -d -c ',x,' | grep "^#CHROM" | tr "\t" "\n" | wc -l')
                col_n <- system(system_cmd, intern = TRUE) %>% as.numeric()
        }else{
                # Determine the number of lines to skip (with ##)
                system_cmd <- paste0('grep "^##" ',x,'| wc -l')
                n_skip <- system(system_cmd, intern = TRUE) %>% as.numeric()
                # Determine the columns
                system_cmd <- paste0('grep "^#CHROM" ',x,'| tr "\t" "\n" | wc -l')
                col_n <- system(system_cmd, intern = TRUE) %>% as.numeric()
        }
        
        # Detemrine the column class for file loading
        col_class <- c("factor", "integer",rep("factor",3),"double",rep("factor", (col_n - 6) ))
        
        # Load vcf data
        #PL <- read.table(file = x, comment.char = '', skip = n_skip, check.names = FALSE, header = TRUE, sep = '\t') %>% as_tibble()
        PL <- read.table.ffdf(file = x, sep = '\t', skip=n_skip,header=TRUE,
                              colClasses=col_class, VERBOSE = TRUE, first.rows = 1000000 ) %>% as_tibble()
        # Adjust col name
        colnames(PL)[1] <- 'CHROM'
        PL$sample <- colnames(PL)[10]
        
        ##### Format the FORMAT column
        #################################
        format <- unique(as.character(PL$FORMAT))
        form_vec <- list()
        for(i in 1:length(format)){
          form_vec[[i]] <- format[i] %>% str_split(pattern = ':') %>% unlist()
        }
        
        PL_list <- list()
        PL2 <- data.frame()
        for(i in 1:length(format)){
          PL_list[[i]] <- PL %>% filter(FORMAT == format[[i]])
          PL_list[[i]] <- PL_list[[i]] %>% separate(col = colnames(PL)[10], into = form_vec[[i]] , sep = ':') %>% select(-FORMAT)
          PL2 <- dplyr::bind_rows(PL2, PL_list[[i]])
        }
        # Format the VEP column in the INFO column
        ####################################
        PL2 <- PL2 %>%
          separate(col = INFO, into = c('INFO', 'vep1'), sep = "\\|", extra = "merge") %>%
          separate(col = vep1, into = c('1', '2', '3','4','5','6','7','8','9','10'), sep = ",") %>%
          pivot_longer(cols = c('1', '2', '3','4','5','6','7','8','9','10'), names_to = "VEP_ANN_N", values_to = "VEP") %>%
          mutate(VEP_ANN_N = as.numeric(VEP_ANN_N)) %>%
          filter(!is.na(VEP)) %>%
          separate(ALT, into = c('1','2','3'), sep = ',') %>%
          pivot_longer(cols = c('1','2','3'), names_to = "ALT_N", values_to = "ALT") %>%
          mutate(ALT_N = as.numeric(ALT_N)) %>%
          filter(!is.na(ALT)) %>%
          separate(AD, into = c('AD_REF','1','2', '3'), sep = ',') %>%
          group_by(CHROM,POS,REF,ALT) %>%
          pivot_longer(cols = c('1','2','3'), names_to = "AD_ALT_N", values_to = "AD_ALT") %>%
          ungroup() %>% unique() %>%
          filter(AD_ALT_N == ALT_N) %>%
          group_by(CHROM,POS,REF) %>%
          mutate(ALT_SUM = max(ALT_N, na.rm = TRUE)) %>%
          ungroup() %>%
          mutate(VEP = ifelse(VEP_ANN_N == 1, paste0(ALT,'|',VEP), VEP)) %>%
          # group_by(CHROM,POS,REF) %>%
          # mutate(VEP_ANN_N = max(VEP_ANN_N, na.rm = TRUE)) %>%
          # ungroup() %>%
          separate(col= VEP, into = c('var','var_type', 'var_impact','gene_symbol',
                                      'ensembl_geneid','transcript','ensembl_transcript','transcript_type','exon1','exon2',
                                      'transcript_mut','protein_mut', 'transcript_mut_pos2', 'transcript_mut_pos1',
                                      'protein_mut1','delta_AA','delta_codon','c17','mut_pos_alt','strand','c20',
                                      'ann_database', 'HGNC_protein_n', 'c23','c24','c25','ensembl_protein','PaxDb_ID_1',
                                      'PaxDb_ID_2','UK1','SIFT','protein_domains'), sep = "\\|") %>%
          filter(ALT == var) %>% unique()
        
        ###############################
        # Format the INFO column
        ###############################
        # Make another format column
        PL2 <- PL2 %>%
          rowwise() %>%
          mutate(FORMAT = str_replace_all(INFO, c("=.*;AF"= ";AF"))) %>%
          mutate(FORMAT = str_replace_all(FORMAT, c("=.*;AN"= ";AN",
                                                    "=.*;BaseQRankSum"= ";BaseQRankSum",
                                                    "=.*;ClippingRankSum"= ";ClippingRankSum",
                                                    "=.*;DP"= ";DP",
                                                    "=.*;ExcessHet"= ";ExcessHet",
                                                    "=.*;FS"= ";FS",
                                                    "=.*;InbreedingCoeff"= ";InbreedingCoeff",
                                                    "=.*;MQ="= ";MQ=",
                                                    "=.*;MQRankSum="= ";MQRankSum=",
                                                    "=.*;QD"= ";QD",
                                                    "=.*;ReadPosRankSum"= ";ReadPosRankSum",
                                                    "=.*;SOR"= ";SOR",
                                                    "=.*;CSQ"= ";CSQ",
                                                    ";CSQ=[A-Z]"= ";CSQ"))) %>%
          ungroup()
        # Make a second info column
        PL2 <- PL2 %>%
          rowwise() %>%
          mutate(INFO2 = str_remove_all(INFO, "AC=")) %>%
          mutate(INFO2 = str_remove_all(INFO2, "AF=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "AN=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "BaseQRankSum=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "ClippingRankSum=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "DP=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "ExcessHet=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "FS=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "InbreedingCoeff=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "MQRankSum=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "MQ=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "QD=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "ReadPosRankSum=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "SOR=")) %>%
            mutate(INFO2 = str_remove_all(INFO2, "CSQ=")) %>%
          ungroup()
        
        format <- unique(as.character(PL2$FORMAT))
        form_vec <- list()
        for(i in 1:length(format)){
          form_vec[[i]] <- format[i] %>% str_split(pattern = ';') %>% unlist()
        }
        PL_list <- list()
        PL3 <- data.frame()
        i <- 1
        for(i in 1:length(format)){
          PL_list[[i]] <- PL2 %>% filter(FORMAT == format[[i]])
          PL_list[[i]] <- PL_list[[i]] %>% separate(col = INFO2, into = form_vec[[i]] , sep = ';') %>% select(-FORMAT)
          PL3 <- dplyr::bind_rows(PL3, PL_list[[i]])
        }
        ##################################
        # Final adjustments
        PL3 <- PL3 %>%
          filter(ALT == var) %>%
          select(sample,CHROM,POS,ID,REF,ALT,ALT_N,ALT_SUM,AD_REF, AD_ALT, AD_ALT_N, QUAL,FILTER,
                 GT,GQ,PL,PGT,PID,
                 AC,AF,AN,BaseQRankSum,ClippingRankSum,DP,ExcessHet,FS,InbreedingCoeff,
                 MQ,MQRankSum,QD,ReadPosRankSum,CSQ,
                 VEP_ANN_N,var_type,var_impact,gene_symbol,ensembl_geneid,transcript,ensembl_transcript,
                 transcript_type,exon1,exon2,transcript_mut,protein_mut,transcript_mut_pos2,
                 transcript_mut_pos1,protein_mut1,delta_AA,delta_codon,c17,mut_pos_alt,strand,
                 c20,ann_database,HGNC_protein_n,c23,c24,c25,ensembl_protein,PaxDb_ID_1,PaxDb_ID_2,
                 UK1,SIFT,protein_domains)
        bk <-PL3
        
        unique(bk$AF)
        
        # Adjust the classes
        PL3$GQ <- as.numeric(PL3$GQ)
        PL3$AD_REF <- as.numeric(PL3$AD_REF)
        PL3$AD_ALT <- as.numeric(PL3$AD_ALT)
        PL3$AD_ALT_N <- as.numeric(PL3$AD_ALT_N)
        PL3$AN <- as.numeric(PL3$AN)
        PL3$BaseQRankSum <- as.numeric(PL3$BaseQRankSum)
        PL3$ClippingRankSum <- as.numeric(PL3$ClippingRankSum)
        PL3$DP <- as.numeric(PL3$DP)
        PL3$ExcessHet <- as.numeric(PL3$ExcessHet)
        PL3$FS <- as.numeric(PL3$FS)
        PL3$InbreedingCoeff <- as.numeric(PL3$InbreedingCoeff)
        PL3$MQ <- as.numeric(PL3$MQ)
        PL3$MQRankSum <- as.numeric(PL3$MQRankSum)
        PL3$QD <- as.numeric(PL3$QD)
        PL3$ReadPosRankSum <- as.numeric(PL3$ReadPosRankSum)
        PL3$transcript_mut_pos2 <- as.numeric(PL3$transcript_mut_pos2)
        PL3$transcript_mut_pos1 <- as.numeric(PL3$transcript_mut_pos1)
        PL3$protein_mut1 <- as.numeric(PL3$protein_mut1)

        # Return the output
        PL3
}