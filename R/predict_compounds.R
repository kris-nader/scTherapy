######################################################################################################
##                         Get functions for compound prediction                                    ##
######################################################################################################
#
# GNU General Public License v3.0 (https://github.com/kris-nader/scTherapy/blob/main/LICENSE)
#
# Written by Kristen Michelle Nader <kristen.nader@helsinki.fi>
# and Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi> , April 2021
#
# Functions on this page:
# process_DEG,predict_drugs
#
#



#' @title filter DEGs
#' @name process_DEG
#' @description filters DEGs based on avg_log2FC and p_val_adj
#' @details none
#' @param degs_list: Differentially expressed genes(see subclone_DEG for more information)
#' @return list of DEGs and avg_log2fc
#' @export
#' @examples
#' diff_expr_list = process_DEG(degs_list = DEG_subclone)
#' 
process_DEG <- function(diff_expr){
    gene_list="https://raw.githubusercontent.com/kris-nader/TBD/main/geneinfo_beta_input.txt"
    gene_info = data.table::fread(gene_list) %>% as.data.frame()
    diff_expr_df <- diff_expr %>%
        mutate(gene_symbol = rownames(.)) %>% inner_join(gene_info, by = "gene_symbol") %>%
        filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
    diff_expr_list <- setNames(as.list(diff_expr_df$avg_log2FC), diff_expr_df$gene_symbol)
    return(diff_expr_list)
}


#' @title predict compounds based on DEGs
#' @name predict_drugs
#' @description Predicts compounds using pre-trained LightGBM model
#' @details none
#' @param degs_list: Differentially expressed genes(see process_DEG for more information)
#' @param exclude_low_confidence
#' @return dataframe of predicted drugs and doses 
#' @export
#' @examples
#' df = predict_drugs(degs_list = diff_expr_list)
#' 

predict_drugs <- function(degs_list, exclude_low_confidence = !0) {
  
  # check internet connection
  response <- httr::HEAD("http://google.com")
  if (!(http_status(response)$category == "Success")){
    stop("Please check your internet connection.")} else {log_info("Connection successful.")}
  
  # check DEGs
  if(!all(sum(unlist(degs_list) > 0) >= 10, sum(unlist(degs_list) < 0) >= 10)){
    stop("There should be at least 10 up- and down- regulated genes.")}
  
  # Convert the named vector to a named list and encode as a JSON string
  degs_list = lapply(degs_list, round, 4)
  json_value <- toJSON(setNames(as.list(degs_list), names(degs_list)))
  
  # Send a POST request to the API with the JSON string as a parameter
  log_info("Making predictions...")
  response <- POST("https://sctype.app/sctypemodel/mdfull.php", body = list(value = json_value))
  
  if (http_status(response)$category != "Success") {
    stop("API request failed")}
  
  log_info("Preparing output...")
  response_content <- sub(".*cid;Drug_Name", "cid;Drug_Name", content(response, "text"))
  
  # Parse the CSV response into a data frame
  suppressWarnings(df <- read.table(text = response_content, sep = ";", header = !0) %>% 
                     dplyr::select(-dplyr::one_of("V1")) %>% arrange(Response) %>% as.data.frame())

  if(exclude_low_confidence){
    df = df[!grepl("low_confidence", df$Response),]}
  
  log_info("Finished processing.")
  return(df)
}
