#' @title Annotation with LLM for clusters or regions.
#' @description Annotation with LLM for clusters or regions.
#' @param markers The output of FindAllMarkers function in Seurat.
#' @param topN Default is 10. The top N genes were used for anntation. 
#' @param api_key The api key of LLM, like deepseek or chat-gpt.
#' @param tissuename The tissue name.
#' @param ann_type Defualt is region. The annotation types include region, celltype.
#' @param web_search Support web search.
#' @param model Defualt is deepseek-reasoner
#' @return Annotation results.
#' @export 

# set deepseek function
run_aica <- function(markers, topN, api_key, tissuename, ann_type = "celltype", web_search = TRUE, model = "deepseek-reasoner") {
    all_markers_topN <- markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC),.by_group = TRUE) %>%
        slice_head(n = topN)
    result <- tapply(all_markers_topN$gene, all_markers_topN$cluster, function(i) paste0(i[1:topN], collapse=','))
    formatted <- sapply(names(result), function(x) paste0(x, ":", result[[x]]))
    final_input <- paste(formatted, collapse="\n")
    url <- "https://api.deepseek.com/v1/chat/completions"
    headers <- add_headers(
        "Authorization" = paste("Bearer", api_key),
        "Content-Type" = "application/json"
    )
    if(ann_type == "region"){
        body <- list(
            model = model,
            messages = list(list(role = "user", content = paste0('Identify region names in ',tissuename,' of spatial transcriptomics using the following markers separately for each\n row. Only provide the region name(like tumor region,immune region, stromal region) in the output message. Do not show numbers before the name.\n ',final_input))),
            web_search = web_search
        )
    }else if(ann_type == "celltype"){
        body <- list(
            model = model,
            messages = list(list(role = "user", content = paste0('Identify cell type names in ',tissuename,' of single cell transcriptomics using the following markers separately for each\n row. Only provide the cell type name(like T cell, B cell, Fibroblast) in the output message. Do not show numbers before the name.\n ',final_input))),
            web_search = web_search
        )
    }
    response <- POST(url, headers, body = toJSON(body, auto_unbox = TRUE))
    ann_ds <- content(response, "parsed")
    ann_clusters <- unlist(str_split(ann_ds$choices[[1]]$message$content,"\n",simplify = TRUE)[1,])
    return(ann_clusters)
}

