# AICA: AI-powered Cell Annotation
AICA (AI-powered Cell Annotation) is an R function that leverages large language models (LLMs) to provide automated annotation for single-cell and spatial transcriptomics data. This tool analyzes differentially expressed genes and uses advanced LLM models (like DeepSeek) to intelligently identify cell types or tissue regions.

## Features
🤖 AI-Driven Annotation: Utilize large language models for intelligent biological annotation

🎯 Multiple Annotation Types: Support for single-cell transcriptomics (cell types) and spatial transcriptomics (tissue regions)

⚡ Flexible Configuration: Adjustable topN gene count, choice of different LLM models

🌐 Web Search Integration: Optional web search to enhance annotation accuracy

🔧 Easy Integration: Seamlessly integrates with Seurat analysis pipelines

## Installation

    devtools::install_github("KainMiA/AICA")

## Dependencies
### Required R Packages
dplyr (>=1.1.4)

httr (>=1.4.7)

jsonlite (>=1.8.8)

## API Key
You need to obtain a DeepSeek API key:

1.Visit the DeepSeek Platform

2.Register an account and create an API key

3.Ensure your account has sufficient balance

## Usage
    # 1. Perform differential expression analysis
    library(Seurat)
    all_markers <- FindAllMarkers(seurat_object)
    
    # 2. Run AICA for annotation
    annotations <- run_aica(
      markers = all_markers,
      topN = 10,
      api_key = "your_deepseek_api_key_here",
      tissuename = "Mouse Brain",
      ann_type = "celltype",  # or "region" for spatial transcriptomics
      web_search = FALSE
    )
    
    # 3. View annotation results
    print(annotations)

## Parameter Details
Parameter	Type	Default	Description
markers	data.frame	Required	Output from Seurat's FindAllMarkers, must contain cluster, gene, and avg_log2FC columns
topN	integer	10	Number of top marker genes per cluster used for annotation
api_key	character	Required	DeepSeek API key
tissuename	character	Required	Tissue name (e.g., "Human Liver", "Mouse Brain")
ann_type	character	"region"	Annotation type: "region" (spatial regions) or "celltype" (cell types)
model	character	"deepseek-reasoner"	LLM model name
web_search	logical	FALSE	Whether to enable web search enhancement
