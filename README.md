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
* markers
This should be the direct output from Seurat::FindAllMarkers() or a dataframe with identical structure.

Required columns: cluster, gene, avg_log2FC.

The function will rank genes by avg_log2FC (descending) within each cluster.

topN
We recommend values between 5-20.

Higher values provide more gene context but increase prompt size and API costs.

Lower values may lead to less specific annotations.

api_key
Keep your API key secure. Consider using environment variables for production:

r
Sys.setenv(DEEPSEEK_API_KEY = "your_key_here")
# Then in function call: api_key = Sys.getenv("DEEPSEEK_API_KEY")
ann_type
"region": Best for spatial transcriptomics data (e.g., identifying "tumor core", "immune niche", "stromal region").

"celltype": Best for single-cell RNA-seq data (e.g., identifying "T cells", "Fibroblasts", "Endothelial cells").

model
"deepseek-reasoner": Recommended for best performance (default).

"deepseek-chat": General-purpose model.

Other models may be available; check DeepSeek's documentation for updates.

web_search
Enable only when you need cutting-edge knowledge (e.g., newly discovered cell types).

May be useful for cancer or developmental biology research where nomenclature evolves rapidly.

Note: Increases response time and may incur higher API costs.
