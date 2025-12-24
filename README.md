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
dplyr
httr 
jsonlite

## API Key
You need to obtain a DeepSeek API key:

1.Visit the DeepSeek Platform

2.Register an account and create an API key

3.Ensure your account has sufficient balance
