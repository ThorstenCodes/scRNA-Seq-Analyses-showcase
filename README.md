# 🔬 Single-Cell Analysis of COVID-19 Lung Tissue (2025 Update)

This repository contains an updated (2025) analysis of the landmark dataset from  
🧪 **"A molecular single-cell lung atlas of lethal COVID-19"**  
by J.C. Melms et al., *Cell* (2021) – [PMID: 33915568](https://pubmed.ncbi.nlm.nih.gov/33915568/).

---

## 📦 What's Inside

- 🧬 Single-cell RNA-seq analysis using **Scanpy** and **scvi-tools**
- 🐳 A **Docker container** to run memory-intensive parts (e.g. dataset integration) on external machines or cloud environments
- 💻 Designed to **offload heavy processing** from local laptops (e.g. MacBook) which might otherwise crash during full dataset integration
- 📘 Written documentation with step-by-step explanations (ideal for those who prefer reading over watching tutorials)

---

## 🚀 Key Features

- Full scRNA-seq pipeline: filtering, normalization, integration, clustering, visualization
- GPU-compatible with `scvi-tools` for efficient model training
- Docker setup for consistent environments across systems
- Easy-to-follow notebooks with comments and rationale

---

## 🧰 Technologies Used

- Python 3.10
- Scanpy
- scvi-tools
- Docker
- matplotlib / seaborn / numpy / pandas
