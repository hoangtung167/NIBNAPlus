## NIBNAPlus
NIBNA Plus tool is based on NIBNA tool with an addition of passenger genes filtering and alternative community detection algorithms. The original tool can be found at https://github.com/mandarsc/NIBNA. [1]

## Necessary Logistics
To run the program, you need to install `python3` along with the following libaries `seaborn`, `numpy`, `pandas`, `matplotlib`, `scipy`, `python-louvain`, `networkx`, `pygam`. [1]

## Running the program 
1. Change the `DATA_DIR` and `OUT_DIR` in `utils.py`, `network_node_importance.py`, and `passenger_genes_filtering.py` to the appropriate folder in your computer
2. In the folder NIBNA Plus, run the script
```
python3 nibna_cancer_driver_script.py
```
3. The program will ask you about the use of a custom passenger genes dataset. If you provide yoru own dataset, the program will perform analysis on this data. Otherwise, the program will use the default dataset that we pre-analyzed at 'Data/passenger_genes.csv'
```
Do you want to use custom dataset for passenger genes filtering (True/False): 
```

4. File output will be in your custom output folder in part 1. *"The list of files created by the script are as follows,* </br>
*`critical_nodes.csv` contains list of all predicted cancer drivers.* </br>
*`cancer_node_importance.jpg` contains a plot showing the distribution of node importance scores.* </br>
*`top_k_50_validated_genes.csv` contains top-50 predicted coding cancer drivers. Similarly, the remaining file names with same name convention contain predicted cancer drivers     for different values of threshold.* </br>
*`top_k_validated_genes_weighted.csv` contains the number of predicted coding drivers validated using CGC.* </br>
*`coding_candidate_drivers_mutations.csv` contains list of predicted coding drivers with mutations.* </br>
*`coding_candidate_drivers_no_mutations.csv` contains list of predicted coding drivers without mutations.* </br>
*`noncoding_candidate_drivers.csv` contains list of predicted non-coding drivers.* </br>
*`performance_metrics.csv` contains precision, recall and f1-score of the predicted coding cancer drivers."* </br> (Mandar, Bioinformatics, 2021)

The results are saved in a csv file saved in `Output` directory where each row indicates the number of top-k coding genes found by this approach.
## Citation
[1] Chaudhary, Mandar S., Vu VH Pham, and Thuc D. Le. "NIBNA: A network-based node importance approach for identifying breast cancer drivers." Bioinformatics (2021).


