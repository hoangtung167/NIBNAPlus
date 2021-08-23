## NIBNAPlus
NIBNA Plus tool is based on NIBNA tool with an addition of passenger genes filtering and alternative community detection algorithms. The original tool can be found at https://github.com/mandarsc/NIBNA.

## Necessary Logistics
To run the program, you need to install 'python3' along with the following libaries 'seaborn' , 'numpy' , 'pandas' , 'matplotlib' , 'scipy' , 'python-louvain' , 'networkx' , 'pygam'. 

## Running the program 
1. Change the 'Data_DIR' and 'OUT_DIR' in 'utils.py', 'network_node_importance.py', and 'passenger_genes_filtering.py' to the appropriate folder in your computer
2. In the folder NIBNA Plus, run the script
```
python3 nibna_cancer_driver_script.py
```
3. The program will ask you about the use of a custom passenger genes dataset. If you provide yoru own dataset, the program will perform analysis on this data. Otherwise, the program will use the default dataset that we pre-analyzed at 'Data/passenger_genes.csv'
```
Do you want to use custom dataset for passenger genes filtering (True/False): 
```





