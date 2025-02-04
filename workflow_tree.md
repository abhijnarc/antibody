
```sh
echo "# Workflow Tree

antibody_nextflow/
│── main.nf                 
│── nextflow.config         
│── data/                   
│   ├── [pref].fa
│   ├── 1ngl_clean_pdb
│   ├── ag_res.txt
│   ├── config_template.txt
│── scripts/                
│   ├── igf_bulk.py
│   ├── clean_ab.sh
│   ├── cdr_parse.py
│   ├── ab_res.py
│   ├── ambig.sh
│   ├── unambig.sh
│   ├── generate_config.sh
│── results/               
│── README.md
" > workflow_tree.md
