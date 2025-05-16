
<-- Pipeline for design evaluation with alphafold metrics + Rosetta metrics --> 

cd evaluation_pipeline
python ../evaluate_predictions.py in/1747032813539-7091-ae5c995d_job -o {name}_confidence_summary.tsv
python ../cif_to_pdb.py in/1747032813539-7091-ae5c995d_job
## Prerequisites: 
### install Rosetta and copy rosetta.source.release-371 into evaluation_pipeline or other Rosetta folder (modify the path in rosetta scripts and config files accordingly)
### create empty directories .exelogs/err and .exelogs/out
clean.sh
### pending to implement: move all pdbs generated in/1747032813539-7091-ae5c995d_job that passed with SUCCESS (see output .tsv) to a directory in/1747032813539-7091-ae5c995d_job/pdbs
sbatch --array=1-$(ls in/1747032813539-7091-ae5c995d_job/pdbs | wc -l) relax_metrics_array.slurm 
python collect_scores.py output/relax_metrics_20250514_174149 [--pattern "*.sc"]  


