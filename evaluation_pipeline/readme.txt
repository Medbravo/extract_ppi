
<-- Pipeline for design evaluation with alphafold metrics + Rosetta metrics --> 
rsync -v -r --no-times --no-perms --exclude='*seed*' /mnt/cerebro.corp/medbravo_backend/1747032813539-7091-ae5c995d_job/alphafold/* in/1747032813539-7091-ae5c995d_job
cd evaluation_pipeline
python ../evaluate_predictions.py ../in/1747032813539-7091-ae5c995d_job -o 1747032813539-7091-ae5c995d_job_confidence_summary.tsv --pae_threshold 15
python ../cif_to_pdb.py ../in/1747032813539-7091-ae5c995d_job
## Prerequisites: 
### install Rosetta and copy rosetta.source.release-371 into evaluation_pipeline or other Rosetta folder (modify the path in rosetta scripts and config files accordingly)
### create empty directories .exelogs/err and .exelogs/out
clean.sh
python copy_success_pdbs.py 1747032813539-7091-ae5c995d_job
sbatch --array=1-$(ls 1747032813539-7091-ae5c995d_job_success_pdbs/*.pdb | wc -l) relax_metrics_array.slurm 1747032813539-7091-ae5c995d_job_success_pdbs
python collect_scores.py 1747032813539-7091-ae5c995d_job_success_pdbs/relax_metrics_20250514_174149 [--pattern "*.sc"]  


