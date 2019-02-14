snakemake   \
        --jobs 40 --use-conda --profile tara-metag   \
                        --cluster "sbatch --parsable  --partition=compute --job-name=TARA.{rule}.{wildcards} --mem=4gb --time=48:00:00 --ntasks=1 --nodes=1" -np

