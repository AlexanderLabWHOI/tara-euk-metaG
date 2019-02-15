snakemake   \
        --jobs 40 --use-conda --profile tara-metag  \
        --cluster-config cluster.yaml --cluster " echo sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=TARA.{rule}.{wildcards} --mem={cluster.mem} --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"

