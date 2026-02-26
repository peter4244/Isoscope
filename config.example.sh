# Preprocessing Configuration
# Copy this file to config.sh and edit paths for your system:
#   cp config.example.sh config.sh
#
# Only two paths are needed; sorted/compressed files are derived automatically.

# Unindexed GENCODE GTF (the script will create .sorted.gtf.gz and .sorted.gtf.gz.tbi next to it)
GENCODE_GTF="/path/to/gencode.v49.primary_assembly.annotation.chrnamesedited.gtf"

# SQANTI output directory (must contain sqanti3_corrected.gtf)
SQANTI_DIR="/path/to/sqanti_runs/isoseq_sqanti3_filtered"
