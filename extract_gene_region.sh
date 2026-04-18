#!/usr/bin/env bash
#
# DEPRECATED: extract_gene_region.sh has been replaced by extract_longread_region.sh
# (renamed in v4.0). This stub forwards invocations to the new script so existing
# call sites keep working; please update your scripts and docs.
#
# The legacy --merge-by-group flag is a no-op in the new script (which always
# produces both individual/ and merged/ subdirectories) and is stripped from
# the forwarded arguments.

echo "WARNING: extract_gene_region.sh is deprecated; use extract_longread_region.sh." >&2
echo "         Forwarding this invocation to extract_longread_region.sh..." >&2

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Strip --merge-by-group (no longer needed; the new script always merges)
FORWARDED_ARGS=()
for arg in "$@"; do
    if [[ "$arg" != "--merge-by-group" ]]; then
        FORWARDED_ARGS+=("$arg")
    fi
done

exec bash "${SCRIPT_DIR}/extract_longread_region.sh" "${FORWARDED_ARGS[@]}"
