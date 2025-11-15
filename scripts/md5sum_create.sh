#!/bin/sh

outfile="md5_all_fqgz.txt"
: > "$outfile"   # empty the output file

# choose md5 utility available on system
if command -v md5sum >/dev/null 2>&1; then
    HASHER="md5sum"
elif command -v md5 >/dev/null 2>&1; then
    HASHER="md5"
else
    echo "No md5 or md5sum command found" >&2
    exit 1
fi

# Safe IFS and find loop (POSIX)
# - IFS= to keep spaces
# - read -r to keep backslashes
# - no -d '' available in POSIX
find . -type f -name '*.fq.gz' -print | sort | \
while IFS= read -r file; do
    "$HASHER" "$file" >> "$outfile"
    echo $file
done
