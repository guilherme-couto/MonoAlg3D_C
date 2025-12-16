#!/usr/bin/env bash

MAX_SIZE=$((50 * 1024 * 1024))
removed=0

echo "Checking for large staged files (>50 MB)..."

git diff --cached --name-only --diff-filter=ACM | while read file; do
    [ -f "$file" ] || continue

    size=$(stat -c%s "$file")

    if [ "$size" -gt "$MAX_SIZE" ]; then
        echo "Removing '$file' from commit ($(numfmt --to=iec "$size"))"
        git reset -q HEAD -- "$file"
        removed=1
    fi
done

if [ "$removed" -eq 1 ]; then
    echo
    echo "Large files were excluded from this commit. They remain in your working tree."
fi

exit 0
