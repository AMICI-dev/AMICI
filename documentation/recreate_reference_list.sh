#!/bin/sh
pandoc --filter=pandoc-citeproc bib.md > references.md
