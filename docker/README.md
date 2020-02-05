# AMICI with Docker

## Create image

In your clone of the AMICI repository, run in the base dir:

```bash
git archive -o docker/amici.tar.gz --format=tar.gz HEAD
cd docker && docker build -t $USER/amici:latest .
```

## Use image

TODO: start Jupyter + example notebooks
