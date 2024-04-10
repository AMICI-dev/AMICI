# Using AMICI inside containers

AMICI docker images are regularly published to
https://hub.docker.com/r/dweindl/amici.

## Docker / Podman

The commands below are for use with Docker. When using Podman, just replace
`docker` by `podman`.

### Create an image

In the AMICI base directory run:

```bash
git archive -o container/amici.tar.gz --format=tar.gz HEAD
cd container && docker build -t $USER/amici:latest .
```
Note that this will include files from the last commit, but no uncommitted
changes.

### Pull published image

```bash
docker pull dweindl/amici:latest
```

## Singularity / Apptainer

### Create an image

In the AMICI base directory run:

```bash
# prepare amici files to be copied to the image
# Note that this will include files from the last commit, but no uncommitted
# changes.
git archive -o container/amici.tar.gz --format=tar.gz HEAD
# install spython if necessary
test -x "$(command -v spython)" || pip install spython
# convert Dockerfile to singularity/apptainer definition using spython
(cd container/ && spython recipe Dockerfile amici.def)
# build image
(cd container/ && singularity build amici.sif amici.def)
```

### Pull published image

```bash
singularity pull docker://dweindl/amici:latest
```
