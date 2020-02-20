# AMICI with Docker

## Create image

In the AMICI base directory run:

```bash
git archive -o docker/amici.tar.gz --format=tar.gz HEAD
cd docker && docker build -t $USER/amici:latest .
```

## Published images

AMICI docker images are regularly published to
https://hub.docker.com/layers/dweindl/amici/.
