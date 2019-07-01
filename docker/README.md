# Build docker image

## build our image
```bash
cd .. && docker build -t nrps-galaxy -f docker/Dockerfile .`
```

## push image to Docker Cloud
```bash
docker login
docker tag nrps-galaxy lcouderc/nrps-galaxy
docker push lcouderc/nrps-galaxy
```