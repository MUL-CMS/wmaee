#!/bin/bash

IMAGE="wmaee"
IMAGE_PROXY="nginx-reverse"
# backup directory in the container
JUPYTERHUB_PORT=8000
HOST_PORT=32654

if [[ "$(docker images -q "$IMAGE:latest" 2> /dev/null)" == "" ]];
then
  echo "[INFO]: Building image \"$IMAGE:latest\""
  docker build -t $IMAGE .
else
  echo "[INFO]: Image \"$IMAGE\" does exist already"
fi

if [[ "$(docker images -q "$IMAGE_PROXY:latest" 2> /dev/null)" == "" ]];
then
  cd nginx-reverse
  echo "[INFO]: Building image \"$IMAGE_PROXY:latest\""
  docker build -t $IMAGE_PROXY .
  cd ..
else
  echo "[INFO]: Image \"$IMAGE_PROXY\" does exist already"
fi


if [ "$( docker container inspect -f '{{.State.Running}}' $IMAGE )" = "true" ]; then
  echo "[INFO]: The container \"$IMAGE\" is already running"
else
  echo "[INFO]: Starting container \"$IMAGE\""
  # https://github.com/s3fs-fuse/s3fs-fuse/issues/1046
  # https://github.com/moby/moby/issues/16233
  docker run --security-opt apparmor:unconfined --cap-add SYS_ADMIN --device /dev/fuse -p $HOST_PORT:$JUPYTERHUB_PORT -d --name $IMAGE "$IMAGE:latest"
  # in previous version I forgot to set the permissions in the Dockerfile, just make sure that this is fixed now
fi


if [ "$( docker container inspect -f '{{.State.Running}}' $IMAGE_PROXY )" = "true" ]; then
  echo "[INFO]: The container \"$IMAGE_PROXY\" is already running"
else
  echo "[INFO]: Starting container \"$IMAGE_PROXY\""
  # https://github.com/s3fs-fuse/s3fs-fuse/issues/1046
  # https://github.com/moby/moby/issues/16233
  docker run -p 80:80 -p 443:443 -d --name $IMAGE_PROXY "$IMAGE_PROXY:latest"
  # in previous version I forgot to set the permissions in the Dockerfile, just make sure that this is fixed now
fi
