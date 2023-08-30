#!/bin/bash

IMAGE="wmaee"
# backup directory in the container
PORT=8000

if [[ "$(docker images -q "$IMAGE:latest" 2> /dev/null)" == "" ]];
then
  echo "[INFO]: Building image \"$IMAGE:latest\""
  docker build -t $IMAGE .
else
  echo "[INFO]: Image \"$IMAGE\" does exist already"
fi


if [ "$( docker container inspect -f '{{.State.Running}}' $IMAGE )" = "true" ]; then
  echo "[INFO]: The container \"$IMAGE\" is already running"
else
  echo "[INFO]: Starting container \"$IMAGE\""
  # https://github.com/s3fs-fuse/s3fs-fuse/issues/1046
  # https://github.com/moby/moby/issues/16233
  docker run --security-opt apparmor:unconfined --cap-add SYS_ADMIN --device /dev/fuse -p $PORT:$PORT -d --name $IMAGE "$IMAGE:latest"
  # in previous version I forgot to set the permissions in the Dockerfile, just make sure that this is fixed now
fi
