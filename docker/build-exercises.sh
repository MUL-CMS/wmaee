#!/bin/bash

IMAGE="wmaee"
# backup directory in the container
BACKUP_DIRECTORY=/media/backup
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
  docker run -p $PORT:$PORT -v $SSHFS_MOUNT_POINT:$BACKUP_DIRECTORY -d --name $IMAGE "$IMAGE:latest"
fi
