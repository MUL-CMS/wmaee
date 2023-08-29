#!/bin/bash

IMAGE="wmaee"
IMAGE_BACKUP_DIR="/media/backup"
SSHFS_REMOTE="fileserver"
SSHFS_REMOTE_DIRECTORY="/share/homes/dgehringer/calculations/wmaee/data/2023"
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
  docker exec -it wmaee bash -c "chmod 600 /root/.ssh/config"
  # mount the SSHFS drive
  docker exec -it wmaee bash -c "sshfs $SSHFS_REMOTE:$SSHFS_REMOTE_DIRECTORY $IMAGE_BACKUP_DIR"
fi
