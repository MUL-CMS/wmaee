#!/bin/bash

IMAGE="wmaee"
# backup directory in the container
BACKUP_DIRECTORY=/media/wmaee-backup
PORT=8000
CIFS_USERNAME=dgehringer
CIFS_DOMAIN=DMWZONE
CIFS_VOLUME_NAME="wmaee-backup"
CIFS_CREDENTIALS="/home/dominik/.ssh/.smbcredentials"
CIFS_HOSTNAME=193.171.82.188
CIFS_SHARE_NAME="wmaee-backup"
CIFS_PASSWORD=$(grep "password" $CIFS_CREDENTIALS | cut -d'=' -f2)

function volume_exists {
  if [ "$(docker volume ls -f name=$1 | awk '{print $NF}' | grep -E '^'$1'$')" ];
  then
    return 0
  else
    return 1
  fi
}

if [[ "$(docker images -q "$IMAGE:latest" 2> /dev/null)" == "" ]];
then
  echo "[INFO]: Building image \"$IMAGE:latest\""
  docker build -t $IMAGE .
else
  echo "[INFO]: Image \"$IMAGE\" does exist already"
fi

if volume_exists $CIFS_VOLUME_NAME;
then
  echo "[INFO]: Backup docker volume \"$CIFS_VOLUME_NAME\" exists already"
else
  echo "[INFO]: Creating docker volume \"$CIFS_VOLUME_NAME\""
  docker volume create \
    --driver local \
    --opt type=cifs \
    --opt device=//$CIFS_HOSTNAME/$CIFS_SHARE_NAME \
    --opt o=addr=$CIFS_HOSTNAME,username=$CIFS_USERNAME,password=$CIFS_PASSWORD \
    --name $CIFS_VOLUME_NAME
fi


if [ "$( docker container inspect -f '{{.State.Running}}' $IMAGE )" = "true" ]; then
  echo "[INFO]: The container \"$IMAGE\" is already running"
else
  echo "[INFO]: Starting container \"$IMAGE\""
  docker run -p $PORT:$PORT -d --name $IMAGE -v "$CIFS_VOLUME_NAME:$BACKUP_DIRECTORY" "$IMAGE:latest"
fi