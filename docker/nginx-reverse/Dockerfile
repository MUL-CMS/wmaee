
FROM nginx:latest

ENV DOMAIN=modelling.unileoben.ac.at

# install GPG to decrypt the certificates
RUN apt-get update --yes && apt-get install --yes --no-install-recommends gpg && apt-get clean

COPY uncrypt.sh .
COPY passwd .
RUN chmod +x uncrypt.sh


COPY "$DOMAIN.cert.pem.gpg" .
RUN bash -c "./uncrypt.sh $DOMAIN.cert.pem.gpg"
RUN mv "$DOMAIN.cert.pem" /etc/ssl/certs
RUN rm "$DOMAIN.cert.pem.gpg"

COPY "$DOMAIN.key.gpg" .
RUN bash -c "./uncrypt.sh $DOMAIN.key.gpg"
RUN mv "$DOMAIN.key" /etc/ssl/private
RUN rm "$DOMAIN.key.gpg"


COPY dhparam.pem.gpg .
RUN bash -c "./uncrypt.sh dhparam.pem.gpg"
RUN mv dhparam.pem /etc/ssl/certs
RUN rm dhparam.pem.gpg

RUN rm uncrypt.sh passwd

RUN apt-get remove --yes gpg && apt-get --yes autoremove

COPY default.conf /etc/nginx/conf.d/default.conf
