
FROM nginx:latest

ENV DOMAIN=modelling.unileoben.ac.at

COPY default.conf /etc/nginx/conf.d/default.conf
COPY "$DOMAIN.cert.pem" /etc/ssl/certs
COPY "$DOMAIN.key" /etc/ssl/private
COPY dhparam.pem /etc/ssl/certs/dhparam.pem