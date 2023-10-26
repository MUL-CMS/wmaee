#!/bin/bash
# this is a demo script do not generate a new certificate
# You need a certificate that is signed by an authority. You have to apply for that at ZID via the CCDB app in Lotus notes
# https://confluence.unileoben.ac.at:8443/pages/viewpage.action?pageId=36405321
DOMAIN="modelling.unileoben.ac.at"

openssl req -new -newkey rsa:4096 -nodes -keyout "$DOMAIN.key" -out "$DOMAIN.csr" -subj "/C=AT/ST=Steiermark/L=Leoben/O=Montanuniversität Leoben/OU=Department Materials Science/CN=$DOMAIN" -addext "subjectAltName = DNS:localhost,DNS:$DOMAIN,DNS:www.$DOMAIN"

