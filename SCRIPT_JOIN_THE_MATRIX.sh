#!/bin/sh

echo "$(cat ./outputs/tmp/indxs_*)" > ./outputs/tmp/INDXS
echo "$(cat ./outputs/tmp/values_*)" > ./outputs/tmp/COEFS
echo "$(cat ./outputs/tmp/Alp_cff_*)" > ./outputs/tmp/ALP
echo "$(cat ./outputs/tmp/Bta_cff_*)" > ./outputs/tmp/BTA

touch ./outputs/tmp/tr_INDXS
touch ./outputs/tmp/tr_COEFS
