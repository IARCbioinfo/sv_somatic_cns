perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1 | sed -e 's/\(.\)/\1 /g' > $2
