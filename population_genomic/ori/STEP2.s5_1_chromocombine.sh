DIR="parts_Ne_est"
IDFILE="idfile.txt"

/home-fn/users/nscc1082/software/chromocombine-0.0.4/chromocombine -o my_res.unnamed -d $DIR
/home-fn/users/nscc1082/software/fs-2.0.7/scripts/chromopainterindivrename.pl $IDFILE my_res.unnamed my_res

