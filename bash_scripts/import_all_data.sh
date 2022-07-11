source ~/.research_config

cd $LOCAL_CODE_DIR"/import-schraivogel-2020"
qsub $LOCAL_CODE_DIR"/import-schraivogel-2020/run_all.sh"

cd $LOCAL_CODE_DIR"/import-papalexi-2021"
qsub $LOCAL_CODE_DIR"/import-papalexi-2021/run_all.sh"

cd $LOCAL_CODE_DIR"/import-liscovitch-2021"
qsub $LOCAL_CODE_DIR"/import-liscovitch-2021/run_all.sh"

cd $LOCAL_CODE_DIR"/import-frangieh-2021"
qsub $LOCAL_CODE_DIR"/import-frangieh-2021/run_all.sh"
