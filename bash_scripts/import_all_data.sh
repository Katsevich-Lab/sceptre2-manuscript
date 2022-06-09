source ~/.research_config

cd $LOCAL_CODE_DIR"/import-schraivogel-2020"
qsub -q short.q -l m_mem_free=16G $LOCAL_CODE_DIR"/import-schraivogel-2020/run_all.sh"

cd $LOCAL_CODE_DIR"/import-papalexi-2021"
qsub -q short.q -l m_mem_free=16G $LOCAL_CODE_DIR"/import-papalexi-2021/run_all.sh"

cd $LOCAL_CODE_DIR"/import-liscovitch-2021"
qsub -q short.q -l m_mem_free=16G $LOCAL_CODE_DIR"/import-liscovitch-2021/run_all.sh"