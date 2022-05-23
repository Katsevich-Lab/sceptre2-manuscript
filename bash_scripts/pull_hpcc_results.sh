source ~/.research_config
LOCAL_DIR=$LOCAL_SCEPTRE2_DATA_DIR"results/"
REMOTE_DIR=$REMOTE_SCEPTRE2_DATA_DIR"results/"
rsync -rltvP $REMOTE_DIR $LOCAL_DIR
